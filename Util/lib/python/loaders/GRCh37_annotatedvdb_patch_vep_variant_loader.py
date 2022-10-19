"""! @brief Annotated VDB VEP Loader """
#!pylint: disable=invalid-name
##
# @package loaders
# temp loader to migrate GRCh37 external AnnotatedVDB into GUS database
# ##

import json 
import traceback

from types import SimpleNamespace
from copy import deepcopy
from io import StringIO

from GenomicsDBData.Util.utils import xstr, warning, print_dict, to_numeric, deep_update, die
from GenomicsDBData.Util.list_utils import qw, is_subset, is_equivalent_list
from GenomicsDBData.Util.postgres_dbi import Database, raise_pg_exception

from AnnotatedVDB.Util.parsers import VcfEntryParser, VepJsonParser, CONSEQUENCE_TYPES
from AnnotatedVDB.Util.variant_annotator import VariantAnnotator
from AnnotatedVDB.Util.loaders import VEPVariantLoader
from AnnotatedVDB.Util.database import VARIANT_ID_TYPES

class PatchVEPVariantLoader(VEPVariantLoader):
    """! functions for loading variants from VEP result """
    
    def __init__(self, datasource, logFileName=None, verbose=False, debug=False):
        """! VEPVariantLoader base class initializer

            @param datasource          datasource description
            @param logFileName         full path to logging file
            @param verbose             flag for verbose output
            @param debug               flag for debug output
            
            @returns                   An instance of the VepVariantLoader class with initialized counters and copy buffer
        """
        super(VEPVariantLoader, self).__init__(datasource, logFileName, verbose, debug)
        self.__vep_parser = None
        self._initialize_counters(['record'])
        self.log((type(self).__name__, "initialized"), prefix="INFO")
        


    def initialize_vep_parser(self, rankingFile, rankConsequencesOnLoad, verbose):
        """! initialize VEP Parser

            @param rankingFile                file containing consequence ranks
            @param rankConsequencesOnLoad     rank consequences on load?
            @param verbose                  
        """
        self.__vep_parser = VepJsonParser(rankingFile, rankConsequencesOnLoad=rankConsequencesOnLoad, verbose=verbose)
        
        
    def vep_parser(self):
        """! @returns VEP Parser """
        return self.__vep_parser
    

    def __get_result_frequencies(self):
        """! get allele frequencies from the VEP result, catching null objects
            @returns                frequency JSON
        """
        matchingVariantId = self._current_variant.ref_snp_id if self.is_dbsnp() else None
        frequencies = self.__vep_parser.get_frequencies(matchingVariantId)
        if frequencies is None:
            return None
        return frequencies['values']     
    
    
    def __get_allele_frequencies(self, allele, frequencies):
        """! given an allele and frequency dict, retrieve and return allele-specific frequencies
            @param allele          allele to be matched
            @param frequencies     result frequencies as JSON/dict
            
            @returns 
        """
        if frequencies is None:
            return None

        if allele in frequencies:
            return frequencies[allele]

        return None
    
    
    def __clean_result(self):
        """! remove attributes from VEP result that have already been extract
        to reduce loading overhead 
            @returns a cleaned JSON result
        """
        result = self.__vep_parser.get_annotation(deepCopy=True)
        result.pop('colocated_variants', None) # allele frequencies
        for ctype in CONSEQUENCE_TYPES:
            ctypeKey = ctype + '_consequences'
            result.pop(ctypeKey, None)
            
        return result


    
    def __infer_variant_end_location(self, normRef):
        """! infer span of indels/deletions for a 
        specific alternative allele, modeled off 
        GUS Perl VariantAnnotator & dbSNP normalization conventions

        """
        variant = self._current_variant
        ref = variant.ref_allele
        alt = variant.alt_alleles
        # CT/CTTC - insTC

        rLength = len(ref)
        aLength = len(alt)

        position = int(variant.position)
        rsPosition = variant.rs_position
        if rsPosition is None:
            rsPosition = position
        else:
            rsPosition = int(rsPosition)

        if rLength == 1 and aLength == 1: # SNV
            return position

        if rLength == aLength: # MNV
            if ref == alt[::-1]: #inversion
                return position + rLength - 1

            # substitution
            return position + len(normRef) - 1

        if rLength > aLength: # deletions
            if len(alt) > 1: # indel
                if len(normRef) == 0: # was normalized; adjust
                    return rsPosition + len(ref) - 2
                return rsPosition + len(ref) - 1
            else: # straight up deletion
                return rsPosition + len(normRef) -  1

        if rLength < aLength: # insertion
            return rsPosition + 1


    def __clean_qc_info(self, infoStr):
        ''' remove quotes and escaped quotes / the json encoding & dumping is confusing them 
        easier to remove than resolve '''

        cleanStr = infoStr
        cleanStr = cleanStr.replace('\\"', '')
        cleanStr = cleanStr.replace("\\'", '')
        cleanStr = cleanStr.replace('"', '')
        cleanStr = cleanStr.replace("'", '')

        return cleanStr
        

    def __update_adsp_qc(self, qcJson):
        """ parse the ADSP QC result and adjust to match format of GRCh38 result """
        if qcJson is not None:
            for release, result in qcJson.items():
                if len(result) == 1: # just filter status
                    qcJson[release] = {'filter': result['FILTER_STATUS']}
                else:
                    result['INFO_BAYLOR'] = self.__clean_qc_info(result['INFO_BAYLOR'])
                    result['INFO_BROAD'] = self.__clean_qc_info(result['INFO_BROAD'])
                    qcJson[release] = { 'info': result, 
                                        'filter': result['FILTER_STATUS'],
                                        'qual': result['FILTER']}

        return qcJson


    def __parse_alt_alleles(self, record):
        """! iterate over alleles and calculate values to load
            add to copy buffer
            @param vcfEntry             the VCF entry for the current variant
        """
        
        hasVepResult = record['vep_output'] is not None

        # extract result frequencies    
        frequencies = self.__get_result_frequencies() if hasVepResult else None
        variant = self._current_variant # just for ease/simplicity
        variantExternalId = variant.ref_snp_id if hasattr(variant, 'ref_snp_id') else None

        # clean result -- remove elements in the result JSON that are extracted 
        # and saved to other fields
        # does a deep copy so won't affect allele-specific extractions, so do just once
        # as this will remove the info for all alleles (frequencies, consequences)
        cleanResult = self.__clean_result() if hasVepResult else None
        
        for alt in variant.alt_alleles:
            self.increment_counter('variant')
            annotator = VariantAnnotator(variant.ref_allele, alt, variant.chromosome, variant.position)       
                                    
            normRef, normAlt = annotator.get_normalized_alleles()
            binIndex = self._bin_indexer.find_bin_index(variant.chromosome, variant.position,
                                                        self.__infer_variant_end_location(normRef))

            # NOTE: VEP uses left normalized alleles to indicate freq allele, with - for full deletions
            # so need to use normalized alleles to match freq allele / consequence allele
            # to the correct variant alt allele
            normRef, normAlt = annotator.get_normalized_alleles(snvDivMinus=True)
            alleleFreq = None if frequencies is None \
                else self.__get_allele_frequencies(normAlt, frequencies)    
                        
            msConseq = self.__vep_parser.get_most_severe_consequence(normAlt) if hasVepResult else None

            recordPK = self._pk_generator.generate_primary_key(annotator.get_metaseq_id(), variantExternalId)      
    
            # DEFAULT_COPY_FIELDS = qw('chromosome record_primary_key position is_multi_allelic bin_index ref_snp_id 
            # metaseq_id display_attributes allele_frequencies adsp_most_severe_consequence 
            # adsp_ranked_consequences vep_output row_algorithm_id', returnTuple=False)
            # ['is_adsp_variant', 'adsp_qc', 'cadd_scores', 'other_annotation']

            copyValues = ['chr' + xstr(variant.chromosome), 
                        recordPK,
                        xstr(variant.position),
                        xstr(variant.is_multi_allelic, falseAsNull=True, nullStr='NULL'),
                        binIndex,
                        xstr(variant.ref_snp_id, nullStr='NULL'),
                        annotator.get_metaseq_id(),
                        xstr(annotator.get_display_attributes(variant.rs_position), nullStr='NULL'),
                        xstr(alleleFreq, nullStr='NULL'),
                        xstr(msConseq, nullStr='NULL'),
                        xstr(self.__vep_parser.get_allele_consequences(normAlt), nullStr='NULL') if hasVepResult else 'NULL',
                        xstr(cleanResult, nullStr='NULL'), # cleaned result JSON
                        xstr(record['is_adsp_variant'], nullStr='NULL'),
                        xstr(self.__update_adsp_qc(record['adsp_qc']), nullStr='NULL'),
                        xstr(record['cadd_scores'], nullStr='NULL'),
                        xstr(record['other_annotation'], nullStr='NULL'),
                        xstr(self._alg_invocation_id)
                        ]      
 
            if self._debug:
                self.log("Display Attributes - " + print_dict(annotator.get_display_attributes(variant.rs_position)), prefix="DEBUG")
                self.log(("Copy Values - ", copyValues), prefix="DEBUG")
            self.add_copy_str('#'.join(copyValues))



    def __set_variant(self, record, namespace=True):
        """! extract basic variant attributes from a Vthe parsed CF entry 
        @param namespace              if True return SimpleNamespace, else return dict
        @returns attributes as a simple namespace so they can be accessed in dot notation"""

        attributes = {}
        chrom = xstr(record['chromosome'])
        if chrom == 'MT':
            chrom = 'M'
        id = record['metaseq_id']

        hasVepResult = record['vep_output'] is not None
        info = record['vep_output']['input']['info'] if hasVepResult else None 
 
        isAdspVariant = record['is_adsp_variant']
        
        ref = record['ref_allele'] 
        alt = record['alt_allele']
            
        variant =  {
            'id' : id,
            'ref_snp_id' : record['ref_snp_id'],
            'ref_allele' : 'N' if not isAdspVariant and ref == '?' else ref,
            'alt_alleles' : ['N' if not isAdspVariant and alt == '?' else alt],
            'is_multi_allelic' : record['is_multi_allelic'],
            'chromosome' : xstr(chrom).replace('chr', ''),
            'position' : int(record['position']),
            'rs_position' : info['RSPOS'] if info is not None and 'RSPOS' in info else None
        }

        if self._debug:
            self.log("Variant " + print_dict(variant), prefix="DEBUG")
        
        return SimpleNamespace(**variant) if namespace else variant

    
    def parse_variant(self, record):
        """! parse & load record from external database
        @param record   record object from DB query
        @returns copy string for db load 
        """
                
        newConseqCount = self.__vep_parser.get_consequence_parser().get_new_conseq_count()
        
        try:
            self.increment_counter('record')
            vepResult = record['vep_output']
            if self._debug:
                self.log("VEP Output - " + json.dumps(vepResult), prefix="DEBUG")
            
            if vepResult is not None:
                self.__vep_parser.set_annotation(deepcopy(vepResult))
            else:
                self.__vep_parser.set_annotation(None)

            # extract identifying variant info for frequent reference
            self._current_variant = self.__set_variant(record, namespace=True)
    
            # rank consequences
            if vepResult is not None:
                self.__vep_parser.adsp_rank_and_sort_consequences()
            
            # log any new consequences
            if self.__vep_parser.get_consequence_parser().get_new_conseq_count() > newConseqCount:
                newConseqCount = self.__vep_parser.get_consequence_parser().get_new_conseq_count()
                self.log(("New consequence added for", self._current_variant.id, "-",
                    self.__vep_parser.get_consequence_parser().get_added_consequences(mostRecent=True)),
                    prefix="WARNING")
    
            # copying over so mutli-allelic variants already split; so there will be only one allele  per variant
            self.__parse_alt_alleles(record)
            
        except Exception as err:
            self.log(str(err), prefix="ERROR")
            raise err
            
    
    
