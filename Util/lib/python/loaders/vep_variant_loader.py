"""! @brief Annotated VDB VEP Loader """
#!pylint: disable=invalid-name
##
# @package loaders
# @file vep_variant_loader.py
#
# @brief  Annotated VDB VEP Loader
#
# @section vep_variant_loader Description
# provides functions & class for managing variant loads from VEP result 
#
# @section todo_vep_variant_loader TODO
# - create __verify_current_variant() to make sure _current_variant is set before trying to access
# - extract common elements to parent class and make VEPLoader a child
#
# @section libraries_vep_variant_loader Libraries/Modules
# - json - parse/serialize JSON
# - copy - for deep copy functionality 
# - [GenomicsDBData.Util.utils](https://github.com/NIAGADS/GenomicsDBData/blob/master/Util/lib/python/utils.py)
#   + provides variety of wrappers for standard file, string, and logging operations
# - [GenomicsDBData.Util.list_utils](https://github.com/NIAGADS/GenomicsDBData/blob/master/Util/lib/python/list_utils.py)
#   + provides variety list and set operations
# - [GenomicsDBData.Util.postgres_dbi](https://github.com/NIAGADS/GenomicsDBData/blob/master/Util/lib/python/postgres_dbi.py)
#   + wrapper for connecting to database & handling PG errors
# - [AnnotatedVDB.Util.parsers](https://github.com/NIAGADS/AnnotatedVDB/blob/master/Util/lib/python/parsers)
#   + parsers for VCF and VEP JSON
# - [AnnotatedVDB.Util.loaders](https://github.com/NIAGADS/AnnotatedVDB/blob/master/Util/lib/python/loaders)
#   + parent class 
# - [AnnotatedVDB.Util.variant_annotator](https://github.com/NIAGADS/AnnotatedVDB/blob/master/Util/lib/python/variant_annotator.py)
#   + generate standard variant annotations from position and allele information
# @section author_vep_variant_loader Author(s)
# - Created by Emily Greenfest-Allen (fossilfriend) 2022

import json 

from logging import StreamHandler

from copy import deepcopy

from GenomicsDBData.Util.utils import xstr, print_dict, deep_update
from GenomicsDBData.Util.list_utils import qw

from AnnotatedVDB.Util.parsers import VcfEntryParser, VepJsonParser, CONSEQUENCE_TYPES
from AnnotatedVDB.Util.variant_annotator import VariantAnnotator
from AnnotatedVDB.Util.loaders import VariantLoader, JSONB_UPDATE_FIELDS

NBSP = " " # for multi-line sql

COPY_FIELDS = qw('allele_frequencies adsp_most_severe_consequence adsp_ranked_consequences vep_output row_algorithm_id', returnTuple=False)

class VEPVariantLoader(VariantLoader):
    """! functions for loading variants from VEP result """
    
    def __init__(self, datasource, verbose=False, debug=False):
        """! VEPVariantLoader base class initializer

            @param datasource          datasource description
            @param logFileName         full path to logging file
            @param verbose             flag for verbose output
            @param debug               flag for debug output
            
            @returns                   An instance of the VepVariantLoader class with initialized counters and copy buffer
        """
        super(VEPVariantLoader, self).__init__(datasource, verbose, debug)
        self.__vep_parser = None
        self.logger.info(type(self).__name__ + " initialized")
        


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


    def __add_vcf_frequencies(self, alleleFreqs, vcfEntry, allele):
        """! add VCF GMAF from the VCF FREQ INFO field to allele frequencies 

            @param alleleFreqs         allele frequencies from VEP result
            @param vcfEntry            VCF entry
            @param allele              alt allele to match
            @returns                  updated frequency dictionary
        """
        
        alleleGMAFs = vcfEntry.get_frequencies(allele)
        if alleleFreqs is None:
            return alleleGMAFs
        if alleleGMAFs is None:
            return alleleFreqs
        
        # update the alleleFreq dict, taking into account nestedness
        return deep_update(alleleFreqs, alleleGMAFs)
        
    
    def __get_variant_primary_key(self, metaseqId):
        """ fetch variant primary key """
        result = self.is_duplicate(metaseqId, returnMatch=True)
        if result is None:
            raise KeyError("No record for variant %s found in DB.  VEP Variant Loader does updates only.", metaseqId)
        return result[0]['primary_key']
        
        
    def __parse_alt_alleles(self, vcfEntry):
        """! iterate over alleles and calculate values to load
            add to copy buffer
            @param vcfEntry             the VCF entry for the current variant
        """
        
        # extract result frequencies    
        frequencies = self.__get_result_frequencies()
        variant = self._current_variant # just for ease/simplicity
            
        # clean result -- remove elements in the result JSON that are extracted 
        # and saved to other fields
        # does a deep copy so won't affect allele-specific extractions, so do just once
        # as this will remove the info for all alleles (frequencies, consequences)
        cleanResult = self.__clean_result() 
        
        for alt in variant.alt_alleles:
            self.increment_counter('variant')
            update = False
            
            annotator = VariantAnnotator(variant.ref_allele, alt, variant.chromosome, variant.position)    
            recordPK = self.__get_variant_primary_key(annotator.get_metaseq_id())

            if self.has_attribute('vep_output', recordPK, returnVal=False):
                if self.skip_existing(): # otherwise will update existing record
                    self.increment_counter('duplicates')
                    if self._log_skips:
                        self.logger.warning("Existing data found for: %s | %s; SKIPPING", annotator.get_metaseq_id(), variant)
                    continue
                if self._log_skips:
                    self.logger.warning("Existing data found for: %s | %s; UPDATING", annotator.get_metaseq_id(), variant)
                    
            # NOTE: VEP uses left normalized alleles to indicate freq allele, with - for full deletions
            # so need to use normalized alleles to match freq allele / consequence allele
            # to the correct variant alt allele
            normRef, normAlt = annotator.get_normalized_alleles(snvDivMinus=True)
            
            alleleFreq = None if frequencies is None \
                else self.__get_allele_frequencies(normAlt, frequencies)    
            alleleFreq = self.__add_vcf_frequencies(alleleFreq, vcfEntry, alt)
                        
            msConseq = self.__vep_parser.get_most_severe_consequence(normAlt)
    
            updateValues = [
                recordPK,
                'chr' + xstr(variant.chromosome),
                xstr(alleleFreq, nullStr='{}'),
                xstr(msConseq, nullStr='{}'),
                xstr(self.__vep_parser.get_allele_consequences(normAlt), nullStr='{}'),
                xstr(cleanResult, nullStr='{}'), # cleaned result JSON
                int(self._alg_invocation_id)
            ]      
            
            if self.is_adsp():
                updateValues.append(xstr(True)) # is_adsp_variant
            
            if self._debug:
                self.logger.debug("Update Values: %s", updateValues)

            self.increment_counter('update')
            self._update_buffer.append(updateValues)

    
    def initialize_update_sql(self):
        """ generate update sql """
        # chromosome record_primary_key allele_frequencies 
        # adsp_most_severe_consequence adsp_ranked_consequences vep_output row_algorithm_id
        fields = COPY_FIELDS
        if self.is_adsp():
            fields.append('is_adsp_variant')
            
        updateString = ""
        for f in fields:
            if f in JSONB_UPDATE_FIELDS:
                updateString += ' '.join((f,  '=', 'jsonb_merge(COALESCE(v.' + f, ",'{}'::jsonb),", 'd.' + f + '::jsonb)')) + ', ' # e.g. v.gwas_flags = jsonb_merge(v.gwas_flags, d.gwas_flags::jsonb)
            elif f == 'bin_index':
                updateString += ' '.join((f,  '=', 'd.' + f + '::ltree')) + ', '
            else:
                updateString += ' '.join((f,  '=', 'd.' + f)) + ', '
        
        updateString = updateString.rstrip(', ') + NBSP
        
        if self.is_adsp():
            fields.append('is_adsp_variant')
            updateString += ", v.is_adsp_variant = d.is_adsp_variant" + NBSP

        fields = ['record_primary_key', 'chromosome'] + fields
        schema = "AnnotatedVDB"
        sql = "UPDATE " + schema + ".Variant v SET" + NBSP + updateString \
            + "FROM (VALUES %s) AS d(" + ','.join(fields) + ")" + NBSP \
            + "WHERE v.chromosome = d.chromosome" + NBSP \
            + "AND v.record_primary_key = d.record_primary_key"
        
        self.set_update_sql(sql)
        if self._debug:        
            self.logger.debug("Update SQL: " + self._update_sql)
            
    
    def parse_variant(self, line):
        """! parse & load single line from file 
        @param line             the VEP result in string form
        @returns copy string for db load 
        """
        if self._variant_validator is None:
            raise UnboundLocalError("Validator has not been initialized, please execute the initialize_variant_validator method before parsing variants")
                    
        if self.resume_load() is False and self._resume_after_variant is None:
            raise ValueError('Must set VariantLoader resume_afer_variant if resuming load')
            
        newConseqCount = self.__vep_parser.get_consequence_parser().get_new_conseq_count()
                
        try:
            self.increment_counter('line')
            resultJson = json.loads(line)
            self.__vep_parser.set_annotation(deepcopy(resultJson))

            entry = VcfEntryParser(self.__vep_parser.get('input'), identityOnly=self.is_adsp())
                
            if not self.resume_load():
                self._update_resume_status(entry.get('id'))
                return None
            
            # otherwise proceed with the parsing       
            entry.update_chromosome(self._chromosome_map)
            
            # there are json formatting issues w/the input str
            # so replace w/the parsed entry; which is now a dict
            self.__vep_parser.set('input', entry.get_entry())
            
            # extract identifying variant info for frequent reference
            self._current_variant = entry.get_variant(dbSNP=self.is_dbsnp(), namespace=True)
        
            # rank consequences
            self.__vep_parser.adsp_rank_and_sort_consequences()
            
            # log any new consequences
            if self.__vep_parser.get_consequence_parser().get_new_conseq_count() > newConseqCount:
                newConseqCount = self.__vep_parser.get_consequence_parser().get_new_conseq_count()
                self.logger.warning("New consequence added for %s - %s",
                    self._current_variant.id, 
                    self.__vep_parser.get_consequence_parser().get_added_consequences(mostRecent=True)),
                    
            # iterate over alleles
            self.__parse_alt_alleles(entry)
            
        except Exception as err:            
            raise err
            

    