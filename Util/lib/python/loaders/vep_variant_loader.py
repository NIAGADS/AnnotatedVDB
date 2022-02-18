"""! @brief Annotated VDB Loader """
#!pylint: disable=invalid-name
##
# @file vep_variant_loader.py
#
# @brief  Annotated VDB Loader
#
# @section vep_variant_loader Description
# provides functions & class for managing variant loads
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
import traceback

from copy import deepcopy
from io import StringIO

from GenomicsDBData.Util.utils import xstr, warning, print_dict, to_numeric, deep_update
from GenomicsDBData.Util.list_utils import qw, is_subset, is_equivalent_list
from GenomicsDBData.Util.postgres_dbi import Database, raise_pg_exception

from AnnotatedVDB.Util.parsers import VcfEntryParser, VepJsonParser
from AnnotatedVDB.Util.variant_annotator import VariantAnnotator
from AnnotatedVDB.Util.loaders import VariantLoader

class VEPVariantLoader(VariantLoader):
    """! functions for loading variants from VEP result """
    
    def __init__(self, datasource, logFileName=None, verbose=False, debug=False):
        """! VEPVariantLoader base class initializer

            @param copySql             SQL for copy statement
            @param chromosomeMap       chromosome map
            @param datasource          datasource description
            @param logFileName         full path to logging file
            @param verbose             flag for verbose output
            @param debug               flag for debug output
            
            @returns                   An instance of the VepVariantLoader class with initialized counters and copy buffer
        """
        super(VEPVariantLoader, self).__init__(datasource, logFileName, verbose, debug)
        self.__vep_parser = None
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
        
    
    def __parse_vcf_frequencies(self, allele, vcfGMAFs):
        """! retrieve allele frequencies reported in INFO FREQ field
        @param allele          allele to match (not normalized)
        @param vcfGMAF         global minor allele frequencies from the VCF FREQ info field
        @returns               dict of source:frequency for the allele
        """
        
        if vcfGMAFs is None:
            return None
        
        zeroValues = ['.', '0']

        # FREQ=GnomAD:0.9986,0.001353|Korea1K:0.9814,0.01861|dbGaP_PopFreq:0.9994,0.0005901
        # altIndex needs to be incremented as first value is for the ref allele)
        altIndex = self._current_variant.alt_alleles.index(allele) + 1
        populationFrequencies = {pop.split(':')[0]:pop.split(':')[1] for pop in vcfGMAFs.split('|')}
        vcfFreqs = {pop: {'gmaf': to_numeric(freq.split(',')[altIndex])} \
                    for pop, freq in populationFrequencies.items() \
                    if freq.split(',')[altIndex] not in zeroValues}
    
        return None if len(vcfFreqs) == 0 else vcfFreqs


    def __add_vcf_frequencies(self, alleleFreqs, vcfFreqs, allele):
        """! add VCF GMAF from the VCF FREQ INFO field to allele frequencies 

            @param alleleFreqs         allele frequencies from VEP result
            @param vcfFreqs            FREQ string from VCF INFO field
            @param allele             _description_
            @returns                 _description_
        """
        
        alleleGMAFs = self.__parse_vcf_frequencies(allele, vcfFreqs)
        if alleleFreqs is None:
            return alleleGMAFs
        if alleleGMAFs is None:
            return alleleFreqs
        
        # update the alleleFreq dict, taking into account nestedness
        return deep_update(alleleFreqs, alleleGMAFs)
    
        
    def __parse_alt_alleles(self, vcfEntry, resultJson):
        """! iterate over alleles and calculate values to load
             add to copy buffer
            @param vcfEntry             the VCF entry for the current variant
            @param resultJSON           the VEP result / to be added to copyStr
        """
        # extract result frequencies
        if self._debug:
            self.log('Entering ' + type(self).__name__ + '.' + 'parse_alt_alleles', prefix="DEBUG")
            
        frequencies = self.__get_result_frequencies()
        variant = self._current_variant # just for ease/simplicity
        variantExternalId = variant.ref_snp_id if hasattr(variant, 'ref_snp_id') else None
        
        for alt in variant.alt_alleles:
            self.increment_counter('variant')
            annotator = VariantAnnotator(variant.ref_allele, alt, variant.chromosome, variant.position)       
            normRef, normAlt = annotator.get_normalized_alleles()
            binIndex = self._bin_indexer.find_bin_index(variant.chromosome, variant.position,
                                                        vcfEntry.infer_variant_end_location(alt, normRef))
            
            # NOTE: VEP uses left normalized alleles to indicate freq allele
            # so need to use normalized alleles to match freq allele / consequence allele
            # to the correct variant alt allele
            alleleFreq = None if frequencies is None \
                else self.__get_allele_frequencies(normAlt, frequencies)    
            alleleFreq = self.__add_vcf_frequencies(alleleFreq, vcfEntry.get_info('FREQ'), alt)
                        
            msConseq = self.__vep_parser.get_most_severe_consequence(normAlt)

            recordPK = self._pk_generator.generate_primary_key(annotator.get_metaseq_id(), variantExternalId)

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
                        xstr(self.__vep_parser.get_allele_consequences(normAlt), nullStr='NULL'),
                        xstr(resultJson),
                        xstr(self._alg_invocation_id)
                        ]           
       
            if self._debug:
                self.log("Display Attributes " + print_dict(annotator.get_display_attributes(variant.rs_position)), prefix="DEBUG")
                self.log(copyValues, prefix="DEBUG")
            self.add_copy_str('#'.join(copyValues))

    
    def parse_result(self, result):
        """! parse & load single line from file 
        @params result             the VEP result in string form
        @params checkSkip          make checks to see if line should be skipped (after resume)
        @returns copy string for db load 
        """
        if self._debug:
            self.log('Entering ' + type(self).__name__ + '.' + 'parse_result', prefix="DEBUG")
            self.log(('Resume?', self.resume_load()), prefix="DEBUG")
            
        if self.resume_load() is False and self._resume_after_variant is None:
            err = ValueError('Must set VariantLoader result_afer_variant if resuming load')
            self.log(str(err), prefix="ERROR")
            raise err
        
        try:
            self.increment_counter('line')
            resultJson = json.loads(result)
            self.__vep_parser.set_annotation(deepcopy(resultJson))
 
            entry = VcfEntryParser(self.__vep_parser.get('input'))
                
            if not self.resume_load():
                self._update_resume_status(entry.get('id'))
                return None
            
            # otherwise proceed with the parsing            
            entry.update_chromosome(self._chromosome_map)
            
            # there are json formatting issues w/the input str
            # so replace w/the parsed entry; which is now a dict
            resultJson['input'] = entry.get_entry()
            
            # extract identifying variant info for frequent reference
            self._current_variant = entry.get_variant(dbSNP=self.is_dbsnp(), namespace=True)
        
            # rank consequences
            self.__vep_parser.adsp_rank_and_sort_consequences()
    
            # iterate over alleles
            self.__parse_alt_alleles(entry, resultJson)
            
        except Exception as err:
            self.log(str(err), prefix="ERROR")
            raise err
            
    
    