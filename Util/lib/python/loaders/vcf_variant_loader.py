"""! @brief Annotated VDB VCF Loader """
#!pylint: disable=invalid-name

##
# @package loaders
# @file vcf_variant_loader.py
#
# @brief  Annotated VDB VCF Loader
#
# @section vcf_variant_loader Description
# provides functions & class for managing variant loads from a VCF entry
#
# @section todo_vcf_variant_loader TODO
# - create __verify_current_variant() to make sure _current_variant is set before trying to access
# - extract common elements to parent class and make VEPLoader a child
#
# @section libraries_vcf_variant_loader Libraries/Modules
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
# @section author_vcf_variant_loader Author(s)
# - Created by Emily Greenfest-Allen (fossilfriend) 2022

import json 
import traceback

from copy import deepcopy
from io import StringIO

from GenomicsDBData.Util.utils import xstr, warning, print_dict, to_numeric, deep_update
from GenomicsDBData.Util.list_utils import qw, is_subset, is_equivalent_list
from GenomicsDBData.Util.postgres_dbi import Database, raise_pg_exception

from AnnotatedVDB.Util.parsers import VcfEntryParser, CONSEQUENCE_TYPES
from AnnotatedVDB.Util.variant_annotator import VariantAnnotator
from AnnotatedVDB.Util.loaders import VariantLoader, REQUIRED_COPY_FIELDS
from AnnotatedVDB.Util.database import VARIANT_ID_TYPES

class VCFVariantLoader(VariantLoader):
    """! functions for loading variants from a VCF file """
    
    def __init__(self, datasource, logFileName=None, verbose=False, debug=False):
        """! VCFVariantLoader base class initializer

            @param datasource          datasource description
            @param logFileName         full path to logging file
            @param verbose             flag for verbose output
            @param debug               flag for debug output
            
            @returns                   An instance of the VCFVariantLoader class with initialized counters and copy buffer
        """
        super(VCFVariantLoader, self).__init__(datasource, logFileName, verbose, debug)
        self.log((type(self).__name__, "initialized"), prefix="INFO")
    
    
    def initialize_copy_sql(self, copyFields=None):
        copyFields = REQUIRED_COPY_FIELDS
        copyFields.extend(["ref_snp_id", "is_multi_allelic", "display_attributes", "allele_frequencies"])
        if self._debug:
            self.log(copyFields, prefix="DEBUG")
            
        super(VCFVariantLoader, self).initialize_copy_sql(copyFields=copyFields)
        
    
    def __parse_alt_alleles(self, vcfEntry):
        """! iterate over alleles and calculate values to load
            add to copy buffer
            @param vcfEntry             the VCF entry for the current variant
        """
        
        if self._debug:
            self.log('Entering ' + type(self).__name__ + '.' + 'parse_alt_alleles', prefix="DEBUG")
        
        variant = self._current_variant # just for ease/simplicity
        variantExternalId = variant.ref_snp_id if hasattr(variant, 'ref_snp_id') else None
        
        for alt in variant.alt_alleles:   
            
            if alt == '.':
                self.log(("Skipping variant", self._current_variant.id, "no alt allele (alt = .)"), prefix="WARNING")
                self.increment_counter('skipped')
                continue
            
            self.increment_counter('variant')
            annotator = VariantAnnotator(variant.ref_allele, alt, variant.chromosome, variant.position)       
            
            if not self.is_dbsnp() and self.skip_existing():
                if self.is_duplicate(annotator.get_metaseq_id(), 'METASEQ'):
                    self.increment_counter('duplicates')
                    continue
                
            recordPK = self._pk_generator.generate_primary_key(annotator.get_metaseq_id(), variantExternalId)      
                
            if self.is_adsp(): # check for duplicates and update is_adsp_variant flag
                if self.is_duplicate(recordPK, 'PRIMARY_KEY'):
                    self.__add_adsp_update_statement(recordPK, self._current_variant.chromosome)
                    self.increment_counter('update')
                    continue
            
            normRef, normAlt = annotator.get_normalized_alleles()
            binIndex = self._bin_indexer.find_bin_index(variant.chromosome, variant.position,
                                                        vcfEntry.infer_variant_end_location(alt, normRef))

            alleleFreq = vcfEntry.get_frequencies(alt)
        
            copyValues = ['chr' + xstr(variant.chromosome),
                        recordPK,
                        xstr(variant.position),
                        annotator.get_metaseq_id(),
                        binIndex,
                        xstr(self._alg_invocation_id),
                        xstr(variant.ref_snp_id, nullStr='NULL'),                        
                        xstr(variant.is_multi_allelic, falseAsNull=True, nullStr='NULL'),
                        xstr(annotator.get_display_attributes(variant.rs_position), nullStr='NULL'),
                        xstr(alleleFreq, nullStr='NULL')
                        ]      
            
            if self.is_adsp():
                copyValues.append(xstr(True)) # is_adsp_variant
                
            if self._debug:
                self.log("Display Attributes " + print_dict(annotator.get_display_attributes(variant.rs_position)), prefix="DEBUG")
                self.log(copyValues, prefix="DEBUG")
                
            self.add_copy_str('#'.join(copyValues))
        
        
    def parse_variant(self, line):
        """! parse & load single line from file 
        @param line             the VCF entry in string form
        @returns copy string for db load 
        """

        if self._debug:
            self.log('Entering ' + type(self).__name__ + '.' + 'parse_result', prefix="DEBUG")
            self.log(('Resume?', self.resume_load()), prefix="DEBUG")
            
        if self.resume_load() is False and self._resume_after_variant is None:
            err = ValueError('Must set VariantLoader resume_afer_variant if resuming load')
            self.log(str(err), prefix="ERROR")
            raise err
        
        try:
            self.increment_counter('line')
            entry = VcfEntryParser(line)
            
            if not self.resume_load():
                self._update_resume_status(entry.get('id'))
                return None
            
            # otherwise proceed with the parsing       
            entry.update_chromosome(self._chromosome_map)
            
            # extract identifying variant info for frequent reference & logging in the loader calling script
            self._current_variant = entry.get_variant(dbSNP=self.is_dbsnp(), namespace=True)
            
            if self.is_dbsnp() and self.skip_existing():
                if self.is_duplicate(self._current_variant.ref_snp_id, 
                                     'REFSNP', 'chr' + self._current_variant.chromosome):
                    self.increment_counter('duplicates')
                    return None
            
            # iterate over alleles
            self.__parse_alt_alleles(entry)
            
        except Exception as err:
            self.log(str(err), prefix="ERROR")
            raise err