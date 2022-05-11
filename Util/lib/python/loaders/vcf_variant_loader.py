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
from AnnotatedVDB.Util.loaders import VariantLoader, REQUIRED_COPY_FIELDS, JSONB_UPDATE_FIELDS, BOOLEAN_FIELDS
from AnnotatedVDB.Util.database import VARIANT_ID_TYPES

NBSP = " " # for multi-line sql
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
        
        self.__update_value_generator = None
        self.__update_existing = False
        self.__update_fields = None
        self.__vcf_header_fields = None
        self.__chromosome = None
        super(VCFVariantLoader, self).__init__(datasource, logFileName, verbose, debug)
        self.log((type(self).__name__, "initialized"), prefix="INFO")
 
 
    def set_chromosome(self, chrom):
        self.__chromosome = chrom
    

    def chromosome(self):
        return self.__chromosome
 
 
    def set_vcf_header_fields(self, fields):
        self.__vcf_header_fields = fields
        
        
    def vcf_header_fields(self):
        return self.__vcf_header_fields
    
 
    def set_update_existing(self, updateExisting):
        self.__update_existing = updateExisting
        
        
    def update_existing(self):
        return self.__update_existing
    
    
    def initialize_copy_sql(self, copyFields=None):
        fields = REQUIRED_COPY_FIELDS  # "chromosome", "record_primary_key", "position", "metaseq_id", "bin_index", "row_algorithm_id"]   
        fields.extend(["ref_snp_id", "is_multi_allelic", "display_attributes", "allele_frequencies"])
        fields.extend(copyFields)
        
        if self._debug:
            self.log(("COPY fields", fields), prefix="DEBUG")
            
        super(VCFVariantLoader, self).initialize_copy_sql(copyFields=fields)
        
        
    def set_update_fields(self, fields):
        self.__update_fields = deepcopy(fields)
        
    
    def set_update_value_generator(self, func):
        self.__update_value_generator = func
        
    
    def generate_update_values(self, entry, flag=None):
        return self.__update_value_generator(self, entry, flag)
    
    
    def build_update_sql(self):
        """ generate update sql """
        fields = self.__update_fields
        if fields is None:
            self.log("Need to specify update fields using loader.set_update_fields(<array>) to do an update")
            raise ValueError("Need to specify update fields using loader.set_update_fields(<array>) to do an update")
        
        updateString = ""   
        for f in fields:
            if f in JSONB_UPDATE_FIELDS:
                updateString += ' '.join((f,  '=', 'COALESCE(v.' + f, ",'{}'::jsonb)", '||', 'd.' + f + '::jsonb')) + ',' + NBSP 
                # e.g. v.gwas_flags = v.gwas_flags || d.gwas_flags::jsonb
            elif f in BOOLEAN_FIELDS:
                updateString += ' '.join((f,  '=', 'd.' + f + '::BOOLEAN')) + ',' + NBSP
            else:
                updateString += ' '.join((f,  '=', 'd.' + f)) + ',' + NBSP
        
        updateString = updateString.rstrip(', ') + NBSP
        
        if self.is_adsp() and 'is_adsp_variant' not in fields: # specified as datasource
            fields.append('is_adsp_variant')
            updateString += ", v.is_adsp_variant = d.is_adsp_variant" + NBSP

        fields = ['record_primary_key', 'chromosome'] + fields
        schema = "AnnotatedVDB"
        # partition = '_chr' + xstr(self.chromosome()) if self.chromosome() is not None else ''
        sql = "UPDATE " + schema + ".Variant v SET" + NBSP + updateString \
            + "FROM (VALUES %s) AS d(" + ','.join(fields) + ")" + NBSP \
            + "WHERE v.chromosome = d.chromosome" + NBSP \
            + "AND v.record_primary_key = d.record_primary_key"
        
        self.set_update_sql(sql)
        if self._debug:        
            self.log("Update SQL: " + self._update_sql, prefix="DEBUG")
            
            
            
    def __buffer_update_values(self, entry, flags):
        """ save udpate values to value list """
        if self._debug:
            self.log('Entering ' + type(self).__name__ + '.' + '__buffer_update_values', prefix="DEBUG")
            
        # TODO: check for duplicate has to happen in generate_update_values if not accounted for by flags
        recordPK, uFlags, uValues = self.generate_update_values(entry, flags)
        
        isAdspVariant = None
        if uFlags:
            if 'is_adsp_variant' in 'uFlags':
                isAdspVariant = uFlags['is_adsp_variant']
            
        if recordPK is None:
            if self._debug:
                self.log(("Variant", self._current_variant.id, "not in DB, inserting."), prefix="DEBUG")
            return 'INSERT'
        
        if uFlags is not None:
            if 'update' in uFlags:
                if uFlags['update'] is False:
                    if self._debug:
                        self.log((recordPK, "already has values for update fields; skipping"), prefix="DEBUG")
                    self.increment_counter('skipped')
                    return 'SKIPPED'
        
        values = [recordPK]
        values.append('chr' + self._current_variant.chromosome)
        
        if self._debug:
            self.log(("Updating variant", recordPK, " (is ADSP = ", xstr(isAdspVariant) + ")", "with", uValues), prefix="DEBUG")
            self.log(("Update fields", self.__update_fields), prefix="DEBUG")
        
        for f in self.__update_fields:
            cValue = uValues[f]
         
            if f in JSONB_UPDATE_FIELDS:
                values.append(xstr(cValue))
            else:
                values.append(uValues[f] if uValues[f] != 'NULL' else None)
            
        if self.is_adsp() and 'is_adsp_variant' not in uValues:
            values.append(True)
            
        if self._debug:
            self.log(("Buffered Values:", values), prefix="DEBUG")
        
        self._update_buffer.append(tuple(values))
        self.increment_counter('update')
        return 'UPDATE'
    
        
    def __add_adsp_update_statement(self, recordPK, chromosome):
        """ add update is_adsp_variant to update buffer """
        chrm = 'chr' + xstr(chromosome)
        values = (recordPK, chrm)
        self._update_buffer.append(values)
        # updateStr = "UPDATE AnnotatedVDB.Variant SET " \
        #    + "is_adsp_variant = true WHERE record_primary_key = '"  + recordPK + "' " \
        #    + "AND chromosome = '" + chrm + "'" # incl chromosome so it uses an index
        # self._update_buffer.write(updateStr + ';')
        
    
    def __parse_alt_alleles(self, vcfEntry, flags):
        """! iterate over alleles and calculate values to load
            add to copy buffer
            @param vcfEntry             the VCF entry for the current variant
            @param flags                info passed for script that may qualify update or load
        """
        
        if self._debug:
            self.log('Entering ' + type(self).__name__ + '.' + '__parse_alt_alleles', prefix="DEBUG")
            self.log((print_dict(self._current_variant)), prefix="DEBUG")
            self.log(("Update/Load flags:", flags), prefix="DEBUG")
            
        variant = self._current_variant # just for ease/simplicity
        variantExternalId = variant.ref_snp_id if hasattr(variant, 'ref_snp_id') else None
        
        for alt in variant.alt_alleles:              
            if alt == '.':
                self.log(("Skipping variant", self._current_variant.id, "no alt allele (alt = .)"), prefix="WARNING")
                self.increment_counter('skipped')
                continue
            
            annotator = VariantAnnotator(variant.ref_allele, alt, variant.chromosome, variant.position)       
            
            if self.skip_existing():
                if self.is_duplicate(annotator.get_metaseq_id(), 'METASEQ'):
                    self.increment_counter('skipped')
                    continue
            
            status = None
            if self.update_existing():
                status = self.__buffer_update_values(vcfEntry, flags)
                if status != 'INSERT':
                    continue # skipped or updated
                
            # other wise -- new variant load
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
            
            additionalValues = None
            if self.__update_fields is not None:
                pkPlaceholder, uFlags, additionalValues = self.generate_update_values(vcfEntry, flags)                

            # TODO iterate over copy fields and add as each one is hit so order is preserved
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
            
            if additionalValues:
                for f in self.__update_fields:
                    copyValues.append(xstr(additionalValues[f], nullStr='NULL'))
            
            if self.is_adsp():
                copyValues.append(xstr(True)) # is_adsp_variant
                         
            if self._debug:
                self.log("Display Attributes " + print_dict(annotator.get_display_attributes(variant.rs_position)), prefix="DEBUG")
                self.log(("COPY values", copyValues), prefix="DEBUG")
                
            self.add_copy_str('#'.join(copyValues))
            self.increment_counter('variant')
        
        
    def parse_variant(self, line, flags=None):
        """! parse & load single line from file 
        @param line             the VCF entry in string form
        @param flags            flags that may be used to modify update or copy strings
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
            entry = VcfEntryParser(line, self.vcf_header_fields()) if isinstance(line, str) else line # assume vcfEntry was passed
            
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
                    self.increment_counter('skipped')
                    return None
            
            # iterate over alleles
            self.__parse_alt_alleles(entry, flags)
            
        except Exception as err:
            self.log(str(err), prefix="ERROR")
            raise err