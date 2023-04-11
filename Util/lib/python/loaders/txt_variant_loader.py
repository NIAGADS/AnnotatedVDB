"""! @brief Annotated VDB Text/CSV Loader """
#!pylint: disable=invalid-name

##
# @package loaders
# @file txt_variant_loader.py
#
# @brief  Annotated VDB Text/CSV Loader
#
# @section txt_variant_loader Description
# provides functions & class for managing variant loads from a delimited text file
# where column headers match fields in the AnnotatedVDB.Variant table
#
# @section todo_txt_variant_loader TODO
#
# @section libraries_txt_variant_loader Libraries/Modules
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
# - [AnnotatedVDB.Util.database](https://github.com/NIAGADS/AnnotatedVDB/blob/master/Util/lib/python/database)
#   + variant lookups
# - [AnnotatedVDB.Util.variant_annotator](https://github.com/NIAGADS/AnnotatedVDB/blob/master/Util/lib/python/variant_annotator.py)
#   + generate standard variant annotations from position and allele information
# @section author_txt_variant_loader Author(s)
# - Created by Emily Greenfest-Allen (fossilfriend) 2022

import json
from types import SimpleNamespace

from GenomicsDBData.Util.utils import xstr, warning, print_dict, to_numeric, deep_update
from GenomicsDBData.Util.list_utils import qw, is_subset, is_equivalent_list
from GenomicsDBData.Util.postgres_dbi import Database, raise_pg_exception

from AnnotatedVDB.Util.parsers import VcfEntryParser, CONSEQUENCE_TYPES
from AnnotatedVDB.Util.variant_annotator import VariantAnnotator
from AnnotatedVDB.Util.loaders import VariantLoader, REQUIRED_COPY_FIELDS, ALLOWABLE_COPY_FIELDS, JSONB_UPDATE_FIELDS
from AnnotatedVDB.Util.database import VARIANT_ID_TYPES

NBSP = " " # for multi-line sql

class TextVariantLoader(VariantLoader):
    """! functions for loading variants from a delimited text file """
    
    def __init__(self, datasource, logFileName=None, verbose=False, debug=False):
        """! VCFVariantLoader base class initializer

            @param datasource          datasource description
            @param logFileName         full path to logging file
            @param verbose             flag for verbose output
            @param debug               flag for debug output
            
            @returns                   An instance of the TextVariantLoader class with initialized counters and fields
        """
        super(TextVariantLoader, self).__init__(datasource, logFileName, verbose, debug)
        self.__update_existing = False
        self.__required_copy_fields = None
        self.__update_fields = None
        self.__variant_id_type = None
        self.__update_existing = False
        self.log((type(self).__name__, "initialized"), prefix="INFO")
 
 
    def set_update_fields(self, fields):
        self.__update_fields = fields       
 
    
    def set_variant_id_type(self, variantIdType):
        self.__variant_id_type = variantIdType
        
        
    def variant_id_type(self):
        return self.__variant_id_type
    
    
    def set_update_existing(self, updateExisting):
        """! set update existing flag"""
        self.__update_existing = updateExisting
        
        
    def update_existing(self):
        """! @returns update existing flag"""
        return self.__update_existing
    
        
    def initialize_copy_sql(self, headerFields):
        """! _summary_
        @param headerFields                        fields to add, from the header
        """
        copyFields = REQUIRED_COPY_FIELDS # "chromosome", "record_primary_key", "position", "metaseq_id", "bin_index", "row_algorithm_id" , display_attributes]
        copyFields.extend(["display_attributes"])

        self.__update_fields = []
        for f in headerFields:
            if f in ALLOWABLE_COPY_FIELDS and f not in copyFields:
                copyFields.append(f)
                self.__update_fields.append(f)

        if self.is_adsp():
            copyFields.append('is_adsp_variant')
        
        self.__required_copy_fields = copyFields   
        
        if self._debug:
            self.log(copyFields, prefix="DEBUG")

        super(TextVariantLoader, self).initialize_copy_sql(copyFields=copyFields)
        
        
    def build_update_sql(self, useLegacyPk=False):
        """ generate update sql """
         
        fields = self.__update_fields
        updateString = ""
        for f in fields:
            if f in JSONB_UPDATE_FIELDS:
                updateString += ' '.join((f,  '=', 'COALESCE(v.' + f, ",'{}'::jsonb)", '||', 'd.' + f + '::jsonb')) + ', ' # e.g. v.gwas_flags = v.gwas_flags || d.gwas_flags::jsonb
            elif f == 'bin_index':
                updateString += ' '.join((f,  '=', 'd.' + f + '::ltree')) + ', '
            else:
                updateString += ' '.join((f,  '=', 'd.' + f)) + ', '
        
        updateString = updateString.rstrip(', ') + NBSP
        
        if self.is_adsp():
            fields.append('is_adsp_variant')
            updateString += ", v.is_adsp_variant = d.is_adsp_variant" + NBSP

        fields = ['record_primary_key', 'chromosome'] + fields
        schema = "Public" if useLegacyPk else "AnnotatedVDB"
        sql = "UPDATE " + schema + ".Variant v SET" + NBSP + updateString \
            + "FROM (VALUES %s) AS d(" + ','.join(fields) + ")" + NBSP \
            + "WHERE v.chromosome = d.chromosome" + NBSP
        if useLegacyPk:
            sql += "AND LEFT(v.metaseq_id, 50) = split_part(d.record_primary_key, '_', 1)" + NBSP \
                + "AND (v.ref_snp_id = split_part(d.record_primary_key, '_', 2)" + NBSP \
                + "OR (v.ref_snp_id IS NULL AND split_part(d.record_primary_key, '_', 2) = ''))"
        else:
            sql += "AND v.record_primary_key = d.record_primary_key"
        
        self.set_update_sql(sql)
        if self._debug:        
            self.log("Update SQL: " + self._update_sql, prefix="DEBUG")


    def __get_current_variant_id(self):
        """ get id and idType """
        if hasattr(self._current_variant, 'record_primary_key'):
            return (self._current_variant.record_primary_key, 'PRIMARY_KEY')
        
        if hasattr(self._current_variant, 'metaseq_id'):
            return (self._current_variant.metaseq_id, 'METASEQ')
        
        else:
            return (self._current_variant.ref_snp_id, 'REFSNP')
        

    def set_current_variant(self, record):
        """ variantInfo should be a simple namespace or a dict"""   
        
        variantId = record['variant']
        variantInfo = {'id': variantId}
        idElements = variantId.split(':') if ':' in variantId else None

        if self.variant_id_type() == 'METASEQ':
            variantInfo.update({'metaseq_id': variantId})
            
        if self.variant_id_type() == 'REFSNP':
            variantInfo.update({'ref_snp_id': record['ref_snp_id']})
                
        if self.variant_id_type() == 'PRIMARY_KEY':
            variantInfo.update({'record_primary_key': variantId})
            
        if idElements is not None:
            variantInfo.update({'chromosome': 'chr' + xstr(idElements[0])})
                        
        self._current_variant = SimpleNamespace(**variantInfo)
    
    
    def __is_updatable_json_element(self, updateStr, dbVals):
        if dbVals is None: # field is empty in DB
            return True
        
        updateVals = json.loads(updateStr)  
        for uKey in updateVals.keys():
            if uKey in dbVals:
                uVal = updateVals[uKey]
                dVal = dbVals[uKey] if uKey in dbVals else None
                if dVal is None: # not in db
                    return True
                else:
                    if uVal != dVal: # in db but not the same
                        self.log(("u", uVal, "db", dVal), prefix="DEBUG")
                        self.log(("match", uVal == dVal), prefix="DEBUG")
                        return True
            
        return False   

    
    def __buffer_update_values(self, recordPK, record):
        """ save udpate values to value list """
        
        values = [recordPK]
        if hasattr(self._current_variant, 'chromosome'):
            values.append(self._current_variant.chromosome)
        else:
            pkArray = recordPK.split(':')
            values.append('chr' + xstr(pkArray[0]))

        for f in self.__update_fields:
            uValue = record[f] 
            values.append(uValue if uValue != 'NULL' else None)
            
        if self.is_adsp():
            values.append(True)
        
        self._update_buffer.append(tuple(values))
        self.increment_counter('update')
        
    
    def __add_to_copy_buffer(self, record):
        """ generate copy buffer string and add to buffer """

        variant = self._current_variant # just for ease/simplicity
        variantExternalId = variant.ref_snp_id if hasattr(variant, 'ref_snp_id') else None
        annotator = VariantAnnotator(variant.ref_allele, variant.alt_allele, variant.chromosome, variant.position)       
        
        binIndex = self._bin_indexer.find_bin_index(variant.chromosome, variant.position,
                                                    annotator.infer_variant_end_location())

        recordPK = self._pk_generator.generate_primary_key(annotator.get_metaseq_id(), variantExternalId)      
    
        # "chromosome", "record_primary_key", "position", "metaseq_id", "bin_index", "row_algorithm_id" , display_attributes]
        copyValues = ['chr' + xstr(variant.chromosome),
                    recordPK,
                    xstr(variant.position),
                    annotator.get_metaseq_id(),
                    binIndex,
                    xstr(self._alg_invocation_id),
                    xstr(annotator.get_display_attributes(variant.rs_position), nullStr='NULL')     
                    ]      
        
        for f in self.__update_fields:
            copyValues.append(xstr(record[f], nullStr='NULL'))
        
        if self.is_adsp():
            copyValues.append(xstr(True)) # is_adsp_variant
            
        if self._debug:
            self.log("Display Attributes " + print_dict(annotator.get_display_attributes(variant.rs_position)), prefix="DEBUG")
            self.log(copyValues, prefix="DEBUG")
            
        self.add_copy_str('#'.join(copyValues))
        
        
    def parse_variant(self, record):
        """! update & load single record from file / expects to be read using CSV reader, so should be a dict
        @param record                      record from the file in dict form
        """
        
        if isinstance(record, str):
            raise TypeError("Please use the csv package to parse the file so that values can be accessed by column names")
            
        if self.resume_load() is False and self._resume_after_variant is None:
            err = ValueError('Must set VariantLoader resume_afer_variant if resuming load')
            self.log(str(err), prefix="ERROR")
            raise err
        
        try:
            self.increment_counter('line')
            self.set_current_variant(record)
            
            if not self.resume_load():
                self._update_resume_status(self._current_variant.id)
                return None

            recordPK, recordIdType = self.__get_current_variant_id()
            
            if recordIdType != 'PRIMARY_KEY': # get the primary key
                recordPK = self.is_duplicate(recordPK, recordIdType, returnPK=True)

            if recordPK is not None: # duplicate
                self.increment_counter('variant')
                self.increment_counter('duplicates')
                if self.skip_existing():
                    self.increment_counter('skipped')
                    return None
                if self.update_existing():
                    # self.increment_counter('update') -- now done in __buffer_update_values
                    self.__buffer_update_values(recordPK, record)
            
            else:
                self.increment_counter('variant')
                self.__add_to_copy_buffer(record)

        except Exception as err:
            self.log(str(err), prefix="ERROR")
            raise err