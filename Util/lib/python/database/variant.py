"""! @brief AnnotatedVDB Variant Lookups """

#!pylint: disable=invalid-name
##
# @package database
# @file variant.py
#
# @brief  Annotated VDB Variant Lookup
#
# @section variant Description
# provides functions for looking up variants in the database, to validate or find duplicates
#
# @section todo_variant TODO
# - type -> sql mappings as an enum?
#
# @section libraries_variant Libraries/Modules
# 
# - [GenomicsDBData.Util.utils](https://github.com/NIAGADS/GenomicsDBData/blob/master/Util/lib/python/utils.py)
#   + provides variety of wrappers for standard file, string, and logging operations
# - [GenomicsDBData.Util.list_utils](https://github.com/NIAGADS/GenomicsDBData/blob/master/Util/lib/python/list_utils.py)
#   + provides variety list and set operations
# - [GenomicsDBData.Util.postgres_dbi](https://github.com/NIAGADS/GenomicsDBData/blob/master/Util/lib/python/postgres_dbi.py)
#   + wrapper for connecting to database & handling PG errors
# 
# @section author_variant Author(s)
# - Created by Emily Greenfest-Allen (fossilfriend) 2022

import json
import logging
from niagads.db.postgres import Database, raise_pg_exception
from niagads.utils.string import xstr

VARIANT_ID_TYPES = ['REFSNP', 'METASEQ', 'PRIMARY_KEY']

PRIMARY_KEY_LOOKUP_SQL = "record_primary_key = %s AND chromosome = 'chr' || split_part(%s, ':', 1)"
LEGACY_PRIMARY_KEY_LOOKUP_SQL = """LEFT(metaseq_id, 50) = split_part(%s, '_', 1)
    AND (ref_snp_id = split_part(%s, '_', 2) 
    OR (ref_snp_id IS NULL AND split_part(%s, '_', 2) = ''))"""
    
BULK_LOOKUP_FULL_SQL = "SELECT * from get_variant_primary_keys_and_annotations(%s, %s, %s)"; # second param is "firstValueOnly flag"
BULK_LOOKUP_SQL = "SELECT * from map_variants(%s, %s, %s)"; # second param is "firstValueOnly flag", then checkAltVariants

class VariantRecord(object):
    """! functions finding and validating against variants in the DB    
    creates one database cursor for each type of lookups to avoid locks
    """
    
    def __init__(self, gusConfigFile=None, connectionString=None, verbose=False, debug=False):
        """! VariantRecord base class initializer
        @param gusConfigFile          file containing connection info for the DB, following gus.config format
        @param verbose                flag for verbose output
        @param debug                  flag for debug output
        @returns                      An instance of the VariantLoader class with initialized database cursor
        """
        self._debug = debug
        self._verbose = verbose
        self.logger = logging.getLogger(__name__)
        self.__legacy_pk = False # use legacy dynamic PK
        self.__database = None
        self.__cursors  = {}                # possibly one for each ID TYPE
                
        self.__initialize_database(gusConfigFile)
        
        
    def __initialize_database(self, gusConfigFile):
        """! establish databse connection
            @param gusConfigFile             file containing connection info for the DB, following gus.config format
            @exception database exceptions handled by raise_pg_exception
        """
        try:
            self.__database = Database(gusConfigFile)
            self.__database.connect()
        except Exception as err:
            raise_pg_exception(err)
            
    
    def close_cursor(self, cursorKey):
        """
        close cursor by key
        """
        self.__cursors[cursorKey].close()
    
    
    def close_cursors(self):
        """close any open database cursors"""
        for cursor in self.__cursors.values():
            cursor.close()
            
        
    def close_database(self):
        """! close the database connection """
        self.close_cursors()
        self.__database.close()

        
    def close(self):
        """! close the VariantRecord object """
        self.close_database()
    
    
    def rollback(self):
        self.__database.rollback()
        
        
    def __exit__(self, exc_type, exc_value, traceback):
        """! garbage collection """
        self.close()
        
        
    def __validate_id_type(self, idType):
        """! validate provided id type
            @param idType             ID type string
            @exception ValueError          an attribute error if invalid
        """
        if idType.upper() not in VARIANT_ID_TYPES:
            raise ValueError(idType + ' not a valid variant ID type: ' + VARIANT_ID_TYPES)
        
        
    def initialize_named_cursor(self, name: str, realDict=False, withhold=False):
        if ' ' in name:
            raise ValueError("Invalid name " + name + " for cursor; no spaces allowed")
        cursorFactory = 'RealDictCursor' if realDict else None
        self.__cursors[name] = self.__database.named_cursor(name, cursorFactory=cursorFactory, withhold=withhold)
        
        
    def initialize_cursor(self, cursorKey, realDict=False):
        """! initialize a database cursor, and reference by key
            @param cursorKey             string key to identify the cursor in the cursor dictionary
            @param realDict              whether cursor should be a psycopg2 'RealDictCursor'
        """
        cursorFactory = 'RealDictCursor' if realDict else None
        self.__cursors[cursorKey] = self.__database.cursor(cursorFactory)


    def use_legacy_pk(self, useLegacyPK):
        self.__legacy_pk = useLegacyPK


    def legacy_pk(self):
        return self.__legacy_pk

    
    def get_cursor(self, cursorKey, initializeIfMissing=False, realDict=False):
        """! get / initialize a database cursor

            @param cursorKey             the dict key for this cursor
            @exception ValueError            if initializeIfMissing is False, will raise error if key does not exist
            @returns                     the database cursor
        """
        if cursorKey not in self.__cursors:
            if initializeIfMissing:
                self.initialize_cursor(cursorKey, realDict=realDict)
            else:
                raise ValueError("Dictionary Cursor: " + cursorKey + " not initialized")
    
        return self.__cursors[cursorKey]
    

    def bulk_lookup(self, variants, firstHitOnly=True, fullAnnotation=True, checkAltVariants=True):
        """! lookups up a list of variants in the DB / takes metaseq_ids or refsnps 
        and returns the folling JSON: (note this eg is mapped against the GRCh37 DB, but the idea is the same)
            "1:1510801:C:T": {
                "bin_index": "chr1.L1.B1.L2.B1.L3.B1.L4.B1.L5.B1.L6.B1.L7.B2.L8.B2.L9.B1.L10.B1.L11.B1.L12.B1.L13.B1",
                "annotation": {
                    "GenomicsDB": ["ADSP_WGS","NG00027_STAGE1"],
                    "mapped_coordinates": null,
                    "ADSP_QC" (w/out allele counts & freqs)
                } ,
                "match_rank": 1,
                "match_type": "exact",
                "metaseq_id": "1:1510801:C:T",
                "ref_snp_id": "rs7519837",
                "is_adsp_variant": true,
                "record_primary_key": "1:1510801:C:T_rs7519837" 
            },
            
        """
        cursorKey = 'BULK_VARIANT_LOOKUP'
        sql = BULK_LOOKUP_FULL_SQL if fullAnnotation else BULK_LOOKUP_SQL
        cursor = self.get_cursor(cursorKey, initializeIfMissing=True, realDict=False)
        if not isinstance(variants, str): # assume array or tuple
            variants = ','.join(variants)
            
        params = [variants, firstHitOnly, checkAltVariants]
        cursor.execute(sql, tuple(params))
        result = cursor.fetchone()[0]
        
        if self._debug:
            self.logger.debug(xstr({'variants': variants, 'result': result}))
            
        return json.loads(result) if isinstance(result, str) else result
    
    
    def has_json_attr(self, field, key, variantPK, returnVal=True):
        """! checks to see if a JSONB field contains a specific key, returns value if true

            @param field            field/column to query
            @param key              json key fiekd
            @param variantPK           variant identifier

            @returns tuple that is (primary key, value of the json attribute/ None if null) or (primary_key, boolean)
        """       
        
        raise NotImplementedError("TODO: update code to reflect PK only lookups")
        idType = idType.upper()
        self.__validate_id_type(idType)

        try: 
            lookupIdType = 'LEGACY_PRIMARY_KEY' if idType == 'PRIMARY_KEY' and self.legacy_pk() else idType
            cursorKey = 'HAS_JSON_ATTR_' + field + '_' + lookupIdType
            sql = "SELECT record_primary_key, " + field + "->>'" + key \
                + "' FROM AnnotatedVDB.Variant WHERE " + LOOKUP_SQL[lookupIdType]
                
            if self.legacy_pk():
                sql = sql.replace('AnnotatedVDB', 'Public')
                sql = sql.replace('SELECT record_primary_key',
                    "SELECT COALESCE(LEFT(metaseq_id, 50) || '_' || ref_snp_id', LEFT(metaseq_id, 50)) AS record_primary_key")
                
            cursor = self.get_cursor(cursorKey, initializeIfMissing=True, realDict=False)
            numBindParams = sql.count('%s')
            params = []
            for i in range(0, numBindParams):
                params.append(variantId)
            
            if chromosome is not None: # for refsnp lookups
                sql += ' AND chromosome = %s'
                chrm = str(chromosome) if 'chr' in str(chromosome) else 'chr' + str(chromosome)
                params.append(chrm)
                
            if chromosome is None and lookupIdType != 'REFSNP':
                sql += " AND chromosome = 'chr' || split_part(%s, ':', 1)"
                params.append(variantId)
                
            cursor.execute(sql, tuple(params))
            result = cursor.fetchone()

            if result is None:
                return None # variant not in db
            
            return list(result) if returnVal else (result[0], result is not None)

        except Exception as err:
            raise_pg_exception(err)
            
        return None
    
    
    def has_attr(self, field, variantPK, returnVal=True):
        """! checks to see the variant has a value for one or more fields
        usually used before an update so returning PK incase lookup is not the PK
        NOTE: will only return True if all values are None, so getting the value list best if 
        (returnVal = True) if assessing multiple values simultaneously
        
            @param fields                 one or more fields (string or array) /column to query
            @param variantPK             variant primary key
            @param returnVal             return the value if True else return boolean flag
            @returns    value of the attribute or bool if returnVal = False
        """
        
        queryFields = field
        if not isinstance(field, str): # then array of strings
            queryFields = ','.join(field)
            returnVal = True # must return the values if requesting multiple fields
            
        try: 
            cursorKey = 'HAS_ATTR_LOOKUP'
            sql = "SELECT " + queryFields + " FROM AnnotatedVDB.Variant WHERE record_primary_key = %s"
                
            cursor = self.get_cursor(cursorKey, initializeIfMissing=True, realDict=False)
            cursor.execute(sql, (variantPK, ))
            result = cursor.fetchone()
        
            if result is None:
                raise KeyError("No record found for variant %s in database.")
            
            if len(result) == 1:
                return result[0] if returnVal else result[0] is not None
            else:
                return list(result) if returnVal else all(v is not None for v in result)


        except Exception as err:
            raise err
            
    

    def exists(self, variantId, returnMatch=False):
        """! check if against the AnnotatedVDB to see if a variant with the specified id is already present

            @param variantId             variant id to lookup
            @param idType                type of variant identifier, must match one of VARIANT_ID_TYPES
            @param returnMatch           return the match (pk & bin index)
            @exception PG / DatabaseError   if issue w/database lookup, handled by raise_pg_exception
            @returns                     True if the variant exists in the database, else False
        """
        try: 
            cursorKey = 'variant_exists'
            cursor = self.get_cursor(cursorKey, initializeIfMissing=True, realDict=False)

            cursor.execute(BULK_LOOKUP_SQL, (variantId, True, True)) # firstHitOnly, checkAltAlleles
            result = cursor.fetchone()
            
            if result[0][variantId] is None:
                return None if returnMatch else False # variant not in db
            
            return result[0][variantId] if returnMatch else True

        except Exception as err:
            raise err


