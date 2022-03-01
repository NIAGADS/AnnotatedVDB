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

import sqlite3
from GenomicsDBData.Util.postgres_dbi import Database, raise_pg_exception

VARIANT_ID_TYPES = ['REFSNP', 'METASEQ', 'PRIMARY_KEY']

LOOKUP_SQL="SELECT record_primary_key, metaseq_id, ref_snp_id FROM AnnotatedVDB.Variant WHERE"
EXISTS_SQL="SELECT record_primary_key FROM AnnotatedVDB.Variant WHERE"
METASEQ_LOOKUP_SQL="LEFT(metaseq_id, 50) = LEFT(%s, 50) AND metaseq_id = %s AND chromosome = 'chr' || split_part(%s, ':', 1)" 
REFSNP_LOOKUP_SQL="ref_snp_id = %s"
PRIMARY_KEY_LOOKUP_SQL="record_primary_key = %s AND chromosome = 'chr' || split_part(%s, ':', 1)"

LOOKUP_SQL = {
    'REFSNP' : REFSNP_LOOKUP_SQL,
    'METASEQ' : METASEQ_LOOKUP_SQL,
    'PRIMARY_KEY' : PRIMARY_KEY_LOOKUP_SQL
}

class VariantRecord(object):
    """! functions finding and validating against variants in the DB    
    creates one database cursor for each type of lookups to avoid locks
    """
    
    def __init__(self, gusConfigFile=None, verbose=False, debug=False):
        """! VariantRecord base class initializer
        @param gusConfigFile          file containing connection info for the DB, following gus.config format
        @param verbose                flag for verbose output
        @param debug                  flag for debug output
        @returns                      An instance of the VariantLoader class with initialized database cursor
        """
        self.__debug = debug
        self.__verbose = verbose
        self.__database = None
        self.__cursors = {}                # possibly one for each ID TYPE
                
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
            
    
    def close_cursors(self):
        """! close any open database cursors"""
        for cursor in self.__cursors.values():
            cursor.close()
            
        
    def close_database(self):
        """! close the database connection """
        self.close_cursors()
        self.__database.close()

        
    def close(self):
        """! close the VariantRecord object """
        self.close_database()
    
        
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
        
        
    def initialize_cursor(self, cursorKey, realDict=False):
        """! initialize a database cursor, and reference by key
            @param cursorKey             string key to identify the cursor in the cursor dictionary
            @param realDict              whether cursor should be a psycopg2 'RealDictCursor'
        """
        cursorFactory = 'RealDictCursor' if realDict else None
        self.__cursors[cursorKey] = self.__database.cursor(cursorFactory)

    
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
    
    
    def has_attr(self, field, variantId, idType, returnVal=True):
        """! checks to see the variant has a value for that field
        usually used before an update so returning PK incase lookup is not the PK
        
            @param field                field/column to query
            @param variantId             variant identifier
            @param idType                type of variant identifier, must match one of VARIANT_ID_TYPES
            @param returnVal             return the value if True else return boolean flag
            @returns    tuple that is (primary key, value of the attribute/ None if null) or (primary_key, boolean)
        """
        idType = idType.upper()
        self.__validate_id_type(idType)
      
        try: 
            cursorKey = 'HAS_ATTR_' + field + '_' + idType
            sql = "SELECT record_primary_key, " + field + " FROM AnnotatedVDB.Variant WHERE " + LOOKUP_SQL[idType]
            cursor = self.get_cursor(cursorKey, initializeIfMissing=True, realDict=False)
            numBindParams = sql.count('%s')
            params = []
            for i in range(0, numBindParams):
                params.append(variantId)
            
            cursor.execute(sql, tuple(params))
            result = cursor.fetchone()
            
            if result is None:
                return None # variant not in db
            
            return result if returnVal else (result[0], result is not None)

        except Exception as err:
            raise_pg_exception(err)
            
        return None
        
        
    
        
    def exists(self, variantId, idType, chromosome=None):
        """! check if against the AnnotatedVDB to see if a variant with the specified id is already present

            @param variantId             variant id to lookup
            @param idType                type of variant identifier, must match one of VARIANT_ID_TYPES
            @exception PG / DatabaseError   if issue w/database lookup, handled by raise_pg_exception
            @returns                     True if the variant exists in the database, else False
        """
        idType = idType.upper()
        self.__validate_id_type(idType)
       
        try:
            cursorKey = 'EXISTS_' + idType
            cursor = self.get_cursor('EXISTS', initializeIfMissing=True, realDict=False)
            sql = EXISTS_SQL + ' ' + LOOKUP_SQL[idType]

            numBindParams = sql.count('%s')
            params = []
            for i in range(0, numBindParams):
                params.append(variantId)
            
            if chromosome is not None: # for refsnp lookups
                sql += ' AND chromosome = %s'
                params.append(chromosome)
            
            cursor.execute(sql, tuple(params))
            for record in cursor:
                return True # if there is a hit, the variant is already in the database

        except Exception as err:
            raise_pg_exception(err)

        return False

