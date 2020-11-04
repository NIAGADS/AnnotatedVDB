#!/usr/bin/env python
#pylint: disable=multiple-statements,invalid-name

from GenomicsDBData.Util.utils import xstr, warning
from GenomicsDBData.Util.postgres_dbi import Database
from psycopg2.extras import NumericRange
from psycopg2 import DatabaseError, ProgrammingError

BIN_INDEX_SQL = """
WITH bi AS (SELECT find_bin_index(%s, %s, %s) AS global_bin_path)
SELECT br.chromosome, br.global_bin_path, br.location, nlevel(br.global_bin_path) AS bin_level
FROM BinIndexRef br, bi
WHERE br.global_bin_path = bi.global_bin_path
"""

class BinIndex(object):
    ''' class to fetch and map locations to index bins 
    creates its own DB handle

    b/c variants are usually incremental; stores current bin
    and checks if variant is in that; if not fetches bin from DB
    '''

    def __init__(self, gusConfigFile=None, verbose=True):
        self._gusConfigFile = gusConfigFile
        self._verbose = verbose
        self._currentBin = {}
        self.__get_db_handle()
        self._cursor = self._database.cursor("RealDictCursor")


    def __get_db_handle(self):
        ''' create database connection '''
        self._database = Database(self._gusConfigFile)
        self._database.connect()


    def close(self):
        ''' close db connection '''
        self._cursor.close()
        self._database.close()

    def _update_current_bin_index(self, chrm, start, end):
        ''' query against database to get minimum enclosing bin;
        set as new current bin'''
        if self._verbose: warning("Updating current bin")
        result = None
        self._cursor.execute(BIN_INDEX_SQL, (chrm, start, end))
        try:
            self._currentBin = self._cursor.fetchone()

        except ProgrammingError:
            raise ProgrammingError('Could not map ' + chrm + ':' + xstr(start) + '-' + xstr(end) + ' to a bin.')

        if self._verbose: warning(self._currentBin)
        return result


    def find_bin_index(self, chrm, start, end=None):
        ''' finds the bin index for the position;
        if end == None assume SNV, and set to start'''

        if end is None: end = start
        if 'chr' not in chrm: chrm = 'chr' + xstr(chrm)

        if bool(self._currentBin): # if a current bin exists and the variant falls in it, return it
            if self._currentBin['bin_level'] >= 27: # otherwise may be a broad bin b/c of a indel; so need to do a lookup
                brange = self._currentBin['location']
                if self._currentBin['chromosome'] == chrm \
                  and start in brange and end in brange:
                    return self._currentBin['global_bin_path']

        # otherwise, find & return the new bin
        self._update_current_bin_index(chrm, start, end) # not in current bin, so update bin
        return self._currentBin['global_bin_path']
