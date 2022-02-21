'''
utils for parsing CADD tbx file
'''
#pylint: disable=line-too-long,invalid-name,method-hidden

##
# @package util

import pysam
import json
from os import path
from GenomicsDBData.Util.utils import qw, die, warning, xstr

CADD_INDEL_FILE = "InDels.tsv.gz"
CADD_SNV_FILE = "whole_genome_SNVs.tsv.gz"


class CADDUpdater(object):
    ''' utils for slicing CADD tbx database '''

    def __init__(self, logFileName, databasePath):
        self._path = databasePath
        self._snv = self._connect(CADD_SNV_FILE)
        self._indel = self._connect(CADD_INDEL_FILE)
        self._sql_buffer = []
        self._snv_update_count = 0
        self._indel_update_count = 0
        self._chrm = None
        self._log_file = logFileName
        self._lfh = None 
        self._update_count = 0 
        self._not_matched_count = 0


    def setChrm(self, chrm):
        self._chrm = chrm;


    def buffered_variant_count(self):
        return len(self._sql_buffer)

    
    def get_update_count(self, indels=False):
        if indels:
            return self._indel_update_count
        return self._snv_update_count


    def get_total_update_count(self):
        return self._update_count + self._not_matched_count


    def increment_not_matched_count(self):
        self._not_matched_count = self._not_matched_count + 1


    def get_not_matched_count(self):
        return self._not_matched_count


    def increment_update_count(self, indels=False):
        if indels:
            self._indel_update_count = self._indel_update_count + 1
        else:
            self._snv_update_count = self._snv_update_count + 1
        self._update_count = self._update_count + 1

    def lfh(self):
        if self._lfh is None:
            self._lfh = open(self._log_file, 'w')
        
        return self._lfh


    def close_lfh(self):
        if self._lfh is not None:
            self._lfh.close()


    def _connect(self, fileName):
        tbxFile = path.join(self._path, fileName)
        return pysam.TabixFile(tbxFile)


    def svn(self):
        return self._snv


    def indel(self):
        return self._indel


    def buffer_update_sql(self, metaseqId, evidence):
        chrm, pos, ref, alt = metaseqId.split(':')
        if self._chrm is not None:
            chrm = self._chrm
        else:
            chrm = 'chr' + xstr(chrm)

        sql = "UPDATE Variant_" + chrm \
          + " v SET cadd_scores = '" + json.dumps(evidence) + "'" \
          + " WHERE left(v.metaseq_id, 50) = left('" + metaseqId + "',50)" \
          + " AND chromosome = '" + chrm + "'" \
          + " AND v.metaseq_id = '" + metaseqId + "'"
        self._sql_buffer.append(sql)


    def clear_update_sql(self):
        self._sql_buffer = []


    def sql_buffer_str(self):
        return ';'.join(self._sql_buffer)


    def _slice(self, chrm, start, end, indels=False):
        ''' slice CADD file; assume using this more
        often with snv file; returns tuple'''
        tbxFh =  self._indel if indels else self._snv
        caddChr = 'MT' if chrm == 'M' else xstr(chrm)
        return tbxFh.fetch(caddChr, start, end, parser=pysam.asTuple())

  
    def match(self, chrm, position, indels=True):
        ''' match single position in CADD file; default indels to true
        b/c will be using this more with the indel file'''
        tbxFh = self._indel if indels else self._snv
        caddChr = 'MT' if chrm == 'M' else xstr(chrm)
        try:
            return tbxFh.fetch(caddChr, int(position) - 1, int(position), parser=pysam.asTuple())
        except ValueError as e: # happens sometimes on chrm M/MT
            warning("WARNING", e)
            return []
