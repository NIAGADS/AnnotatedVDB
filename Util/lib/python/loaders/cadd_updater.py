"""
utils for parsing CADD tbx file
"""
#pylint: disable=line-too-long,invalid-name,method-hidden

##
# @package loaders

import pysam
import json
from os import path
from types import SimpleNamespace
from GenomicsDBData.Util.utils import die, warning, xstr
from AnnotatedVDB.Util.loaders import VariantLoader

CADD_INDEL_FILE = "gnomad.genomes.r3.0.indel.tsv.gz"
CADD_SNV_FILE = "whole_genome_SNVs.tsv.gz"
NBSP = " " # for multi-line sql

CADD_UPDATE_SQL = """UPDATE AnnotatedVDB.Variant v SET cadd_scores = d.cadd_scores 
  FROM (VALUES %s) AS d(record_primary_key, chromosome, cadd_scores)
  WHERE v.record_primary_key = d.record_primary_key and v.chromosome = d.chromosome"""

class CADDUpdater(VariantLoader):
    """ utils for slicing CADD tbx database """

    def __init__(self, databasePath, logFileName=None, verbose=False, debug=False):
        super(CADDUpdater, self).__init__('CADD', logFileName, verbose, debug)
        self.__path = databasePath
        self.__snv = self.__connect(CADD_SNV_FILE)
        self.__indel = self.__connect(CADD_INDEL_FILE)     
        self.__chromosome = None
        self._initialize_counters(['snv', 'indel', 'not_matched'])
        self.set_update_sql(CADD_UPDATE_SQL)
        self.log((type(self).__name__, "initialized"), prefix="INFO")


    def set_current_variant(self, variantInfo):
        """ variantInfo should be a simple namespace or a dict"""
        self._current_variant = SimpleNamespace(**variantInfo) \
            if  isinstance(variantInfo, dict) else variantInfo


    def set_chromosome(self, chrm):
        self.__chromosome = xstr(chrm).replace('chr', '');


    def get_update_count(self, indels=False):
        if indels:
            return self.get_count('indel')
        return self.get_count('snv')


    def get_total_update_count(self):
        return self.get_count('update') + self.get_count('not_matched')


    def increment_update_count(self, indels=False):
        if indels:
            self.increment_counter('indel')
        else:
            self.increment_counter('snv')
        self.increment_counter('update')


    def __connect(self, fileName):
        tbxFile = path.join(self.__path, fileName)
        return pysam.TabixFile(tbxFile)


    def snv(self):
        return self.__snv


    def indel(self):
        return self.__indel


    def buffer_update_values(self, recordPK, evidence):
        """ save update values to value list """ 
        if self.__chromosome is not None:
            chrm = self.__chromosome
        else:
            values = recordPK.split(':')
            chrm = 'chr' + xstr(values[0])
        
        values = (recordPK, chrm, json.dumps(evidence))
        self._update_buffer.append(values)


    def buffer_update_statement(self, recordPK, evidence):
        """ generate the update statement and add to buffer """
        if self.__chromosome is not None:
                chrm = self.__chromosome
        else:
            values = recordPK.split(':')
            chrm = 'chr' + xstr(values[0])
            
        updateSql = "UPDATE AnnotatedVDB.Variant" + NBSP \
            + "SET cadd_scores = '" + json.dumps(evidence) + "'" + NBSP \
            + "WHERE record_primary_key = '" + recordPK + "'"  + NBSP \
            + "AND chromosome = '" + chrm + "'"
            
        if self._debug:
            self.log("UpdateSQL: " + updateSql, prefix="DEBUG")    
            
        self._update_buffer.write(updateSql + ";")


    def is_cadd_annotated(self):
        """ lookup variant in db
        return None if variant does not exists, set primary key and return flag if value is set"""
        result = self._variant_validator.has_attr('cadd_scores', self._current_variant.metaseq_id, 'METASEQ')
        if result is None:
            return None 
        self._current_variant.record_primary_key = result[0] # recordPK
        return result[1] # flag indicating whether value is assinged


    def slice(self, chrm, start, end, indels=False):
        """ slice CADD file; assume using this more
        often with snv file; returns tuple"""
        tbxFh =  self.__indel if indels else self.__snv
        caddChr = 'MT' if chrm == 'M' else xstr(chrm)
        return tbxFh.fetch(caddChr, start, end, parser=pysam.asTuple())

  
    def match(self, chrm, position, indels=True):
        """ match single position in CADD file; default indels to true
        b/c will be using this more with the indel file"""
        tbxFh = self.indel() if indels else self.snv()
        caddChr = 'MT' if chrm == 'M' else xstr(chrm)
        try:
            return tbxFh.fetch(caddChr, int(position) - 1, int(position), parser=pysam.asTuple())
        except ValueError as e: # happens sometimes on chrm M/MT
            warning("WARNING", e)
            return []
        
        
    def buffer_variant(self):
        chrm, position, ref, alt = self._current_variant.metaseq_id.split(':')
        isIndel = True if len(ref) > 1 or len(alt) > 1 else False
        self._current_variant
        if self.skip_existing():
            validationResult = self.is_cadd_annotated() # will update current_variant.record_primary_key if found
            if validationResult is None:
                self.log(("Variant", self._current_variant.metaseq_id, "not in DB, SKIPPING"), prefix="WARNING")
                self.increment_counter('skipped')
                self.increment_counter('not_matched')

            elif validationResult: # if not none, expect a boolean -- true already set
                self.log(("Variant", self._current_variant.metaseq_id, "/", self._current_variant.record_primary_key, "already updated, SKIPPING"), prefix="WARNING")
                self.increment_counter('skipped')

        else:
            matched = False
            for match in self.match(chrm, position, isIndel):
                matchedAlleles = [match[2], match[3]]

                if ref in matchedAlleles and alt in matchedAlleles:
                    evidence = {'CADD_raw_score': float(match[4]), 'CADD_phred': float(match[5])}
                    self.buffer_update_values(self._current_variant.record_primary_key, evidence)
                    # buffer_update_statement(self._current_variant.record_primary_key, evidence)  
                    self.increment_update_count(isIndel)
                    matched = True
                    
                    if self._debug:
                        self.log((self._current_variant.metaseq_id, "matched:", matchedAlleles, evidence), prefix="DEBUG")

                    break # no need to look for more matches

            if not matched:
                self.buffer_update_values(self._current_variant.record_primary_key, {})  
                # self.buffer_update_statement(self._current_variant.record_primary_key, {})  
                self.increment_counter('not_matched')
                if self._debug:
                    self.log((self._current_variant.metaseq_id, "not matched / loading placeholder"), prefix="DEBUG")
