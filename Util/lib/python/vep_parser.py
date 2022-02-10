'''
utils for parsing VEP JSON output
'''
# pylint: disable=line-too-long,invalid-name

from GenomicsDBData.Util.utils import warning, to_numeric, die
from GenomicsDBData.Util.list_utils import qw
from GenomicsDBData.Util.auto_viv_dict import AutoVivificationDict
from AnnotatedVDB.Util.adsp_consequence_parser import ConsequenceParser

from operator import itemgetter
from datetime import datetime

CONSEQUENCE_TYPES = qw('transcript regulatory_feature motif_feature intergenic')

CODING_CONSEQUENCES= qw('synonymous_variant missense_variant inframe_insertion inframe_deletion stop_gained stop_lost stop_retained_variant start_lost frameshift_variant coding_sequence_variant')

def is_coding_consequence(conseqStr):
    ''' check term against list of coding consequences and return 
    True if found '''
    matches = [value for value in conseqStr.split(',') 
                   if value in CODING_CONSEQUENCES]

    return len(matches) > 0


class VepJsonParser(object):
    ''' class to organize utils for parsing VEP JSON output '''

    def __init__(self, rankingFileName, verbose=True):
        self._verbose = verbose
        self._consequenceParser = ConsequenceParser(rankingFileName)
        # self._consequenceRankings = self.parse_ranking_file(rankingFileName)
        self._annotation = None


    def assign_adsp_consequence_rank(self, conseq):
        ''' assemble consequence terms into ordered comma delim string
        and lookup in ranking table; add rank to consequence annotation;
        backup old impact information
        '''

        terms = conseq['consequence_terms']
        matchingConseq = self.find_matching_consequence(terms) 

        conseq['vep_impact'] =  conseq['impact']
        conseq['impact'] = matchingConseq['impact']
        conseq['rank'] = matchingConseq['rank']
        conseq['consequence_is_coding'] = matchingConseq['coding']

        return conseq
        

    def __verify_annotation(self):
        ''' check that annotation is set '''
        assert self._annotation is not None, \
          "DEBUG - must set value of _annotation in the VEP parser to access it"



    def adsp_rank_and_sort_consequences(self):
        ''' applies ADSP ranking to consequence, re orders the array and
        and saves newly ordered list
        '''
        for ctype in CONSEQUENCE_TYPES:
            rankedConseqs = self.__adsp_rank_consequences(ctype)
            if rankedConseqs is not None:
                self.set(ctype + '_consequences', rankedConseqs)

            
    # =========== modifiers ==================
    def set_annotation(self, annotation):
        ''' set the annotation json '''
        self._annotation = annotation

    
    def set(self, key, value):
        ''' set a value to annotation json '''
        self.__verify_annotation()
        self._annotation[key] = value


    # =========== accessors ==================
   
    def get_conseq_rank(self, conseq):
        ''' return value from consequence rank map for the specified
        consequence '''
        if conseq in self._consequenceRankings:
            return self._consequenceRankings[conseq]
        else:
            raise IndexError('Consequence ' + conseq + ' not found in ranking file.')


    def _adsp_rank_consequences(self, conseqType):
        ''' extract consequences and apply ranking and sort,
        convert from list to list of dicts, keyed on allele
        to ensure consequences are sorted per allele'''

        result = None
        localConseqMap = {} # so that we don't have to look up repeated consequences
        consequences = self.get(conseqType + '_consequences')
        if consequences is None:
            return None

        result = {}
        for index, conseq in enumerate(consequences):
            # need to build the hash pulling the key out of variant allele &
            # appending to the list
            va = conseq['variant_allele']
            conseq['vep_consequence_order_num'] = index

            if va in result:
                result[va].append(self.assign_adsp_consequence_rank(conseq))
            else:
                result[va] = [self.assign_adsp_consequence_rank(conseq)]

        for va in result: # sort by rank within each variant allele
            # result[va] = sorted(result[va], key = lambda x: (x['rank'], x['vep_consequence_order_num']))
            # itemgetter supposed to be orders of magnitude faster than lambda
            result[va] = sorted(result[va], key=itemgetter('rank', 'vep_consequence_order_num'))

        return result


    def get_frequencies(self, matchingVariantId=None):
        ''' extract frequencies from colocated_variants section
        for dbSNP VEP run -- sometimes VEP matches to variants with
        rsId different from current rsId; use matchingVariantId to ensure
        extracting the correct frequency (e.g. overlapping indels & snvs)
        '''

        self.__verify_annotation()
        if 'colocated_variants' not in self._annotation:
            return None

        frequencies = None
        cv = self._annotation['colocated_variants']

        if len(cv) > 1:
            # return first non-cosmic mutation that matches the expected variant id (if supplied)
            fCount = 0
            for covar in cv:
                if covar['allele_string'] != 'COSMIC_MUTATION':
                    if 'frequencies' in covar:
                        if matchingVariantId is not None:
                            if covar['id'] == matchingVariantId: 
                                frequencies = self.__extract_frequencies(covar)
                        else:
                            frequencies = self.__extract_frequencies(covar)
                            fCount = fCount + 1
            if fCount > 1: 
                # based on experience, when this happens, involves multiple refsnps mapped to location,
                # so all frequencies should be equal
                # let's just print a warning
                inputVariant = self._annotation['input'].replace('\t', ' ')
                warning("WARNING", "INDEL " + inputVariant + "mapped to multiple refSNPs/frequencies based on location not alleles")
            # else:
            return frequencies # which may be None
        
        elif 'frequencies' in cv[0]:
            frequencies = self.__extract_frequencies(cv[0])

        return frequencies


    def __extract_frequencies(self, covar):
        ''' extract frequencies and update minor allele to include 
        1000 Genomes global allele frequency if present (stored in
        colocated_variant field not frequencies array)
        '''

        frequencies = {}
        if 'minor_allele' in covar:
            frequencies['minor_allele'] = covar['minor_allele']
            if 'minor_allele_freq' in covar:
                frequencies['minor_allele_freq'] = covar['minor_allele_freq']

        frequencies['values'] = self.__group_frequencies_by_source(covar['frequencies'])
        return frequencies
            

    def __group_frequencies_by_source(self, frequencies):
        ''' group 1000Genomes and gnomAD frequencies '''

        if frequencies is None:
            return None
        
        result = AutoVivificationDict() # needed to later add in minor allele freqs
        espKeys = ['aa', 'ea']
        for allele in frequencies:
            gnomad = [{key: value} for key, value in frequencies[allele].items() if 'gnomad' in key ]
            esp = [{key: value} for key, value in frequencies[allele].items() if key in espKeys]
            genomes = [{key: value} for key, value in frequencies[allele].items() if 'gnomad' not in key and key not in espKeys]
            if bool(gnomad):
                result[allele]['gnomAD'] = gnomad
            if bool(genomes):
                result[allele]['1000Genomes'] = genomes
            if bool(esp):
                result[allele]['ESP'] = esp

        return result


    def __get_consequences(self, key):
        ''' special getter for consequences b/c fields may be missing; don't want
        to throw error '''
        if key in self._annotation:
            return self._annotation[key]
        else:
            return None


    def get(self, key):
        ''' get the annotation value associated with the key '''
        self.__verify_annotation()

        if key == 'frequencies':
            return self.get_frequencies()
        if 'consequences' in key:
            return self.__get_consequences(key)
        else:
            return self._annotation[key]


    def get_annotation(self):
        ''' return updated annotation '''
        return self._annotation


