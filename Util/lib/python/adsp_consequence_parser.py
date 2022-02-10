# TODO: save updated conseques
# Check Alpha/numeric sort from Nick's code

'''
Some python functions for ranking consequence combinations
Edited for increased readability

Modified from code by Nick Wheeler, Case Western (W. Bush Group)

modifed by EGA/fossilfriend (2021-2022) as follows:
  * wrap in class
  * moved file/load parsing & consequence matching / fetching from vep_parser.py
  * added save option
  * removed unnecessary string parsing which accounts for by JSON & parsing UTF-8 encoding on file reads
  * introduced readable variable names to follow conventions of GenomicsDB coding
  * replace original code dependency on pandas w/set operations & OrderedDict

NOTE:
  `terms`: a single consequence combination of one or more terms; in list form
  `conseq`: a single consequence combination, in comma separated string form
  `conseqs`: a list of consequence combinations

'''

# pylint: disable=line-too-long,invalid-name,no-self-use
from __future__ import print_function

import csv
import json
from collections import Counter, OrderedDict
from datetime import date

from GenomicsDBData.Util.utils import warning, to_numeric, die, int_to_alpha, verify_path
import GenomicsDBData.Util.list_utils as lu
from AnnotatedVDB.Util.conseq_group_enum import ConseqGroup


class ConsequenceParser(object):
    ''' class to organize utils for parsing and re-ranking ADSP ranked consequences of VEP output '''

    def __init__(self, rankingFileName, saveOnAddConseq=False, rankOnLoad=False, verbose=False, debug=False):
        self._verbose = verbose
        self._debug = debug
        self._rankingFileName = rankingFileName
        self._consequenceRankings = self._parse_ranking_file(rankingFileName)
        self._addedConsequences = []
        self._saveOnAddConsequence = saveOnAddConseq
        
        if rankOnLoad: # re-rank file on load (e.g., first time using a ranking file)
            if self._verbose:
                warning("INFO:", "rankOnLoad = True")
            self._update_rankings()
            
        self._matchedConseqTerms = {} # for already matched/to speed up lookups


    def save_ranking_file(self, fileName=None):
        ''' save ranking if file.  
          if no file name is provided, will upate original file name with today's date
          '''
        header = lu.qw('consequence rank')
        if fileName is None:
            fileName = self._rankingFileName.split('.')[0] + "_" \
              + date.today().strftime("%m-%d-%Y") + ".txt"
              
        if verify_path(fileName): # file already exists / add versioning for each newly added conseq
            fileName = self._rankingFileName.split('.')[0] + "_v" + self.get_new_conseq_count() + ".txt"

        with open(fileName, 'w') as ofh:
            print('\t'.join(header), file=ofh)
            for conseq, rank in self._consequenceRankings.items():
                print(conseq, rank, sep='\t', file=ofh)
            

    def _parse_ranking_file(self, fileName):
        ''' parse ranking file and save as dictionary lookup '''
        if self._verbose:
            warning("Parsing ranking file: ", fileName)

        result = OrderedDict()
        rank = 1
        with open(fileName, 'r') as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                conseq = lu.alphabetize_string_list(row['consequence']) # ensures unique keys
                if 'rank' in row:
                    result[conseq] = to_numeric(row['rank'])
                else: # assume load order is rank order
                    result[conseq] = rank
                    rank = rank + 1

                
        return result


    # =========== accessors ==================

    def get_new_conseq_count(self):
        '''return count of newly added consequences'''
        return len(self._addedConsequences)
    

    def new_consequences_added(self):
        ''' flag = true if new consequences were added '''
        return len(self._addedConsequences) > 0

    
    def get_added_consequences(self):
        ''' return list of new consequences '''
        return self._addedConsequences
        
        
    def get_rankings(self):
        ''' return consequence rankings '''
        return self._consequenceRankings


    def get_consequence_rank(self, conseq, failOnError=False):
        ''' return value from consequence rank map for the specified
        consequence '''
        if conseq in self._consequenceRankings:
            return self._consequenceRankings[conseq]
        else:
            if failOnError:
                raise IndexError('Consequence ' + conseq + ' not found in ADSP rankings.')
            return None
        

    def find_matching_consequence(self, terms, failOnMissing=False):
        ''' match list of consequences against those in the
        rankings and return ranking info associated with match

        attempt to integrate into the rankings if not found '''
        
        if len(terms) == 0:
            return self.get_consequence_rank(terms[0])

        conseqKey = '.'.join(terms)
        # storing matches b/c this step is slow (checking equivalent lists etc)
        if conseqKey not in self._matchedConseqTerms:
            match = None
            for conseqStr in self._consequenceRankings:
                conseqList = conseqStr.split(',')
                if lu.is_equivalent_list(terms, conseqList):
                    match = self.get_consequence_rank(conseqStr)
                    break
                
            if match is None:
                # no match found
                if failOnMissing:
                    raise IndexError('Consequence combination ' + ','.join(terms) + ' not found in ADSP rankings.')
                else: # add term & make recursive cal
                    if self._verbose:
                        warning('Consequence combination ' + ','.join(terms) + ' not found in ADSP rankings.')
                    self._update_rankings(terms)
                    return self.find_matching_consequence(terms)
                
            self._matchedConseqTerms[conseqKey] = match

        return self._matchedConseqTerms[conseqKey]

        
    # =========== reranking functions  ==================

    def get_known_consequences(self):
        ''' extract the known consequences/combos 
            replaces `ret_sort_combos`

            returns as set to try and prevent duplicates

            b/c these are keys to a hash, they are unique
        '''
        return list(self._consequenceRankings.keys()) # convert from odict_keys object


    def _add_consequence(self, terms):
        ''' add new consequence (combination) to the list of known consequences '''

        # extract list of known consequences & add the new one
        referenceConseqs = self.get_known_consequences()
        conseqStr = lu.alphabetize_string_list(terms)
        
        if conseqStr in referenceConseqs:
            raise IndexError('Attempted to add consequence combination ' \
                                + conseqStr + ', but already in ADSP rankings.')
        
        referenceConseqs.append(conseqStr)
        self._addedConsequences.append(conseqStr)

        return referenceConseqs

        
    def _update_rankings(self, terms=None):
        '''' update rankings, adding in new term combination if specified,
        where terms is a list of one or more consequences
        
        ranks are applied according to the following logic:

        1. Split all consequence combos into 4 groups:
           - GRP1 - nmd: contains `NMD_transcript_variant` - ConseqGroup.NMD
           - GRP2 - nct: contains `non_coding_transcript_variant` - ConseqGroup.NON_CODING_TRANSCRIPT
           - GRP3 - low: contains ONLY consequences in ConseqGroup.LOW_IMPACT
           - GRP4 - high: includes at least one ConseqGroup.HIGH_IMPACT / should not overlap with groups 1&2

        2. Process GRPS in descending order

        NOTE: iterating over the ConseqGroup enum will allow processing of GRPS in expected order

        TODO: fill in documentation from notes from N. Wheeler

        '''

        if self._verbose:
            warning("Updating consequence rankings")
            
        conseqs = self._add_consequence(terms) if terms is not None else self.get_known_consequences()

        sortedConseqs = []
        for grp in ConseqGroup:
            if self._verbose:
                warning("Ranking", grp.name, "consequences.")
            requireSubset = True if grp.name == 'LOW_IMPACT' else False
            members = grp.get_group_members(conseqs, requireSubset) # extract conseqs belonging to current group
            if self._verbose:
                warning("Found " + str(len(members)) + ":", members)
            if len(members) > 0:
                sortedConseqs += self._sort_consequences(members, grp)
                
        # convert to dict & update
        if self._debug:
            warning("FINAL SORTED CONSEQUENCES", sortedConseqs)
            
        self._consequenceRankings = lu.list_to_indexed_dict(sortedConseqs)

        if terms is not None and self._saveOnAddConsequence:
            if self._verbose:
                warning("Added new consequence `" + ','.join(terms) + "`. Saving version:", self.get_new_conseq_count())
            self.save_ranking_file()

            
    def _sort_consequences(self, conseqs, conseqGrp):
        ''' 
        sort an input list of consequence terms (a consequence combination) according 
        to a numerially indexed dictionary of their rankings
        
        dictionary of their rankings based on the 
        conseqGrp (type ConseqGroup enum) toDict()`

        returns a list of consequence terms sorted as follows:
        
        '''

        # if LOW_IMPACT (GRP 3) use LOW_IMPACT,
        # else use HIGH_IMPACT
        # ranking dicts are retrieved here so that the calculation
        
        grpRankingDict = ConseqGroup.HIGH_IMPACT.toDict() \
          if conseqGrp.name != 'LOW_IMPACT' \
          else conseqGrp.toDict()

        completeRankingDict = ConseqGroup.get_complete_indexed_dict() # needed for non-exclusive grps
        if self._debug:
            warning("GRP Dict:", json.dumps(grpRankingDict, indent=4))
            warning("COMPLETE Dict:", json.dumps(completeRankingDict, indent=4))

        indexedConseqs = []
        for c in conseqs:
            indexedConseqs.append(self._calculate_ranking_indexes(c, grpRankingDict, completeRankingDict))

        if self._debug:
            warning("Indexed", conseqGrp.name, "consequences:", indexedConseqs)

        sortedConseqs = indexedConseqs
        sortedConseqs.sort(key=lambda x: x[0]) # sorts normally by alphabetical order            
        sortedConseqs.sort(key=lambda x:len(x[0]), reverse=True) # sorts by descending length
        sortedConseqs.sort(key=lambda x: x[0][0])# sorts by first character

        if self._debug:
            warning("Sorted", conseqGrp.name, "consequences:", sortedConseqs)

        return [','.join(sc[1]) for sc in sortedConseqs] # (alpha index, term) tuples; return term as str


    def _calculate_ranking_indexes(self, conseq, grpDict, refDict):
        '''  return tuple of
             (alphabetic representation, internally sorted consequence combo as a list)

        given a set of input consequences 
        and a numerically indexed dictionary of their ranking based on grp and reference dict
        for non-group members
    
        the alphabetic representation allows to rank on order & sum of term ranks

        modified from N. Wheeler's code as follows:
        * if terms are not present in the ranking dict in the original code,
          they are alphabetized & appended to the ranking dict, with new indexes
          --> here, a ranking dictionary containing all terms is passed for evaluating
          sets with non-exclusive group membership, and the alphabetized non-member terms
          are ranked according that dict
        '''

        terms = conseq.split(',')
        if self._debug:
            warning(terms)
        memberTerms = [c for c in terms if c in grpDict]
        nonMemberTerms = [c for c in terms if c not in grpDict]
        if self._debug:
            warning("member terms:", memberTerms)
            warning("nonmember terms:", nonMemberTerms)
            
        indexes = self._get_consequence_ranks(memberTerms, grpDict)
        if len(nonMemberTerms) > 0:
            indexes += self._get_consequence_ranks(nonMemberTerms, refDict)

        if self._debug:
            warning("indexes:", indexes)
            
        # get alphabetic representation
        alphaIndexes = [int_to_alpha(x) for x in indexes]
        if self._debug:
            warning("alpha representation:", alphaIndexes)
        alphaIndexes.sort()

        # sort consequences in the combination by their ranking indexes
    
        indexedConseq = self._internal_consequence_sort(memberTerms + nonMemberTerms, indexes, returnStr=False)
        if self._debug:
            warning("indexed", indexedConseq)
            
        return (''.join(alphaIndexes), indexedConseq)


    def _internal_consequence_sort(self, terms, rankings, returnStr = False):
        ''' sort a consequence list by a list of rankings
        rankings are based on ConseqGroup indexes / hence 'internal'
        if returnStr = true, return comma separated string, else return sorted list
        '''
        rankingDict = dict(zip(terms, rankings))
        sortedDict =  OrderedDict(sorted(rankingDict.items(), key=lambda kv: kv[1])) # sort by values
        if returnStr:
            return ','.join(list(sortedDict.keys()))
        else:
            return list(sortedDict.keys())
    

    def _get_consequence_ranks(self, terms, rankingDict):
        '''
        retrieve a list of consequence rankings for each term in 
        the consequence combination according to a numerially indexed 
        dictionary of consequences in the assigned group

        # NOTE: all terms should be in the rankingDict by this time / errors should
        # have been caught when list of consequences was split into groups
        '''
        return [rankingDict[c] for c in terms]
