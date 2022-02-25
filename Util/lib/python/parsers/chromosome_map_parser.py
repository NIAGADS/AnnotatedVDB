"""! @brief Chromosome Map Parser"""

##
# @package parsers
# @file chromosome_map_parser.py
#
# @brief  Parse and Access Chromosome Mappings
#
# @section chromosome_map_parser Description
# parses mappings of third party chromosome source ids (e.g., refseq) to chromosome number
# may also include chromosome length
# - for now assumes the tab-delim with at least the following columns:
# > source_id       chromosome     length
#
# @section todo_chromosome_map_parser TODO
# 
#
# @section libraries_chromosome_map_parser Libraries/Modules
# - cvs - reading and parsing the file
#
# @section author_chromosome_map_parser Author(s)
# - Created by Emily Greenfest-Allen (fossilfriend) 2022

import csv
from GenomicsDBData.Util.utils import warning, xstr

class ChromosomeMap(object):
    """! Generator for chromosome map parser/object
    parses mappings of third party chromosome sequence ids (e.g., refseq) to chromosome number
    may also include chromosome length
    - for now assumes the tab-delim with at least the following columns:
    > source_id       chromosome     length      
    """
    
    def __init__(self, fileName, verbose=False, debug=False):
        """! ChromosomeMap base class initializer
            @param fileName             full path to chromosome mapping file
            @param verbose              verbose output flag
            @param debug                debug flag
            @return                     An instance of a ChromosomeMap with initialized mapping dict
        """
        self.__verbose = verbose
        self.__debug = debug
        self.__fileName = fileName
        self.__map = {}
        self.__parse_mapping()
        
        
    def __parse_mapping(self):
        """! parse chromosome map
        """

        if self.__verbose:
            warning("Loading chromosome map from:", self.__fileName)

        with open(self.__fileName, 'r') as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for row in reader:
                # source_id	chromosome	chromosome_order_num	length
                key = row['source_id']
                value = row['chromosome'].replace('chr', '')
                self.__map[key] = value
                        
    
    def chromosome_map(self):
        """! @returns                 the chromosome map
        """
        return self.__map
        
        
    def get_sequence_id(self, chrmNum):
        """! given a chromosome number, tries to find matching sequence id

            @param chrmNum             chromosome number to match
            @returns                   matching sequence id, None if not found
        """
        for sequenceId, cn in self.__map.items():
            if cn == chrmNum or cn == 'chr' + xstr(chrmNum):
                return sequenceId
            
        return None

            
    def get(self, sequenceId):
        """! return chromosome number mapped to the provided sequence ID

            @param sequenceId           sequence ID to look up
            @returns                    matching chrm number, will fail on error
        """
        # want to raise AttributeError if not in the map, so not checking
        return self.__map[sequenceId]
        
            