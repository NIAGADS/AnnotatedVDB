#!/usr/bin/env python
from AnnotatedVDB.Util.adsp_consequence_parser import ConsequenceParser
from AnnotatedVDB.Util.conseq_group_enum import ConseqGroup
import argparse
import json

def run():
    fileName = args.fileName

    print(ConseqGroup.get_all_terms())

    cparser = ConsequenceParser(fileName, rankOnLoad=True, verbose=True, debug=True)
    for c in cparser.get_known_consequences():
        print(c)

    cparser.save_ranking_file()
    

    #print(json.dumps(cparser.get_rankings(), indent=4)) # pretty print

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fileName', required=True)

    args = parser.parse_args()

    run()
