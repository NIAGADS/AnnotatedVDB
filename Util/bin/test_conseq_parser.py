#!/usr/bin/env python
from AnnotatedVDB.Util.adsp_consequence_parser import ConsequenceParser
from AnnotatedVDB.Util.conseq_group_enum import ConseqGroup
import argparse


def run():
    fileName = args.fileName

    print("#### Test ConseqGroup Enum -- get and print all terms")
    print(ConseqGroup.get_all_terms())
    print()
    
    print("#### Test load and rank consequence groups")
    print("Debug mode? " + str(args.debug))
    print("Verbose mode?" + str(args.verbose))
    cparser = ConsequenceParser(fileName, rankOnLoad=True, verbose=args.verbose, debug=args.debug)
    print()
    
    print("#### Test print ranked, known consequences -- consequences only")
    for c in cparser.get_known_consequences():
        print(c)
    print()
        
    combo = 'splice_acceptor_variant,splice_donor_variant,3_prime_UTR_variant,intron_variant'
    print("#### Test find ranking for a consequence combo: (expect 5)")
    rank = cparser.find_matching_consequence(combo.split(','))
    print("Match rank: " + str(rank))
    print()
            
    newConseqCombo = ['TFBS_amplification','TF_binding_site_variant']
    print("#### Test add term and re-rank: " + ','.join(newConseqCombo))
    print("Fail on missing? True")
    try:
        cparser.find_matching_consequence(newConseqCombo, failOnMissing=True)
    except Exception as err:
        print(str(err))
    print()
    
    print("#### Fail on missing? False")
    rank = cparser.find_matching_consequence(newConseqCombo)
    print("Match rank: " + str(rank))
    print()
    
    print("#### Test saving updated ranking to a file/file name generated by timestamp")
    cparser.save_ranking_file()
    
    print("Done with test")
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('-f', '--fileName', required=True)

    args = parser.parse_args()

    run()
