#!/usr/bin/env python
# pylint: disable=unused-import,invalid-name
"""
Loads AnnotatedVDB from JSON output from running VEP on dbSNP
"""
from __future__ import print_function

import argparse
import threading
import csv
import gzip
import mmap
import json

from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor
from CBILDataCommon.Util.utils import warning, xstr, die, qw, execute_cmd
from CBILDataCommon.Util.postgres_dbi import Database

CHROMOSOMES = qw('1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT');

def parse_ranking():
    ''' parse ranking file and save as dictionary lookup '''
    result = OrderedDict()
    with open(args.rankingFile, 'r') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            conseq = row['consequence']
            result[conseq]['rank'] = row['adsp_ranking']
            result[conseq]['impact'] = row['adsp_impact']

    return result

def process_vcf_input_string(inputStr):
    ''' processes the VCF input string and return map '''
    fields = qw('chrom pos id ref alt qual filter info')
    values = inputStr.split('\t')
    result = dict(zip(fields, values))
    
    # now unpack the info field and save as its own
    values = result['info']
    
    
    # VCF 
    # CHROM POS     ID        REF ALT    QUAL FILTER INFO  
    #X\t605409\trs780063150\tC\tA\t.\t.\tRS=780063150;RSPOS=605409;dbSNPBuildID=144;SSR=0;SAO=0;VP=0x05000088000d000026000100;GENEINFO=SHOX:6473;WGT=1;VC=SNV;U3;INT;CFL;ASP;KGPhase3;CAF=0.9996,0.0003994;COMMON=0;TOPMED=0.99999203618756371,0.00000796381243628

def load_annotation(chromosome):
    ''' parse over a JSON file, extract position, frequencies,
    ids, and ADSP-ranked most severe consequence; bulk load using COPY '''
    fname = chromosome + '.json'
    opener = open
    if args.isGzipped:
        fname = fname + '.gz'
        opener = gzip.open

    with opener(fname, 'r') as fh:
        mapFile = mmap.mmap(fh.fileno(), 0)
        for line in iter(mapFile.readline, ""): 
            # NOTE in python 3 the sentinal parameter needs to be b"" (type bytes)
            annotation = json.loads(line.rstrip())
            
            ids = []
            if args.dbSNP:
                ids = variants_from_dbsnp_input_string(annotation['input'])
            else:
                variant = { 'metaseq_id': annotation['input'],
                            'ref_snp_id': extract_ref_snp_id(annotation) }
                ids.append(variant)

            isMultiAllelic = len(ids > 1)

            if isMultiAllelic:
                frequencies = extract_multiallelic_frequencies(annotation)
            else:
                frequencies = extract_frequencies(annotation)

            # rank consequence





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='load AnnotatedDB from VEP JSON output')
    parser.add_argument('-d', '--dir', help="directory containing VEP results", required=True)
    parser.add_argument('-r', '--rankingFile', help="full path to ADSP VEP consequence ranking file", required=True)
    parser.add_argument('--dbSNP', action='store_true', help="loading dbSNP?")
    parser.add_argument('--maxWorkers', help="max number of threads", type=int, default=5)
    parser.add_argument('--commit', action='store_true', help="run in commit mode", required=False)
    parser.add_argument('--gusConfigFile', '--full path to gus config file, else assumes $GUS_HOME/config/gus.config')
    parser.add_argument('--isGzipped', help="flag to indicate that input files are gzipped")
    parser.add_argument('--test', help="load 100 rows as test", action='store_true')
    parser.add_argument('-c', '--chr', help="comma separated list of chrs, e.g, 1,5,7,M; if not specified will run all")
    args = parser.parse_args()

    chrList = CHROMOSOMES
    if args.chr:
        chrList = args.chr.split(',')

    consequenceRanks = parse_ranking()
    die(consequenceRanks)

    with ThreadPoolExecutor(max_workers=args.maxWorkers) as executor:
        for c in chrList:          
            warning("Create and start thread for chromosome:", xstr(c))
            executor.submit(load_annotation, c)        
