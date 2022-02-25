#!/usr/bin/env python
""" split VCF files by chromosome """

import argparse
import gzip
import mmap
from os import path
from concurrent.futures import ProcessPoolExecutor
from GenomicsDBData.Util.utils import warning, xstr, print_dict, execute_cmd, create_dir, die, get_opener
from AnnotatedVDB.Util.enums import HumanChromosome as Human
from AnnotatedVDB.Util.parsers import ChromosomeMap


def create_fh_dict():
    """ create one fh for each chrm """
    chromosomes = [c.value for c in Human]
    header = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'] 
    for chrm in chromosomes:
        chrmS = xstr(chrm)
        ofName = path.join(args.outputDir, 'chr' + chrmS + ".vcf")
        fhDict[chrmS] = open(ofName, 'w')
        print('\t'.join(header), file=fhDict[chrmS])
    
        
def close_fh_dict():
    for key, fh in fhDict.items():
        fh.close()    

def run():
    """ main function """ 
    count = 0
    currentSequenceId = None
    with open(args.fileName, 'r') as fhandle:
        mappedFile = mmap.mmap(fhandle.fileno(), 0, prot=mmap.PROT_READ) # put file in swap
        with gzip.GzipFile(mode='r', fileobj=mappedFile) as gfh:
            for line in gfh:
                line = line.decode("utf-8").rstrip()  # in python 3, gzip.GzipFile seems to treat the strings a binary
                if line.startswith('#'):
                    continue
                count += 1
                values = line.split('\t')
                sequenceId = values[0]

                fhKey = sequenceId if chrmMap is None else chrmMap.get(sequenceId)    
                
                if currentSequenceId != sequenceId:
                    warning("New Sequence Found:", sequenceId, "- chr" + fhKey + ".vcf", prefix="WARNING")
                    currentSequenceId = sequenceId
                    
                print(line, file=fhDict[fhKey])
                if count % 5000000 == 0:
                    warning("Processed", '{:,}'.format(count), "lines", prefix="INFO")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fileName', help="full path to input file", required=True)
    parser.add_argument('-o', '--outputDir', help="output directory, will create if doesn't exist", required=True)
    parser.add_argument('-c', '--chromosomeMap', help="full path to chromosome map file")
    args = parser.parse_args()
    
    chrmMap = ChromosomeMap(args.chromosomeMap) if args.chromosomeMap else None
    create_dir(args.outputDir)
    
    fhDict = {} 
    create_fh_dict()
    run()
    close_fh_dict()
    