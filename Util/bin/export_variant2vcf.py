#!/usr/bin/env python3 
#pylint: disable=invalid-name

""" 
Export AnnotatedVDB.Variant to VCF files to 
support bulk processing and updating of table 
"""

import argparse
import os.path as path
import random

from concurrent.futures import ProcessPoolExecutor, as_completed
from psycopg2 import DatabaseError

from niagads.utils.string import xstr
from niagads.utils.dict import print_dict
from niagads.utils.sys import warning, print_args, execute_cmd
from niagads.utils.list import qw
from niagads.utils.reg_ex import matches
from niagads.db.postgres import Database

from AnnotatedVDB.Util.enums import HumanChromosome as Human

SELECT_SQL = "SELECT record_primary_key, metaseq_id, row_algorithm_id FROM AnnotatedVDB.Variant WHERE chromosome = %s"
VARIANTS_PER_FILE = 10000000
ITERATION_SIZE = 500000
VCF_HEADER = qw('#CHRM POS ID REF ALT QUAL FILTER INFO', returnTuple=False)
INVALID_ALLELES = 'I|R|D|N' # there are some "N's" in the DB, but the appear to all come from dbSNP directly

def print_file(variants, directory, chromosome, count):
    vcfFile = path.join(directory, chromosome + '_' + xstr(count) + '.vcf')
    with open(vcfFile, 'w') as fh:
        print(*VCF_HEADER, sep='\t', file=fh) 
        for recordPk, metaseqId in variants.items():
            chr, pos, ref, alt = metaseqId.split(':')
            print(chr, pos, recordPk, ref, alt, '.', '.', '.', sep='\t', file=fh)


def run(directory:str, chromosome:str):
    logFileName = path.join(directory, chromosome + ".log")
    invalidFileName = path.join(directory, chromosome + "_invalid.txt")
    # duplicateFileName = path.join(directory, chromosome + "_duplicates.txt")
    warning("Starting " + chromosome + "; logging to " + logFileName)

    duplicates = {}
    try:
        database = Database(args.gusConfigFile)
        database.connect()
        cname = 'select_' + chromosome
        
        fileCount = 1
        with database.named_cursor(cname, cursorFactory="RealDictCursor") as selectCursor, \
            open(logFileName, 'w') as lfh, \
            open(invalidFileName, 'w') as ifh:
            # open(duplicateFileName, 'w') as dfh, \

            warning("Processing " + chromosome, file=lfh, flush=True)
            warning("Executing query: " + SELECT_SQL.replace('%s', "'" + chromosome + "'"), file=lfh, flush=True)
            
            selectCursor.itersize = 500000 
            selectCursor.execute(SELECT_SQL, [chromosome])
            validCount = 0
            variants = {}
            for record in selectCursor:     
                primaryKey = record['record_primary_key']
                metaseqId = record['metaseq_id']
                rowAlgId = record['row_algorithm_id']
                
                """ it appears all (hopefully) duplicates are invalid variants or from an early NHGRI GWAS Catalog Load
                if primaryKey in variants or duplicates.get(primaryKey) != True: # not in the hash or False
                    print(primaryKey, rowAlgId, sep='\t', file=dfh, flush=True)
                    duplicates[primaryKey] = True
                    continue
                """
                
                if matches(INVALID_ALLELES, metaseqId):
                    print(primaryKey, rowAlgId, sep='\t', file=ifh, flush=True)
                    continue
                else:
                    validCount += 1
                    variants[primaryKey] = metaseqId
                
                if validCount % ITERATION_SIZE == 0:
                    warning("Fetched " + xstr(validCount) + " valid records from " + chromosome + ".", file=lfh, flush=True)
                    
                if validCount % VARIANTS_PER_FILE == 0:
                    warning("Fetched " + xstr(validCount) + " records from " + chromosome + ". Writing to file " + xstr(fileCount), file=lfh, flush=True)
                    print_file(variants, directory, chromosome, fileCount)
                    
                    # del duplicates
                    # duplicates = { k: False for k in variants.keys()}   

                    variants = {}
                    fileCount += 1
            
            # residuals
            warning("Fetched " + xstr(validCount) + " valid records from " + chromosome + ". Writing residuals to file " + xstr(fileCount), file=lfh, flush=True)
            print_file(variants, directory, chromosome, fileCount)
                
    except Exception as err:
        raise(err)
    finally:
        database.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="export AnnotatedVDB.Variant to VCF files", allow_abbrev=False)
    parser.add_argument('-o', '--outputDir', required=True,
                        help="output directory")
    parser.add_argument('--chr', default='all',
                        help="chromosome; comma separated list of one or more chromosomes or 'all' or 'autosome'; ignored if --file is specified")
    parser.add_argument('--maxWorkers', type=int, default=10)
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--veryVerbose', action='store_true')
    parser.add_argument('--gusConfigFile',
                        help="GUS config file. If not provided, assumes default: $GUS_HOME/conf/gus.config")
    args = parser.parse_args()
    
    
    chrList = None
    if args.chr not in ['all', 'autosome']:
        chrList = args.chr.split(',')
    elif args.chr == 'autosome':
        chrList = [c.value for c in Human if c.value not in ['X', 'Y', 'M', 'MT']] 
    else: 
        chrList = [c.value for c in Human]
        
    random.shuffle(chrList) # so that not all large chrms are done at once if all is selected
    with ProcessPoolExecutor(args.maxWorkers) as executor:
        futureUpdate = {executor.submit(run, directory=args.outputDir,  chromosome='chr' + xstr(c)) : c for c in chrList}
        for future in as_completed(futureUpdate): # this should allow catching errors 
            try:
                future.result()
            except Exception as err:
                raise(err)       