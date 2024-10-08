#!/usr/bin/env python3
# pylint: disable=invalid-name,not-an-iterable,unused-import,too-many-locals
"""
Loads AnnotatedVDB variants from from VCF File
Loads each chromosome in parallel
"""

from __future__ import print_function

import argparse
import gzip
import mmap
import json
import sys
import csv

from datetime import datetime
from os import path
from sys import stdout
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from psycopg2 import DatabaseError

from GenomicsDBData.Util.utils import xstr, warning, print_dict, print_args, die, get_opener
from GenomicsDBData.Util.postgres_dbi import Database, raise_pg_exception

from AnnotatedVDB.Util.loaders import VCFVariantLoader 
from AnnotatedVDB.Util.parsers import ChromosomeMap
from AnnotatedVDB.Util.enums import HumanChromosome as Human


def initialize_loader(logFilePrefix):
    """! initialize loader """

    lfn = xstr(logFilePrefix)
    if args.resumeAfter:
        lfn += '_resume_' + args.resumeAfter 
    if args.failAt:
        lfn += '_fail_' + args.failAt.replace(':', '-')
    lfn += '.log'
    warning("Logging to", lfn)
    
    try:
        loader = VCFVariantLoader(args.datasource, logFileName=lfn, verbose=args.verbose, debug=args.debug)
        
        if args.verbose:
            loader.log(("Parameters:", print_dict(vars(args), pretty=True)), prefix="INFO")

        loader.set_algorithm_invocation('load_vcf_result', xstr(logFilePrefix) + '|' + print_args(args, False))
        loader.log('Algorithm Invocation Id = ' + xstr(loader.alg_invocation_id()), prefix="INFO")
        
        loader.initialize_pk_generator(args.genomeBuild, args.seqrepoProxyPath)
        loader.initialize_bin_indexer(args.gusConfigFile)
        loader.initialize_copy_sql() # use default copy fields
        loader.set_chromosome_map(chrmMap)
        
        loader.set_skip_existing(args.skipExisting, args.gusConfigFile) # initialize validator db connection
        
        if args.failAt:
            loader.set_fail_at_variant(args.failAt)
            loader.log("Fail at variant: " + loader.fail_at_variant(), prefix="INFO")
        
    except Exception as err:
        warning("ERROR", "Problem initializing Loader")
        raise(err)
        
    return loader


def load_annotation(fileName, logFilePrefix):
    """! parse over a VCF file; bulk load using COPY """
 
    loader = initialize_loader(logFilePrefix)
    loader.log('Parsing ' + fileName, prefix="INFO")
    loader.log('Writing metaseq_id -> primary_key mapping to ' + fileName + '.mapping', prefix='INFO')
    
    resume = args.resumeAfter is None # false if need to skip lines
    if not resume:
        loader.log(("--resumeAfter flag specified; Finding skip until point", args.resumeAfter), prefix="INFO")
        loader.set_resume_after_variant(args.resumeAfter)

    testComplete = False # flag if test is complete
    try: 
        database = Database(args.gusConfigFile)
        database.connect()
        opener = get_opener(args.fileName)
        with opener(fileName, 'r') as fhandle, database.cursor() as cursor, \
            open(fileName + ".mapping", 'w') as mfh:
            loader.set_cursor(cursor)
            mappedFile = mmap.mmap(fhandle.fileno(), 0, prot=mmap.PROT_READ) # put file in swap
            lineCount = 0
            for line in iter(mappedFile.readline, b""):
                line = line.decode("utf-8")  # in python 3, mapped file is a  binary IO object
                if line.startswith("#"): # skip comments
                    continue    
                
                if (lineCount == 0 and not loader.resume_load()) \
                    or (lineCount % args.commitAfter == 0 and loader.resume_load()): 
                    if args.debug:
                        loader.log('Processing new copy object', prefix="DEBUG")
                    tstart = datetime.now()
                        
                primaryKeyMapping = loader.parse_variant(line.rstrip())
                for metaseqId, pk in primaryKeyMapping.items():
                    print(metaseqId, pk, sep='\t', file=mfh, flush=True)
               
                lineCount += 1
                    
                if not loader.resume_load():
                    if lineCount % args.logAfter == 0:
                        loader.log(('{:,}'.format(lineCount), "lines"), 
                                    prefix="SKIPPED")
                    continue
                    
                if loader.resume_load() != resume: # then you are at the resume cutoff
                    resume = True
                    loader.log(('{:,}'.format(lineCount), "lines"), 
                                prefix="SKIPPED")
                    continue
                
                if lineCount % args.logAfter == 0 \
                    and lineCount % args.commitAfter != 0:
                    loader.log((lineCount, "lines (", 
                                loader.get_count('variant'), " variants)"), 
                                prefix="PARSED")

                if lineCount % args.commitAfter == 0:
                    if args.debug:
                        tendw = datetime.now()
                        message = 'Copy object prepared in ' + str(tendw - tstart) + '; ' + \
                            str(loader.copy_buffer(sizeOnly=True)) + ' bytes; transfering to database'
                        loader.log(message, prefix="DEBUG") 
                    
                    loader.load_variants()
                    if args.datasource.lower() == 'adsp':
                        loader.update_variants()

                    message = '{:,}'.format(loader.get_count('variant')) + " variants"
                    messagePrefix = "COMMITTED"
                    if args.commit:
                        database.commit()
                    else:
                        database.rollback()
                        messagePrefix = "LOADED"
                        message += " -- rolling back"

                    if lineCount % args.logAfter == 0:
                        message += "; up to = " + loader.get_current_variant_id()
                        loader.log(message, prefix=messagePrefix)
                        
                        if loader.get_count('update') > 0:
                            message = '{:,}'.format(loader.get_count('update')) + " variants"
                            loader.log(message, prefix="UPDATED")
                        
                        if loader.get_count('duplicates') > 0:
                            message = '{:,}'.format(loader.get_count('duplicates')) + " variants"    
                            loader.log(message, prefix="SKIPPED")

                        if args.debug:
                            tend = datetime.now()
                            loader.log('Database copy time: ' + str(tend - tendw), prefix="DEBUG")
                            loader.log('        Total time: ' + str(tend - tstart), prefix="DEBUG")

                    if args.test:
                        break

            mappedFile.close()
            
            # ============== end mapped file ===================
            
            # commit anything left in the copy buffer
            loader.load_variants();

            message = '{:,}'.format(loader.get_count('variant')) + " variants"
            messagePrefix = "COMMITTED"
            if args.commit:
                database.commit()
            else:
                database.rollback()
                messagePrefix = "LOADED"
                message += " -- rolling back"
            message += "; up to = " + loader.get_current_variant_id()
            loader.log(message, prefix=messagePrefix)
            
            if loader.get_count('update') > 0:
                message = '{:,}'.format(loader.get_count('update')) + " variants"
                loader.log(message, prefix="UPDATED")
                            
            if loader.get_count('duplicates') > 0:
                message = '{:,}'.format(loader.get_count('duplicates')) + " variants"    
                loader.log(message, prefix="SKIPPED")
                
            loader.log("DONE", prefix="INFO")
                        
            if args.test:
                loader.log("DONE - TEST COMPLETE" , prefix="WARNING")
                
        # ============== end with open, cursor ===================
        
    except DatabaseError as err:
        loader.log("Problem submitting COPY statement, error involves any variant from last COMMIT until "  \
            + loader.get_current_variant_id(), prefix="ERROR")
        raise(err)
    except IOError as err:
        raise(err)
    except Exception as err:
        loader.log("Problem parsing variant: " + loader.get_current_variant_id(), prefix="ERROR")
        loader.log((lineCount, ":", print_dict(line)), prefix="LINE")
        loader.log(str(err), prefix="ERROR")
        raise(err)
    finally:
        mappedFile.close()
        database.close()
        loader.close()
        print(loader.get_algorithm_invocation_id(), file=stdout)


def validate_args():
    """! validate the parameters, print warnings, and update some values as necessary """
    if args.failAt:
        warning("--failAt option provided / running in NON-COMMIT mode")
        args.commit = False

    if not args.logAfter:
        args.logAfter = args.commitAfter

    if args.fileName and args.chr:
        die("Both --fileName and --chr supplied, please specify only one, see --usage")

    if not args.fileName:
        if not args.chr:
            warning("Neither --fileName of --chr arguments supplied, assuming parallel load of all chromosomes")
            args.chr = 'all'
        else:        
            if not args.dir:
                die("Loading by chromosome, please supply path to directory containing the VEP results")
            if not args.extension:
                die("Loading by chromosome, please provide file extension. Assume file names are like chr1.extension / TODO: pattern matching")
            

if __name__ == "__main__":
    parser = argparse.ArgumentParser(allow_abbrev=False, # otherwise it can substitute --chr for --chromosomeMap
                                    description='load AnnotatedDB from a VCF file, specify either a file or one or more chromosomes')
    parser.add_argument('-d', '--dir',
                        help="directory containing VCF files / only necessary for parallel load")
    parser.add_argument('-e', '--extension', 
                        help="file extension (e.g., vcf.gz) / required for parallel load")
    parser.add_argument('-g', '--genomeBuild', default='GRCh38', help="genome build: GRCh37 or GRCh38")
    parser.add_argument('-s', '--seqrepoProxyPath', required=True,
                        help="full path to local SeqRepo file repository")
    parser.add_argument('-m', '--chromosomeMap', required=False,
                        help="chromosome map")
    parser.add_argument('--commit', action='store_true', help="run in commit mode", required=False)
    parser.add_argument('--gusConfigFile',
                        '--full path to gus config file, else assumes $GUS_HOME/config/gus.config')
    parser.add_argument('--test', action='store_true',
                        help="load 'commitAfter' rows as test")
    parser.add_argument('--resumeAfter',
                        help="variantId after which to resume load (log lists lasts committed variantId)")
    parser.add_argument('-c', '--chr', 
                        help="comma separated list of one or more chromosomes to load, e.g., 1, 2, M, X, `all`, `allNoM` / required for parallel load"),
    parser.add_argument('--fileName',
                        help="full path of file to load, if --dir option is provided, will be ignored")
    parser.add_argument('--commitAfter', type=int, default=500,
                        help="commit after specified inserts")
    parser.add_argument('--maxWorkers', default=10, type=int)
    parser.add_argument('--logAfter', type=int,
                        help="number of inserts to log after completion; will work best if factor/multiple of commitAfter")
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--debug',  action='store_true',
                        help="log database copy time / may print development debug statements")
    parser.add_argument('--failAt', 
                        help="fail on specific variant and log output; if COMMIT = True, COMMIT will be set to False")
    parser.add_argument('--skipExisting', action='store_true',
                        help="check each variant against the database, load non-duplicates only -- time consuming")
    parser.add_argument('--datasource', choices=['dbSNP', 'DBSNP', 'dbsnp', 'ADSP', 'NIAGADS', 'EVA'],
                        default='dbSNP',
                        help="variant source: dbSNP, NIAGADS, ADSP, or EVA (European Variant Archive")
    args = parser.parse_args()

    validate_args()

    chrmMap = ChromosomeMap(args.chromosomeMap) if args.chromosomeMap else None

    if args.fileName:
        load_annotation(args.fileName, args.fileName + "-vcf-variant-loader")
        
    else:
        chrList = args.chr.split(',') if not args.chr.startswith('all') \
            else [c.value for c in Human]

        numChrs = len(chrList)
        with ProcessPoolExecutor(args.maxWorkers) as executor:
            if numChrs > 1:
                for c in chrList:
                    if args.chr == 'allNoM' and c == 'M':
                        continue 
                    warning("Create and start thread for chromosome:", xstr(c))
                    inputFile = path.join(args.dir, 
                        'chr' + xstr(c) + "." + args.extension)
                    executor.submit(load_annotation, fileName=inputFile, logFilePrefix='chr' + xstr(c))
            else: # debugs better w/out thread overhead, so single file -- no threading
                inputFile = path.join(args.dir, 
                        'chr' + xstr(chrList[0]) + "." + args.extension)
                load_annotation(inputFile, 'chr' + xstr(chrList[0]))
