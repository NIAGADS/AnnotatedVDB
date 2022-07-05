#!/usr/bin/env python
# pylint: disable=invalid-name,not-an-iterable,unused-import,too-many-locals
"""
Loads AnnotatedVDB from JSON output from running VEP on dbSNP
Loads each chromosome in parallel
"""

from __future__ import print_function

import argparse
import gzip
import mmap
import json
import glob

from datetime import datetime
from os import path
from concurrent.futures import ProcessPoolExecutor
from sys import stdout
from psycopg2 import DatabaseError

from GenomicsDBData.Util.utils import xstr, warning, print_dict, print_args, die
from GenomicsDBData.Util.postgres_dbi import Database, raise_pg_exception

from AnnotatedVDB.Util.loaders import VEPVariantLoader 
from AnnotatedVDB.Util.parsers import ChromosomeMap
from AnnotatedVDB.Util.enums import HumanChromosome as Human

# ADSP_VARIANT_UPDATE_SQL = "UPDATE AnnotatedVDB.Variant SET is_adsp_variant = true WHERE record_primary_key = %s and chromosome = %s"
ADSP_VARIANT_UPDATE_SQL = """UPDATE AnnotatedVDB.Variant v SET is_adsp_variant = true 
  FROM (VALUES %s) AS d(record_primary_key, chromosome)
  WHERE v.record_primary_key = d.record_primary_key and v.chromosome = d.chromosome"""

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
        loader = VEPVariantLoader(args.datasource, logFileName=lfn, verbose=args.verbose, debug=args.debug)
                
        loader.log(("Parameters:", print_dict(vars(args), pretty=True)), prefix="INFO")

        loader.set_algorithm_invocation('load_vep_result', xstr(logFilePrefix) + '|' + print_args(args, False))
        loader.log('Algorithm Invocation Id = ' + xstr(loader.alg_invocation_id()), prefix="INFO")
        
        loader.initialize_pk_generator(args.genomeBuild, args.seqrepoProxyPath)
        loader.initialize_vep_parser(args.rankingFile, True, args.verbose)
        loader.initialize_bin_indexer(args.gusConfigFile)
        loader.initialize_copy_sql() # use default copy fields
        
        loader.set_chromosome_map(chrmMap)
        
        if loader.is_adsp():
            loader.initialize_variant_validator(args.gusConfigFile)
            loader.set_update_sql(ADSP_VARIANT_UPDATE_SQL)
        
        if args.skipExisting:
            loader.set_skip_existing(True, args.gusConfigFile) # initialize validator db connection
            if args.logSkips:
                loader.log_skips()
            
        if args.failAt:
            loader.set_fail_at_variant(args.failAt)
            loader.log("Fail at variant: " + loader.fail_at_variant(), prefix="INFO")
        
    except Exception as err:
        warning("ERROR", "Problem initializing Loader")
        raise(err)
        
    return loader


def load_annotation(fileName, logFilePrefix):
    """! parse over a JSON file, extract position, frequencies,
    ids, and ADSP-ranked most severe consequence; bulk load using COPY """
 
    loader = initialize_loader(logFilePrefix)
    loader.log('Parsing ' + fileName, prefix="INFO")
    
    resume = args.resumeAfter is None # false if need to skip lines
    if not resume:
        loader.log(("--resumeAfter flag specified; Finding skip until point", args.resumeAfter), prefix="INFO")
        loader.set_resume_after_variant(args.resumeAfter)

    testComplete = False # flag if test is complete
    try: 
        database = Database(args.gusConfigFile)
        database.connect()
        with open(fileName, 'r') as fhandle, database.cursor() as cursor:
            loader.set_cursor(cursor)
            mappedFile = mmap.mmap(fhandle.fileno(), 0, prot=mmap.PROT_READ) # put file in swap
            with gzip.GzipFile(mode='r', fileobj=mappedFile) as gfh:
                lineCount = 0
                for line in gfh:
                    if (lineCount == 0 and not loader.resume_load()) \
                        or (lineCount % args.commitAfter == 0 and loader.resume_load()): 
                        if args.debug:
                            loader.log('Processing new copy object', prefix="DEBUG")
                        tstart = datetime.now()
                            
                    loader.parse_variant(line.rstrip())
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

                        message = 'INSERTED = ' + '{:,}'.format(loader.get_count('variant') - loader.get_count('update') - loader.get_count('duplicates')) 
                        messagePrefix = "COMMITTED"
                        if args.commit:
                            database.commit()
                        else:
                            database.rollback()
                            messagePrefix = "ROLLING BACK"

                        if lineCount % args.logAfter == 0:   
                            if loader.get_count('update') > 0:
                                message += '; UPDATED = ' + '{:,}'.format(loader.get_count('update')) 
                            
                            if loader.get_count('duplicates') > 0:
                                message += '; SKIPPED = ' + '{:,}'.format(loader.get_count('duplicates'))
                                
                            message += "; up to = " + loader.get_current_variant_id()
                            loader.log(message, prefix=messagePrefix)
                            
                            if args.debug:
                                tend = datetime.now()
                                loader.log('Database copy time: ' + str(tend - tendw), prefix="DEBUG")
                                loader.log('        Total time: ' + str(tend - tstart), prefix="DEBUG")
                   
                        if args.test:
                            break

            # ============== end with gzip.GzipFile ===================
            
            # commit anything left in the copy buffer
            loader.load_variants();

            message = 'INSERTED = ' + '{:,}'.format(loader.get_count('variant') - loader.get_count('update') - loader.get_count('duplicates')) 
            messagePrefix = "COMMITTED"
            if args.commit:
                database.commit()
            else:
                database.rollback()
                messagePrefix = "ROLLING BACK"  
            
            if loader.get_count('update') > 0:
                message += '; UPDATED = ' + '{:,}'.format(loader.get_count('update')) 
                            
            if loader.get_count('duplicates') > 0:
                message += '; SKIPPED = ' + '{:,}'.format(loader.get_count('duplicates'))
                
            message += "; up to = " + loader.get_current_variant_id()  
            loader.log(message, prefix=messagePrefix)
                  
            loader.log("DONE", prefix="INFO")
            
            # summarize new consequences
            loader.log(("Consequences added during load:", 
                        loader.vep_parser().get_added_conseq_summary()), prefix="INFO")
            
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
        loader.log((lineCount, ":", print_dict(json.loads(line))), prefix="LINE")
        loader.log(str(err), prefix="ERROR")
        raise(err)
    finally:
        mappedFile.close()
        database.close()

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


def get_input_file_name(chrm):
    """ find the file that matches the chromosome & specified extension """  
    pattern = path.join(args.dir, '*chr' + xstr(chrm) + args.extension)
    fileName = glob.glob(pattern)   # *chr b/c there may be a prefix
    return path.join(args.dir, fileName[0])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(allow_abbrev=False, # otherwise it can substitute --chr for --chromosomeMap
                                    description='load AnnotatedDB from JSON output of VEP against dbSNP, specify either a file or one or more chromosomes')
    parser.add_argument('-d', '--dir',
                        help="directory containing VEP results / only necessary for parallel load")
    parser.add_argument('-e', '--extension', 
                        help="file extension (e.g., json.gz) / required for parallel load")
    parser.add_argument('-g', '--genomeBuild', default='GRCh38', help="genome build: GRCh37 or GRCh38")
    parser.add_argument('-s', '--seqrepoProxyPath', required=True,
                        help="full path to local SeqRepo file repository")
    parser.add_argument('-m', '--chromosomeMap', required=False,
                        help="chromosome map")
    parser.add_argument('-r', '--rankingFile', required=True,
                        help="full path to ADSP VEP consequence ranking file")
    parser.add_argument('--commit', action='store_true', help="run in commit mode", required=False)
    parser.add_argument('--gusConfigFile',
                        '--full path to gus config file, else assumes $GUS_HOME/config/gus.config')
    parser.add_argument('--test', action='store_true',
                        help="load 'commitAfter' rows as test")
    parser.add_argument('--resumeAfter',
                        help="variant after which to resume load (log lists lasts committed variant)")
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
    parser.add_argument('--logSkips', action='store_true',
                        help="log skipped variants")
    parser.add_argument('--datasource', choices=['dbSNP', 'DBSNP', 'dbsnp', 'ADSP', 'NIAGADS', 'EVA'],
                        default='dbSNP',
                        help="variant source: dbSNP, NIAGADS, ADSP, or EVA (European Variant Archive")
    args = parser.parse_args()

    validate_args()

    chrmMap = ChromosomeMap(args.chromosomeMap) if args.chromosomeMap else None

    if args.fileName:
        load_annotation(args.fileName, args.fileName + "-vep-variant-loader")
        
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
                    inputFile = get_input_file_name(c)
                    executor.submit(load_annotation, fileName=inputFile, logFilePrefix='chr' + xstr(c))
            else: # debugs better w/out thread overhead, so single file -- no threading
                inputFile = get_input_file_name(chrList[0])
                load_annotation(inputFile, 'chr' + xstr(chrList[0]))
