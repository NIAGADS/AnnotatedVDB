#!/usr/bin/env python3
# pylint: disable=invalid-name,not-an-iterable,unused-import,too-many-locals
"""
Loads AnnotatedVDB variants from from VCF File
Loads each chromosome in parallel
"""

import argparse
import mmap
import logging

from datetime import datetime
from os import path
from sys import stdout

from concurrent.futures import ProcessPoolExecutor

from niagads.utils.logging import ExitOnCriticalExceptionHandler, ExitOnExceptionStreamHandler
from niagads.utils.string import xstr
from niagads.utils.dict import print_dict
from niagads.utils.sys import warning, die, get_opener, print_args, verify_path
from niagads.db.postgres import Database, DatabaseError
from niagads.reference.chromosomes import Human

from AnnotatedVDB.Util.loaders import VCFVariantLoader 

LOGGER = logging.getLogger(__name__)

def initialize_logger():
    for handler in logging.root.handlers[:]: # vrs-logging is getting the way
        logging.root.removeHandler(handler)
    
    logFileName = args.fileName + '-load-vcf.log' if args.fileName else path.join(args.dir, 'load-vcf.log')
    
    logHandler = logging.StreamHandler() if args.log2stderr \
        else ExitOnCriticalExceptionHandler(
                filename=logFileName,
                mode='w',
                encoding='utf-8',
            )
    logging.basicConfig(
        handlers=[logHandler],
        format='%(asctime)s %(funcName)s %(levelname)-8s %(message)s',
        level=logging.DEBUG if args.debug else logging.INFO
    )
    LOGGER.info("LOGGING Initialized")
    LOGGER.debug("DEBUGGING enbaled") 
    

def initialize_loader(fileName):
    """! initialize loader """

    try:
        loader = VCFVariantLoader(args.datasource, verbose=args.verbose, debug=args.debug)
        
        if args.verbose:
            LOGGER.info("Parameters: %", print_dict(vars(args), pretty=True))

        loader.set_algorithm_invocation('load_vcf_result', print_args(args, False), commit=args.commit)
        LOGGER.info('Algorithm Invocation Id = ' + xstr(loader.alg_invocation_id()))
        
        loader.initialize_pk_generator(args.genomeBuild, args.seqrepoProxyPath)
        loader.initialize_bin_indexer(args.gusConfigFile)
        loader.initialize_copy_sql() # use default copy fields
        
        if args.logExisting: loader.log_skips()
        loader.set_skip_existing(args.skipExisting, args.gusConfigFile) # initialize validator db connection
        
        if args.failAt:
            loader.set_fail_at_variant(args.failAt)
            LOGGER.info("Fail at variant: " + loader.fail_at_variant())
        
    except Exception as err:
        LOGGER.critical("Problem initializing Loader")
        raise RuntimeError("Problem initializing Loader")
        
    return loader


def load(fileName):
    """! parse over a VCF file; bulk load using COPY """

    loader = initialize_loader(fileName)
    LOGGER.info('Parsing ' + fileName)
    LOGGER.info('Writing metaseq_id -> primary_key mapping to ' + fileName + '.mapping')
    
    resume = args.resumeAfter is None # false if need to skip lines
    if not resume:
        LOGGER.info("--resumeAfter flag specified; Finding skip until point %", args.resumeAfter)
        loader.set_resume_after_variant(args.resumeAfter)
        
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
                    tstart = datetime.now()
                    if args.debug:
                        LOGGER.debug('Processing new copy object')

                primaryKeyMapping = loader.parse_variant(line.rstrip())

                if args.debug:
                    LOGGER.debug("PKM-%s: %s", lineCount, primaryKeyMapping)
                for metaseqId, pk in primaryKeyMapping.items():
                    print(metaseqId, pk, sep='\t', file=mfh, flush=True)

                lineCount += 1
                    
                if not loader.resume_load():
                    if lineCount % args.logAfter == 0:
                        LOGGER.info('SKIPPED = {:,} lines'.format(lineCount))
                    continue
                    
                if loader.resume_load() != resume: # then you are at the resume cutoff
                    resume = True
                    LOGGER.info('SKIPPED = {:,} lines'.format(lineCount))
                    continue
                
                if lineCount % args.logAfter == 0 and lineCount % args.commitAfter != 0:
                    LOGGER.info("Parsed %s lines (%s variants)",
                        lineCount, loader.get_count('variant'))

                if lineCount % args.commitAfter == 0:
                    if args.debug:
                        tendw = datetime.now()
                        message = 'Copy object prepared in ' + str(tendw - tstart) + '; ' + \
                            str(loader.copy_buffer(sizeOnly=True)) + ' bytes; transfering to database'
                        LOGGER.debug(message)
                    
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
                        message += "; PARSED: " + xstr(lineCount) + "; up to = " + loader.get_current_variant_id()

                        if loader.get_count('update') > 0:
                            message += '; UPDATED =  {:,}'.format(loader.get_count('update')) + " variants"
                        
                        if loader.get_count('skipped') > 0:
                            message += '; SKIPPED = {:,}'.format(loader.get_count('skipped')) + " variants"  

                        LOGGER.info("%s: %s", messagePrefix, message)
                        if args.debug:
                            tend = datetime.now()
                            LOGGER.debug('Database copy time: ' + str(tend - tendw))
                            LOGGER.debug('        Total time: ' + str(tend - tstart))

                    if args.test:
                        break

            mappedFile.close()
            
            # ============== end mapped file ===================
            
            # commit anything left in the copy buffer
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
            message += "; up to = " + loader.get_current_variant_id()

            
            if loader.get_count('update') > 0:
                message += '; UPDATED =  {:,}'.format(loader.get_count('update')) + " variants"
            
            if loader.get_count('skipped') > 0:
                message += '; SKIPPED = {:,}'.format(loader.get_count('skipped')) + " variants"    
    
                
            LOGGER.info("%s: %s", messagePrefix, message)
            LOGGER.info("DONE")
                        
            if args.test:
                LOGGER.info("DONE - TEST COMPLETE")
                
        # ============== end with open, cursor ===================
        
    except DatabaseError as err:
        LOGGER.critical("Problem submitting COPY statement, error involves any variant from last COMMIT until "  \
            + loader.get_current_variant_id())
        raise(err)
    except IOError as err:
        raise(err)
    except Exception as err:
        LOGGER.critical("Problem parsing variant: %s; line: %s", loader.get_current_variant_id(), print_dict(line))
        raise(err)
    finally:
        mappedFile.close()
        database.close()
        print(loader.get_algorithm_invocation_id(), file=stdout)
        loader.close()
        

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
    parser.add_argument('--datasource', choices=['dbSNP', 'DBSNP', 'dbsnp', 'ADSP', 'ADSP-FunGen', 'NIAGADS', 'EVA'],
                        default='dbSNP',
                        help="variant source: dbSNP, NIAGADS, ADSP, ADSP-FunGen or EVA (European Variant Archive")
    parser.add_argument('--logExisting', action='store_true')
    parser.add_argument('--log2stderr', action="store_true")
    
    args = parser.parse_args()

    validate_args()
    
    initialize_logger()
    
    if args.fileName:
        if verify_path(args.fileName):
            load(args.fileName)
        else:
            LOGGER.info(f'Input file {args.fileName} not found.  EXITING.')
        
    else:
        chrList = args.chr.split(',') if not args.chr.startswith('all') \
            else [c.value for c in Human]

        if len(chrList) == 1:
            inputFile = path.join(args.dir, 'chr' + xstr(chrList[0]) + "." + args.extension)
            load(inputFile)
                
        else:
            with ProcessPoolExecutor(args.maxWorkers) as executor:
                for c in chrList:
                    if args.chr == 'allNoM' and c == 'M':
                        continue 
                    inputFile = path.join(args.dir, 
                        'chr' + xstr(c) + "." + args.extension)
                    executor.submit(load, fileName=inputFile)

