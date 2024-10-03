#!/usr/bin/env python3 
#pylint: disable=invalid-name

'''
Generates load file from CADD download (SNVs)
Looks up database variants by chromome and uses pysam
to map against the CADD tab-indexed SNV file
'''

import argparse
import os.path as path
import random
import sys
import logging

from concurrent.futures import ProcessPoolExecutor, as_completed

from niagads.utils.string import xstr
from niagads.utils.dict import print_dict
from niagads.utils.sys import warning, print_args
from niagads.db.postgres import Database, DatabaseError
from niagads.utils.logging import ExitOnCriticalExceptionHandler
from niagads.reference.chromosomes import Human

from AnnotatedVDB.Util.parsers import VcfEntryParser
from AnnotatedVDB.Util.loaders import CADDUpdater


SELECT_SQL = "SELECT record_primary_key, metaseq_id, cadd_scores FROM AnnotatedVDB.Variant WHERE chromosome = %s"

LOGGER = logging.getLogger(__name__)

def initialize_logger(fileName):
    
    for handler in logging.root.handlers[:]: # vrs-logging is getting the way
        logging.root.removeHandler(handler)
    logFileName = fileName + '-load-vef.log' 
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

    return logHandler


def initialize_loader(fileName):
    """! initialize loader """

    try:
        initialize_logger(fileName)
        loader = CADDUpdater(args.databaseDir, verbose=args.verbose, debug=args.debug)

        LOGGER.info("Parameters: %s", print_dict(vars(args), pretty=True))

        loader.set_algorithm_invocation('load_cadd_scores', print_args(args, False), commit=args.commit)
        LOGGER.info('Algorithm Invocation Id = ' + xstr(loader.alg_invocation_id()))
        
        loader.initialize_pk_generator(args.genomeBuild, args.seqrepoProxyPath)
        
        if args.skipExisting:
            loader.set_skip_existing(True, args.gusConfigFile) # initializes validator db connection  in addition to setting skip flag
        else:
            loader.initialize_variant_validator(args.gusConfigFile) #otherwise just initialize validator
        
    except Exception as err:
        warning("ERROR", "Problem initializing Loader")
        raise(err)
        
    return loader
    

def update_cadd_scores_by_query(chromosome=None, querySql=None):
    loader = initialize_loader("cadd-updated-by-query")
    
    if querySql is None:
        querySql = SELECT_SQL
    LOGGER.info('Updating by query: ' + querySql)

    try: 
        database = Database(args.gusConfigFile)
        database.connect()
        recordCount = 0
        cname = 'update'
        
        chrm = chromosome
        if chromosome is not None:
            chrm = chromosome if 'chr' in xstr(chromosome) else 'chr' + xstr(chromosome)
            cname += '_' + chrm
            loader.set_chromosome(chrm)
            
        with database.cursor() as updateCursor, database.named_cursor(cname, cursorFactory="RealDictCursor") as selectCursor:
            selectCursor.itersize = 500000 #args.commitAfter
                        
            loader.set_cursor(updateCursor)            
            if chromosome is None:
                LOGGER.info("Executing query: " + querySql)
                selectCursor.execute(querySql)
            else:            
                LOGGER.info("Executing query: " + querySql.replace('%s', "'" + chrm + "'"))
                selectCursor.execute(querySql, [chrm])
                
            for record in selectCursor:     
                recordCount += 1
                if args.debug and args.veryVerbose:
                    LOGGER.debug("Record: %s", print_dict(record))
                    
                if record['cadd_scores'] is not None:
                    loader.increment_counter('skipped')
                    if args.debug and args.veryVerbose:
                        LOGGER.debug("Skipped %s", record['record_primary_key'])
                else: 
                    loader.set_current_variant(record)
                    loader.buffer_variant()
                    if recordCount % args.commitAfter == 0:
                        if loader.update_buffer(sizeOnly=True) > 0:
                            loader.update_variants()
                            if args.commit:
                                database.commit()
                            else:
                                database.rollback()
                        
                if recordCount % args.logAfter == 0:
                    message = ' '.join(("Parsed", '{:,}'.format(recordCount), "variants",
                                "- SNVS:", '{:,}'.format(loader.get_count('snv')),
                                "- INDELs:", '{:,}'.format(loader.get_count('indel')),
                                "- No match", '{:,}'.format(loader.get_count('not_matched')),
                                "- Skipped", '{:,}'.format(loader.get_count('skipped'))))
                    
                    messagePrefix = "COMMITTED" if args.commit else "ROLLING BACK"
                    
                    LOGGER.info("%s: %s", messagePrefix, message)

                    if args.test:
                        break
                    
            # ======================= end iterate over records ==================================
            
            # clear out buffer
            if loader.update_buffer(sizeOnly=True) > 0:
                loader.update_variants() 
            message = ' '.join(("Updated", '{:,}'.format(recordCount), "variants",
                        "- SNVS:", '{:,}'.format(loader.get_count('snv')),
                        "- INDELs:", '{:,}'.format(loader.get_count('indel')),
                        "- No match", '{:,}'.format(loader.get_count('not_matched')),
                        "- Skipped", '{:,}'.format(loader.get_count('skipped'))))
            
            messagePrefix = "COMMITTED"
            if args.commit:
                database.commit()
            else:
                database.rollback()
                messagePrefix = "ROLLING BACK"
            
            LOGGER.info("%s: %s", messagePrefix, message)
            LOGGER.info("DONE")
            
            if args.test:
                LOGGER.info("DONE - TEST COMPLETE")
            
    except DatabaseError as err:
        LOGGER.critical("Problem updating variant: "  \
            + loader.get_current_variant(toStr=True))
        raise(err)
    except Exception as err:
        LOGGER.critical("Problem parsing variant: " + loader.get_current_variant(toStr=True))
        raise(err)

    finally:
        database.close()
        loader.close()


def update_cadd_scores_by_vcf(): 
    loader = initialize_loader(args.vcfFile)

    try: 
        database = Database(args.gusConfigFile)
        database.connect()
        lineCount = 0
        with database.cursor() as updateCursor, open(args.vcfFile, 'r') as fh:
            loader.set_cursor(updateCursor)    

            for line in fh:
                lineCount += 1
                if line.startswith("#"): continue
                entry = VcfEntryParser(line.rstrip()) 
                loader.set_current_variant(entry.get_variant(), 'id')
                loader.buffer_variant()
                
                if lineCount % args.commitAfter == 0:
                    loader.update_variants()         
                    message = ' '.join(("Updated", '{:,}'.format(lineCount), "variants",
                                "- SNVS:", '{:,}'.format(loader.get_count('snv')),
                                "- INDELs:", '{:,}'.format(loader.get_count('indel')),
                                "- No match", '{:,}'.format(loader.get_count('not_matched'))))
                    
                    messagePrefix = "COMMITTED"
                    if args.commit:
                        database.commit()
                    else:
                        database.rollback()
                        messagePrefix = "ROLLING BACK"
                    
                    if loader.get_count('skipped') > 0:
                        message += "; SKIPPED " + '{:,}'.format(loader.get_count('skipped'))
                
                    LOGGER.info("%s: %s", messagePrefix, message)
                
                    if args.test:
                        break

            # ======================= end iterate over records ==================================

            # clear out buffer
            if (loader.update_buffer(sizeOnly=True)):
                loader.update_variants() 
            
            message = ' '.join(("Updated", '{:,}'.format(lineCount), "variants",
                        "- SNVS:", '{:,}'.format(loader.get_count('snv')),
                        "- INDELs:", '{:,}'.format(loader.get_count('indel')),
                        "- No match", '{:,}'.format(loader.get_count('not_matched'))))
            
            messagePrefix = "COMMITTED"
            if args.commit:
                database.commit()
            else:
                database.rollback()
                messagePrefix = "ROLLING BACK"
                
            if loader.get_count('skipped') > 0:
                message += "; SKIPPED " + '{:,}'.format(loader.get_count('skipped'))
                
            LOGGER.info("%s: %s", messagePrefix, message)

            LOGGER.info("DONE")
            
            if args.test:
                LOGGER.info("DONE - TEST COMPLETE")
            
    except DatabaseError as err:
        LOGGER.critical("Problem updating variant: "  \
            + loader.get_current_variant(toStr=True), prefix="ERROR")
        raise(err)
    except Exception as err:
        LOGGER.critical("Problem parsing variant: %s", loader.get_current_variant(toStr=True))
        raise(err)
    finally:
        database.close()
        loader.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="update CADD scores in DB specified chromosome", allow_abbrev=False)
    parser.add_argument('-d', '--databaseDir', required=True,
                        help="directory containg CADD database +tabindex files")
    parser.add_argument('--vcfFile', 
                        help="if file is specified, updates only the listed variants; otherwise updates all variants in the database for the specified chromosome; expects full path to a VCF file")
    parser.add_argument('--chr', default='all',
                        help="chromosome; comma separated list of one or more chromosomes or 'all' or 'autosome'; ignored if --file is specified")
    parser.add_argument('--maxWorkers', type=int, default=5)
    parser.add_argument('--commitAfter', type=int, default=500)
    parser.add_argument('--logAfter', type=int, default=500000)
    parser.add_argument('--logFilePath', help='generate formulaic log files and store in the specified path')
    parser.add_argument('-g', '--genomeBuild', default='GRCh38', help="genome build: GRCh37 or GRCh38")
    parser.add_argument('-s', '--seqrepoProxyPath', required=True,
                        help="full path to local SeqRepo file repository")
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--veryVerbose', action='store_true')
    parser.add_argument('--test', action='store_true')
    parser.add_argument('--skipExisting', action='store_true')
    parser.add_argument('--gusConfigFile',
                        help="GUS config file. If not provided, assumes default: $GUS_HOME/conf/gus.config")
    parser.add_argument('--commit', action='store_true')
    parser.add_argument('--log2stderr', action="store_true")
    args = parser.parse_args()
    
    if args.vcfFile:
        update_cadd_scores_by_vcf()

    else:
        chrList = None
        if args.chr not in ['all', 'autosome']:
            chrList = args.chr.split(',')
        elif args.chr == 'autosome':
            chrList = [c.value for c in Human if c.value not in ['X', 'Y', 'M', 'MT']] 
        else: 
            chrList = [c.value for c in Human]

        if len(chrList) == 1:
            update_cadd_scores_by_query(xstr(chrList[0]))

        else:
            random.shuffle(chrList) # so that not all large chrms are done at once if all is selected
            with ProcessPoolExecutor(args.maxWorkers) as executor:
                futureUpdate = {executor.submit(update_cadd_scores_by_query, chromosome=xstr(c)) : c for c in chrList}
                for future in as_completed(futureUpdate): # this should allow catching errors 
                    try:
                        future.result()
                    except Exception as err:
                        raise(err)        


    print("SUCCESS", file=sys.stdout)
