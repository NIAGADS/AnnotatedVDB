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

from concurrent.futures import ProcessPoolExecutor, as_completed
from psycopg2 import DatabaseError

from GenomicsDBData.Util.utils import xstr, warning, die, print_dict, print_args
from GenomicsDBData.Util.postgres_dbi import Database, raise_pg_exception

from AnnotatedVDB.Util.parsers import VcfEntryParser
from AnnotatedVDB.Util.loaders import CADDUpdater
from AnnotatedVDB.Util.enums import HumanChromosome as Human

SELECT_SQL = "SELECT record_primary_key, metaseq_id, cadd_scores FROM AnnotatedVDB.Variant WHERE chromosome = %s"

def initialize_loader(logFilePrefix):
    """! initialize loader """

    lfn = path.join(args.logFilePath, xstr(logFilePrefix + ".log"))
    warning("Logging to", lfn)
    
    try:
        loader = CADDUpdater(args.databaseDir, logFileName=lfn, verbose=args.verbose, debug=args.debug)

        loader.log(("Parameters:", print_dict(vars(args), pretty=True)), prefix="INFO")

        loader.set_algorithm_invocation('load_cadd_scores', xstr(logFilePrefix) + '|' + print_args(args, False))
        loader.log('Algorithm Invocation Id = ' + xstr(loader.alg_invocation_id()), prefix="INFO")
        
        loader.initialize_pk_generator(args.genomeBuild, args.seqrepoProxyPath)
        
        if args.skipExisting:
            loader.set_skip_existing(True, args.gusConfigFile) # initializes validator db connection  in addition to setting skip flag
        else:
            loader.initialize_variant_validator(args.gusConfigFile) #otherwise just initialize validator
        
    except Exception as err:
        warning("ERROR", "Problem initializing Loader")
        raise(err)
        
    return loader
    

def update_cadd_scores_by_query(logFilePrefix, chromosome=None, querySql=None):
    loader = initialize_loader(logFilePrefix)
    if querySql is None:
        querySql = SELECT_SQL
    loader.log('Updating by query: ' + querySql, prefix="INFO")
  
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
                loader.log("Executing query: " + querySql, prefix="INFO")
                selectCursor.execute(querySql)
            else:            
                loader.log("Executing query: " + querySql.replace('%s', "'" + chrm + "'"), prefix="INFO")
                selectCursor.execute(querySql, [chrm])
                
            for record in selectCursor:     
                recordCount += 1
                if args.debug and args.veryVerbose:
                    loader.log(("Record:", print_dict(record)), prefix="DEBUG")
                    
                if record['cadd_scores'] is not None:
                    loader.increment_counter('skipped')
                    if args.debug and args.veryVerbose:
                        loader.log(("Skipped", record['record_primary_key']), prefix="DEBUG")
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
                    
                    loader.log(message, prefix=messagePrefix)

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
            
            loader.log(message, prefix=messagePrefix)
            loader.log("DONE", prefix="INFO")
            
            if args.test:
                loader.log("DONE - TEST COMPLETE" , prefix="WARNING")
            
    except DatabaseError as err:
        loader.log("Problem updating variant: "  \
            + loader.get_current_variant(toStr=True), prefix="ERROR")
        raise(err)
    except Exception as err:
        loader.log("Problem parsing variant: " + loader.get_current_variant(toStr=True), prefix="ERROR")
        # loader.log((recordCount, ":", print_dict(record)), prefix="LINE")
        loader.log(str(err), prefix="ERROR")
        raise(err)


    finally:
        database.close()
        loader.close()


def update_cadd_scores_by_vcf(): 
    loader = initialize_loader(args.vcfFile + "-cadd-loader")
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
                    
                    loader.log(message, prefix=messagePrefix)
                    loader.log("Skipped " + '{:,}'.format(loader.get_count('skipped')), prefix="INFO")
                
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
            
            loader.log(message, prefix=messagePrefix)
            loader.log("Skipped " + '{:,}'.format(loader.get_count('skipped')), prefix="INFO")
            loader.log("DONE", prefix="INFO")
            
            if args.test:
                loader.log("DONE - TEST COMPLETE" , prefix="WARNING")             
            
    except DatabaseError as err:
        loader.log("Problem updating variant: "  \
            + loader.get_current_variant(toStr=True), prefix="ERROR")
        raise(err)
    except Exception as err:
        loader.log("Problem parsing variant: " + loader.get_current_variant(toStr=True), prefix="ERROR")
        loader.log((lineCount, ":", line), prefix="LINE")
        loader.log(str(err), prefix="ERROR")
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
    args = parser.parse_args()
    
    if args.vcfFile:
        update_cadd_scores_by_vcf()

    else:
        chrList = args.chr.split(',') if not args.chr.startswith('all') \
            else [c.value for c in Human if c.value not in ['X', 'Y', 'M', 'MT']] if args.chr == 'autosome' \
                else [c.value for c in Human]

        if len(chrList) == 1:
            update_cadd_scores_by_query('db_chr' + xstr(chrList[0]), chromosome=xstr(chrList[0]))

        else:
            random.shuffle(chrList) # so that not all large chrms are done at once if all is selected
            with ProcessPoolExecutor(args.maxWorkers) as executor:
                futureUpdate = {executor.submit(update_cadd_scores_by_query, logFilePrefix='db_chr' + xstr(c),  chromosome=xstr(c)) : c for c in chrList}
                for future in as_completed(futureUpdate): # this should allow catching errors 
                    try:
                        future.result()
                    except Exception as err:
                        raise(err)        


    print("SUCCESS", file=sys.stdout)
