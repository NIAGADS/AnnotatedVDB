#!/usr/bin/env python 
#pylint: disable=invalid-name

'''
Moves external GRCh37 AnnotatedVDB database into the GRCh37 GUS database,
'''

import argparse
import os.path as path
import random
from datetime import datetime
import json

from sys import stdout
from concurrent.futures import ProcessPoolExecutor, as_completed
from psycopg2 import DatabaseError

from GenomicsDBData.Util.utils import xstr, warning, die, print_dict, print_args
from GenomicsDBData.Util.postgres_dbi import Database, raise_pg_exception

from AnnotatedVDB.Util.loaders import PatchVEPVariantLoader 
from AnnotatedVDB.Util.enums import HumanChromosome as Human

# create connection to external AnnotatedVDB --> need it's gus.config in a param
# iterate over chromosomes

SELECT_SQL = """SELECT ref_snp_id, 
metaseq_id, 
REPLACE(chromosome, 'chr', '') AS chromosome, 
location AS position, 
split_part(metaseq_id, ':', 3) AS ref_allele,
split_part(metaseq_id, ':', 4) AS alt_allele,
CASE WHEN is_multi_allelic IS FALSE THEN NULL ELSE is_multi_allelic END AS is_multi_allelic,
CASE WHEN is_adsp_variant IS FALSE THEN NULL ELSE is_adsp_variant END AS is_adsp_variant,
cadd_scores,
other_annotation - 'GenomicsDB' - 'GRCh38' AS adsp_qc, 
CASE WHEN other_annotation->'GRCh38' IS NOT NULL THEN jsonb_build_object('GRCh38', other_annotation->'GRCh38') ELSE NULL END AS other_annotation,
vep_output 
FROM Variant 
WHERE chromosome = %s"""

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
        loader = PatchVEPVariantLoader("patch", logFileName=lfn, verbose=args.verbose, debug=args.debug)
                
        loader.log(("Parameters:", print_dict(vars(args), pretty=True)), prefix="INFO")
        
        loader.set_algorithm_invocation('load_vep_result', xstr(logFilePrefix) + '|' + print_args(args, False))
        loader.log('Algorithm Invocation Id = ' + xstr(loader.alg_invocation_id()), prefix="INFO")
        
        
        loader.initialize_pk_generator(args.genomeBuild, args.seqrepoProxyPath)
        warning("initialized pk generator")
        loader.initialize_vep_parser(args.rankingFile, True, args.verbose)
        warning("initialized vep parser")
        loader.initialize_bin_indexer(args.gusConfigFile)
        warning("initialized bin indexer")
        loader.initialize_copy_sql(['is_adsp_variant', 'adsp_qc', 'cadd_scores', 'other_annotation'], appendToDefault=True) # use default copy fields
        warning("initialized bin copy sql")
        if args.failAt:
            loader.set_fail_at_variant(args.failAt)
            loader.log("Fail at variant: " + loader.fail_at_variant(), prefix="INFO")

    except DatabaseError as err:
        warning("ERROR", "Problem initializing Loader")
        print(err)
        raise(err)    
    except Exception as err:
        warning("ERROR", "Problem initializing Loader")
        warning(err)
        raise(err)
        
    return loader


def validate_args():
    """! validate the parameters, print warnings, and update some values as necessary """
    if args.failAt:
        warning("--failAt option provided / running in NON-COMMIT mode")
        args.commit = False

    if not args.logAfter:
        args.logAfter = args.commitAfter

    if not args.chr:
        warning("--chr argument not supplied, assuming parallel load of all chromosomes")
        args.chr = 'all'

# original fields
'''chromosome
location
is_multi_allelic
is_adsp_variant
ref_snp_id
metaseq_id
bin_index
allele_frequencies
cadd_scores
adsp_most_severe_consequence
adsp_ranked_consequences
loss_of_function
vep_output
other_annotation
row_algorithm_id '''

# new fields 


def patch_annotation(logFilePrefix, chromosome):
    warning("Running patch")
    loader = initialize_loader(logFilePrefix)
    loader.log('Patching chr' + xstr(chromosome), prefix="INFO")

    # resume = args.resumeAfter is None # false if need to skip lines
    # if not resume:
    #    loader.log(("--resumeAfter flag specified; Finding skip until point", args.resumeAfter), prefix="INFO")
    #    loader.set_resume_after_variant(args.resumeAfter)

    try: 
        annotatedVDB = Database(args.annotatedVDBgusConfigFile)
        annotatedVDB.connect()
        database = Database(args.gusConfigFile)
        database.connect()

        chrm = chromosome
        cname = 'insert'
        if chromosome is not None:
            chrm = chromosome if 'chr' in xstr(chromosome) else 'chr' + xstr(chromosome)
            cname += '_' + chrm
            # loader.set_chromosome(chrm)

        with database.cursor() as insertCursor, annotatedVDB.named_cursor(cname, cursorFactory="RealDictCursor") as selectCursor:
            selectCursor.itersize = 500000 #args.commitAfter        
            loader.set_cursor(insertCursor)         
            recordCount = 0
            
            if args.debug:
                loader.log("Executing query: " + SELECT_SQL.replace('%s', "'" + chrm + "'"), prefix="DEBUG")
            else:
                loader.log("Retrieving variants from old AnnotatedVDB for " + chrm, prefix="INFO")
            selectCursor.execute(SELECT_SQL, [chrm])
            for record in selectCursor:     
                if (recordCount == 0 or recordCount % args.commitAfter == 0): 
                    if args.debug:
                        loader.log('Processing new copy object', prefix="DEBUG")
                    tstart = datetime.now()

                if args.debug and args.verbose:
                    loader.log(("Record:", print_dict(record)), prefix="DEBUG")
                recordCount += 1
                    
                loader.parse_variant(record)

                if recordCount % args.logAfter == 0 \
                    and recordCount % args.commitAfter != 0:
                    loader.log((recordCount, "records (", 
                                loader.get_count('variant'), " variants)"), 
                                prefix="PARSED")

                if recordCount % args.commitAfter == 0:
                    if args.debug:
                        tendw = datetime.now()
                        message = 'Copy object prepared in ' + str(tendw - tstart) + '; ' + \
                            str(loader.copy_buffer(sizeOnly=True)) + ' bytes; transfering to database'
                        loader.log(message, prefix="DEBUG") 
                    
                    loader.load_variants()

                    message = 'INSERTED = ' + '{:,}'.format(loader.get_count('variant') - loader.get_count('update') - loader.get_count('duplicates')) 
                    messagePrefix = "COMMITTED"
                    if args.commit:
                        database.commit()
                    else:
                        database.rollback()
                        messagePrefix = "ROLLING BACK"

                    if recordCount % args.logAfter == 0:   
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

            # ============== end SELECT ===================
            
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
        loader.log((recordCount, ":", print_dict(record)), prefix="LINE")
        loader.log(str(err), prefix="ERROR")
        raise(err)
    finally:
        database.close()
        annotatedVDB.close()

    print(loader.get_algorithm_invocation_id(), file=stdout)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(allow_abbrev=False, # otherwise it can substitute --chr for --chromosomeMap
                                    description='load AnnotatedDB from JSON output of VEP against dbSNP, specify either a file or one or more chromosomes')
    parser.add_argument('-g', '--genomeBuild', default='GRCh37', help="genome build: GRCh37 or GRCh38")
    parser.add_argument('-s', '--seqrepoProxyPath', required=True,
                        help="full path to local SeqRepo file repository")
    parser.add_argument('-r', '--rankingFile', required=True,
                        help="full path to ADSP VEP consequence ranking file")
    parser.add_argument('--commit', action='store_true', help="run in commit mode", required=False)
    parser.add_argument('--gusConfigFile',
                        help="full path to gus config file, else assumes $GUS_HOME/config/gus.config")
    parser.add_argument('--annotatedVDBgusConfigFile', required=True,
                        help="full path to external AnnotatedVDB gus config file, else assumes")
    parser.add_argument('--test', action='store_true',
                        help="load 'commitAfter' rows as test")
    parser.add_argument('--resumeAfter',
                        help="variant after which to resume load (log lists lasts committed variant)")
    parser.add_argument('-c', '--chr', 
                        help="comma separated list of one or more chromosomes to load, e.g., 1, 2, M, X, `all`, `allNoM` / required for parallel load"),
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
    args = parser.parse_args()

    validate_args()

    chrList = args.chr.split(',') if not args.chr.startswith('all') \
        else [c.value for c in Human]

    if len(chrList) == 1:
        patch_annotation('chr' + xstr(chrList[0]), chrList[0])
    else:
        with ProcessPoolExecutor(args.maxWorkers) as executor:      
            for c in chrList:
                if args.chr == 'allNoM' and c == 'M':
                    continue 
                warning("Create and start thread for chromosome:", xstr(c))
                executor.submit(patch_annotation, logFilePrefix='chr' + xstr(c), chromosome=c)
        