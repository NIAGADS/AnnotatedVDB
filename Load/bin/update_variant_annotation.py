#!/usr/bin/env python
# pylint: disable=invalid-name,not-an-iterable,unused-import,too-many-locals
"""
Loads AnnotatedVDB from variants in a txt file, which includes the following fields
variant -- record_primary_key, metaseq_id, or ref_snp_id
one or more fields matching columns in the table AnnotatedVDB.Variant
If the variant exists in the database, will update the fields (overwrite, text, numeric fields & concatenate JSONB fields)
If the variant does not exist will either throw error or load depending on options supplied at run time
"""

from __future__ import print_function

import argparse
import csv

from copy import deepcopy
from os import path
from sys import stdout
from psycopg2 import DatabaseError

from GenomicsDBData.Util.utils import xstr, warning, print_dict, print_args, die, get_opener
from GenomicsDBData.Util.postgres_dbi import Database, raise_pg_exception

from AnnotatedVDB.Util.database import VARIANT_ID_TYPES
from AnnotatedVDB.Util.loaders import TextVariantLoader 


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
        loader = TextVariantLoader(args.datasource, logFileName=lfn, verbose=args.verbose, debug=args.debug)
        
        if args.verbose:
            loader.log(("Parameters:", print_dict(vars(args), pretty=True)), prefix="INFO")

        loader.set_algorithm_invocation('update_variant_annotation', xstr(logFilePrefix) + '|' + print_args(args, False))
        loader.log('Algorithm Invocation Id = ' + xstr(loader.alg_invocation_id()), prefix="INFO")
        
        loader.initialize_bin_indexer(args.gusConfigFile)

        loader.set_variant_id_type(args.variantIdType)
        loader.initialize_variant_validator(args.gusConfigFile, args.useDynamicPkSql)

        if args.failAt:
            loader.set_fail_at_variant(args.failAt)
            loader.log("Fail at variant: " + loader.fail_at_variant(), prefix="INFO")
        
    except Exception as err:
        warning("ERROR", "Problem initializing Loader")
        raise(err)
        
    return loader


def update_annotation(fileName):
    """! parse over a text file; bulk update using execute_values """
 
    loader = initialize_loader(fileName + "-update-annotation")
    loader.log('Parsing ' + fileName, prefix="INFO")
    loader.set_update_existing(True)
    
    resume = args.resumeAfter is None # false if need to skip lines
    if not resume:
        loader.log(("--resumeAfter flag specified; Finding skip until point", args.resumeAfter), prefix="INFO")
        loader.set_resume_after_variant(args.resumeAfter)

    try: 
        database = Database(args.gusConfigFile)
        database.connect()
        opener = get_opener(args.fileName)
        with opener(fileName, 'r') as fhandle, database.cursor() as cursor:
            loader.set_cursor(cursor)
            lineCount = 0
            
            reader = csv.DictReader(fhandle, delimiter='\t')      
            updateFields = deepcopy(reader.fieldnames)
            updateFields.remove('variant')
            loader.set_update_fields(updateFields)      
            
            for row in reader:
                loader.parse_variant(row)
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
                    
                    loader.build_update_sql(args.useDynamicPkSql)
                    loader.update_variants()
                    
                    message = '{:,}'.format(loader.get_count('variant')) + " variants"
                    messagePrefix = "COMMITTED"
                    if args.commit:
                        database.commit()
                    else:
                        database.rollback()
                        messagePrefix = "ROLLING BACK"

                    if lineCount % args.logAfter == 0:
                        message += "; up to = " + loader.get_current_variant_id()
                        
                        if loader.get_count('update') > 0:
                            message += '; UPDATED {:,}'.format(loader.get_count('update')) + " variants"
                                                    
                        if loader.get_count('skipped') > 0:
                            message += '; SKIPPED {:,}'.format(loader.get_count('skipped')) + " variants"    

                        loader.log(message, prefix=messagePrefix)
                    if args.test:
                        break


            # ============== end mapped file ===================
            
            # commit anything left in the update buffer
            loader.build_update_sql(args.useDynamicPkSql)
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
            # loader.log(message, prefix=messagePrefix)
            
            if loader.get_count('update') > 0:
                message += '; UPDATED {:,}'.format(loader.get_count('update')) + " variants"
                            
            if loader.get_count('skipped') > 0:
                message += '; SKIPPED {:,}'.format(loader.get_count('skipped')) + " variants"    
                
            loader.log(message, prefix=messagePrefix)
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
        loader.log((lineCount, ":", print_dict(row)), prefix="LINE")
        loader.log(str(err), prefix="ERROR")
        raise(err)
    finally:
        database.close()
        loader.close()
        print('SUCCESS', file=stdout)


def validate_args():
    """! validate the parameters, print warnings, and update some values as necessary """
    if args.failAt:
        warning("--failAt option provided / running in NON-COMMIT mode")
        args.commit = False

    if not args.logAfter:
        args.logAfter = args.commitAfter
            

if __name__ == "__main__":
    parser = argparse.ArgumentParser(allow_abbrev=False, # otherwise it can substitute --chr for --chromosomeMap
                                    description='load AnnotatedDB from a VCF file, specify either a file or one or more chromosomes')

    parser.add_argument('--commit', action='store_true', help="run in commit mode", required=False)
    parser.add_argument('--gusConfigFile',
                        '--full path to gus config file, else assumes $GUS_HOME/config/gus.config')
    parser.add_argument('--test', action='store_true',
                        help="load 'commitAfter' rows as test")
    parser.add_argument('--variantIdType', choices=VARIANT_ID_TYPES)
    parser.add_argument('--resumeAfter',
                        help="variantId after which to resume load (log lists lasts committed variantId)")
    parser.add_argument('--fileName', required=True,
                        help="full path of file to load, if --dir option is provided, will be ignored")
    parser.add_argument('--commitAfter', type=int, default=500,
                        help="commit after specified inserts")
    parser.add_argument('--useDynamicPkSql', action='store_true',
                        help="for legacy version, use sql that dynamically infers PK")
    parser.add_argument('--maxWorkers', default=10, type=int)
    parser.add_argument('--logAfter', type=int,
                        help="number of inserts to log after completion; will work best if factor/multiple of commitAfter")
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--debug',  action='store_true',
                        help="log database copy time / may print development debug statements")
    parser.add_argument('--failAt', 
                        help="fail on specific variant and log output; if COMMIT = True, COMMIT will be set to False")
    parser.add_argument('--datasource', choices=['dbSNP', 'DBSNP', 'dbsnp', 'ADSP', 'NIAGADS', 'EVA'],
                        default='dbSNP',
                        help="variant source: dbSNP, NIAGADS, ADSP, or EVA (European Variant Archive")
    args = parser.parse_args()

    validate_args()

    update_annotation(args.fileName)
   