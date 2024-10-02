#!/usr/bin/env python3
# pylint: disable=invalid-name,not-an-iterable,unused-import,too-many-locals

import argparse

from niagads.reference.chromosomes import Human
from niagads.db.postgres import Database, QueryCanceledError
from niagads.utils.sys import warning
from niagads.utils.string import xstr

if __name__ == "__main__":
    parser = argparse.ArgumentParser(allow_abbrev=False, # otherwise it can substitute --chr for --chromosomeMap
                                    description='undo a variant load')
    parser.add_argument('--commit', action='store_true', help="run in commit mode", required=False)
    parser.add_argument('--gusConfigFile',
                        '--full path to gus config file, else assumes $GUS_HOME/config/gus.config')
    parser.add_argument('--chr', default='all')
    parser.add_argument('--algInvocId', required=True, type=int)
    args = parser.parse_args()
        
    sql = """
    WITH variants AS (
        SELECT record_primary_key FROM AnnotatedVDB.Variant 
        WHERE row_algorithm_id = %s and chromosome = %s
        LIMIT @LIMIT@
    )
    DELETE FROM AnnotatedVDB.Variant 
    WHERE record_primary_key IN (SELECT record_primary_key FROM variants)
    AND chromosome = %s
    """
    
    chromosomes = [args.chr] if args.chr != 'all' else [c.name for c in Human]
    totalDeletions = 0

    try:
        database = Database(args.gusConfigFile)
        database.connect(timeout=5000)
        with database.cursor() as cursor:
            for chr in chromosomes:
                warning("Processing chromosome = " + chr)
                limit = 500
                done = False
                iterationCount = 0
                chrDeletions = 0
                while not done:
                    iterationCount += 1
                    try:
                        cursor.execute(sql.replace("@LIMIT@", xstr(limit)), (args.algInvocId, chr, chr))
                        if cursor.rowcount > 0:
                            chrDeletions += cursor.rowcount
                            message = "Iteration: " + xstr(iterationCount) + " - DELETED: " + xstr(chrDeletions) + " | row_algorithm_id = " + xstr(args.algInvocId) + "; chromosome = " + chr
                            if args.commit:
                                database.commit()
                                warning("COMMIT: " + message)
                            else:
                                warning(message)
                        else:
                            warning("Iteration: " + xstr(iterationCount) + " - DONE")
                            done = True
                    except QueryCanceledError as err:
                        warning("DELETION query timed-out; adjusting limit")
                        database.rollback() # rollback the error transaction
                        if limit == 1: # no results probably
                            warning("No more records found for " + chr)
                            done = True
                        limit = 50 if limit == 500 else 1
                        warning("LIMIT " + xstr(limit))
                    
                
                totalDeletions += chrDeletions
                
                if not args.commit:
                    database.rollback()
                    warning("ROLLING BACK")

        warning("DONE: Deleted " + xstr(totalDeletions) + " rows")
        
    except Exception as err:
        raise err
    finally:
        database.close()