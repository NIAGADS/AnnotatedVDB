#!/usr/bin/env python 
#pylint: disable=invalid-name

'''
Generates load file from CADD download (SNVs)
Looks up database variants by chromome and uses pysam
to map against the CADD tab-indexed SNV file
'''

from __future__ import with_statement
from __future__ import print_function

import argparse
import os.path as path
import gzip
import csv
import json
import pysam
import random
import psycopg2
import sys

from CBILDataCommon.Util.utils import qw, xstr, warning, die
from CBILDataCommon.Util.postgres_dbi import Database

from AnnotatedVDB.Util.cadd_updater import CADDUpdater
from AnnotatedVDB.Util.vcf_parser import VcfEntryParser

import multiprocessing
from concurrent.futures import ProcessPoolExecutor

from AnnotatedVDB.Util.chromosomes import Human


def update_variant_record(record, updater, isIndel): 
    metaseqId = record['metaseq_id']

    # if args.debug:
    #    warning("Matching", metaseqId, file=updater.lfh(), flush=True)

    chrm, position, refAllele, altAllele = metaseqId.split(':')

    matched = False
    for match in updater.match(chrm, position, isIndel):
        matchedAlleles = [match[2], match[3]]

        if refAllele in matchedAlleles and altAllele in matchedAlleles:
            evidence = {'CADD_raw_score': float(match[4]), 'CADD_phred': float(match[5])}
            updater.buffer_update_sql(metaseqId, evidence)
            updater.increment_update_count(isIndel)
            matched = True
            if args.debug: 
                warning(metaseqId, "- matched -", matchedAlleles, evidence, file=updater.lfh(), flush=True)
            break # no need to look for more matches

    if not matched:
        updater.buffer_update_sql(metaseqId, {}) # create a place holder to avoid future lookups
        updater.increment_not_matched_count()
        if args.debug:
            warning(metaseqId, "- not matched", file=updater.lfh(), flush=True)


def update_variant_records_from_vcf():
    ''' lookup and update variant records from a VCF file 
    assuming the load by file was called by a plugin, so this variant has already been
    verified to be new to the resource; no need to check alternative metaseq IDs'''

    cupdater = CADDUpdater(args.logFile, args.databaseDir)

    database = Database(args.gusConfigFile)
    database.connect()

    lineCount = 0
    with database.cursor() as updateCursor, \
      open(args.vcfFile, 'r') as fh:
        try:
            for line in fh:
                if line.startswith("#"):
                    continue

                lineCount = lineCount + 1
                entry = VcfEntryParser(line.rstrip())
                
                refAllele = entry.get('ref')
                altAllele = entry.get('alt')
                chrom = xstr(entry.get('chrom'))
                if chrom == 'MT':
                    chrom = 'M'
                position = int(entry.get('pos'))
                metaseqId = ':'.join((chrom, xstr(position), refAllele, altAllele))

                record = {"metaseq_id" : metaseqId}  # mimic "record"

                if len(refAllele) > 1 or len(altAllele) > 1: # only doing SNVs
                    update_variant_record(record, cupdater, INDEL)
                else:
                    update_variant_record(record, cupdater, SNV)

                if lineCount % args.commitAfter == 0:
                    warning("Processed:", lineCount, "- SNVs:",
                            cupdater.get_update_count(SNV), "- INDELS:", cupdater.get_update_count(INDEL),
                            " - Not Matched:", cupdater.get_not_matched_count(),
                            file=cupdater.lfh(), flush=True)

                    updateCursor.execute(cupdater.sql_buffer_str())

                    if args.commit:
                        database.commit()
                    else:
                        database.rollback()
                    cupdater.clear_update_sql()

            if cupdater.buffered_variant_count() > 0:  # trailing
                updateCursor.execute(cupdater.sql_buffer_str())

                if args.commit:
                    database.commit()
                else:
                    database.rollback()

            warning("DONE - Updated SNVs:", cupdater.get_update_count(SNV), "- Updated INDELS:",
                    cupdater.get_update_count(INDEL), 
                    "- Not Matched:", cupdater.get_not_matched_count(), file=cupdater.lfh(), flush=True)

            cupdater.close_lfh()

        except Exception as e:
            warning(e, entry, file=cupdater.lfh(), flush=True)
            database.rollback()
            database.close()
            print("FAIL", file=sys.stdout)
            raise

        
           
def update_variant_records(chromosome):
    chrLabel = 'chr' + xstr(chromosome)

    logFileName = path.join(args.logFilePath, chrLabel + '.log')
    cupdater = CADDUpdater(logFileName, args.databaseDir)
    cupdater.setChrm(chrLabel)

    selectSQL = "SELECT metaseq_id, cadd_scores FROM Variant_" + chrLabel;

    database = Database(args.gusConfigFile)
    database.connect()

    lineCount = 0
    updateCount = 0
    updateIndelCount = 0
    skipCount = 0

    with database.cursor("RealDictCursor") as selectCursor, \
        database.cursor() as updateCursor:
        try:
            warning("Fetching", chrLabel, "variants", file=cupdater.lfh(), flush=True)
            selectCursor.execute(selectSQL)
            warning("DONE - Fetching", file=cupdater.lfh(), flush=True)

            for record in selectCursor:
                if args.debug and args.veryVerbose:
                    warning(record, file=cupdater.lfh(), flush=True)

                if record['cadd_scores'] is not None:
                    if args.debug and args.veryVerbose:
                        warning("Skipping", record['metaseq_id'], file=cupdater.lfh(), flush=True)
                    skipCount = skipCount + 1
                    continue

                lineCount = lineCount + 1

                metaseqId = record['metaseq_id']

                chrm, position, refAllele, altAllele = metaseqId.split(':')

                if len(refAllele) > 1 or len(altAllele) > 1: # only doing SNVs
                    update_variant_record(record, cupdater, INDEL)
                else:
                    update_variant_record(record, cupdater, SNV)

                if cupdater.get_total_update_count() % args.commitAfter == 0 and cupdater.buffered_variant_count() > 0:
                    if args.commit:
                        if args.debug:
                            warning("Starting Update", file=cupdater.lfh(), flush=True)

                        updateCursor.execute(cupdater.sql_buffer_str())

                        if args.debug:
                            warning("Done", file=cupdater.lfh(), flush=True)

                        cupdater.clear_update_sql()
                        database.commit()
                    else:
                        database.rollback()

                    warning(metaseqId, "- Processed:", lineCount, "- SNVs:",
                            cupdater.get_update_count(SNV), "- INDELS:", cupdater.get_update_count(INDEL),
                            "- Skipped:", skipCount, " - Not Matched:", cupdater.get_not_matched_count(),
                            file=cupdater.lfh(), flush=True)

            if cupdater.buffered_variant_count() > 0:
                updateCursor.execute(cupdater.sql_buffer_str())
                if args.commit: # trailing
                    database.commit()
                else:
                    database.rollback()

            warning("DONE - Updated SNVs:", cupdater.get_update_count(SNV), "- Updated INDELS:",
                    cupdater.get_update_count(INDEL), "- Skipped", skipCount,
                    "- Not Matched:", cupdater.get_not_matched_count(), file=cupdater.lfh(), flush=True)
            cupdater.close_lfh()

        except Exception as e:
            warning(e, file=cupdater.lfh(), flush=True)
            if args.commit:
                database.commit()
            else:
                database.rollback()
            database.close()
            raise

    
    database.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="update CADD scores in DB specified chromosome")
    parser.add_argument('-d', '--databaseDir', help="directory containg CADD database +tabindex files", required=True)
    parser.add_argument('--vcfFile', help="if file is specified, updates only the listed variants; otherwise updates all variants in the database for the specified chromosome; expects full path to a VCF file")
    parser.add_argument('-c', '--chr', help="chromosome; comma separated list of one or more chromosomes or 'all'; ignored if --file is specified", default='all')
    parser.add_argument('--maxWorkers', type=int, default=5)
    parser.add_argument('--commitAfter', type=int, default=500)
    parser.add_argument('--logFilePath', help='generate formulaic log files and store in the specified path')
    parser.add_argument('--logFile', help='specify full path to log file')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--veryVerbose', action='store_true')
    parser.add_argument('--gusConfigFile',
                        help="GUS config file. If not provided, assumes default: $GUS_HOME/conf/gus.config")
    parser.add_argument('--commit', action='store_true')
    args = parser.parse_args()
    
    INDEL = True
    SNV = False
 
    if args.vcfFile:
        update_variant_records_from_vcf()

    else:
        chrList = args.chr.split(',') if args.chr != 'all' \
          else [c.value for c in Human]

        random.shuffle(chrList) # so that not all large chrms are done at once if all is selected

        if len(chrList) == 1:
            update_variant_records(chrList[0])

        else:
            with ProcessPoolExecutor(args.maxWorkers) as executor:
                for c in chrList:
                    warning("Create and start thread for chromosome:", xstr(c))
                    executor.submit(update_variant_records, c)


    print("SUCCESS", file=sys.stdout)
