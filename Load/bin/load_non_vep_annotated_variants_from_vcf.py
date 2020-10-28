#!/usr/bin/env python
# pylint: disable=invalid-name,not-an-iterable,unused-import,too-many-locals
"""
Iterates over VCF files and loads any missing variants that were not annotated by VEP
Run after load_vep_result.py
"""
from __future__ import print_function

import argparse
import gzip
import mmap
import json
import sys
import io
import copy
from datetime import datetime

from os import path

from CBILDataCommon.Util.utils import warning, xstr, die, qw, execute_cmd, xstrN, truncate
from CBILDataCommon.Util.postgres_dbi import Database
from CBILDataCommon.Util.exceptions import print_exception

from AnnotatedVDB.BinIndex.bin_index import BinIndex
from AnnotatedVDB.Util.algorithm_invocation import AlgorithmInvocation
from AnnotatedVDB.Util.vcf_parser import VcfEntryParser

from AnnotatedVDB.Util.chromosomes import Human


def variant_is_missing(database, chromosome, metaseq_id):
    ''' lookup variant id '''
    sql = "SELECT metaseq_id FROM  Variant WHERE LEFT(metaseq_id, 50) = LEFT(%s, 50)"
    sql = sql + " AND chromosome = %s"
    sql = sql + " AND metaseq_id = %s"

    with database.cursor() as cursor:
        cursor.execute(sql, (metaseq_id, chromosome,  metaseq_id)) 
        for record in cursor:
            return False # not missing
    return True # missing


def json2str(value):
    ''' convert JSON to str; accounting for None values '''
    if value is None:
        return 'NULL'
    else:
        return json.dumps(value)


def load_annotation():
    ''' extract basic SNP information from VCF line; check against AnnotatedDB; if missing
    load '''

    indexer = BinIndex(args.gusConfigFile, verbose=False)

    database = Database(args.gusConfigFile)
    database.connect()

    lineCount = 0
    variantCount = 0
    with database.cursor() as cursor:
        with open(args.vcfFile, 'r') as fh:
            with open(args.logFileName, 'w') as lfh:
                warning("Parsing", args.vcfFile, file=lfh, flush=True)
                for line in fh:
                    lineCount = lineCount + 1
                    vepResult = {} # will just be the input string, so it matches the annotated variants

                    if line.startswith('#'):
                        continue

                    entry = VcfEntryParser(line.rstrip())

                    # there are json formatting issues w/the input str
                    # so replace w/the parsed entry; which is now a dict
                    vepResult['input'] = entry.get_entry()

                    refAllele = entry.get('ref')
                    altAllele = entry.get('alt') # assuming not multiallelic

                    if refAllele == '0': # happens sometimes
                        refAllele = '?'
                    if altAllele == '0':
                        altAllele = '?'

                    # truncatedRef = truncate(refAllele, 20)
                    # truncatedAlt = truncate(altAllele, 20)

                    chrom = xstr(entry.get('chrom'))
                    if chrom == 'MT':
                        chrom = 'M'
                    position = int(entry.get('pos'))


                    # metaseqId = ':'.join((chrom, xstr(position), truncatedRef, truncatedAlt))
                    isMultiAllelic = False # assuming b/c of how the VCF files were generated
                    metaseqId = ':'.join((chrom, xstr(position), refAllele, altAllele))

                    try:
                        if variant_is_missing(database, 'chr' + chrom, metaseqId):
                            warning(metaseqId, "- MISSING -- LOADING", file=lfh, flush=True)
                            variantCount = variantCount + 1

                            positionEnd = entry.infer_variant_end_location(altAllele)
                            binIndex = indexer.find_bin_index(chrom, position, positionEnd)

                            cursor.execute(INSERT_SQL, 
                                               ('chr' + chrom,
                                                position,
                                                isMultiAllelic,
                                                binIndex,
                                                metaseqId,
                                                json.dumps(vepResult),
                                                algInvocId))

                            if args.commit:
                                database.commit()
                            else:
                                database.rollback()

                        if lineCount % 50 == 0:
                            warning("Parsed", lineCount, "lines.", file=lfh, flush=True)
                            warning("Loaded", variantCount, "missing variants", file=lfh, flush=True)

                    except Exception:
                        warning("ERROR parsing variant", metaseqId, file=lfh, flush=True)
                        warning(lineCount, ":", line, file=lfh, flush=True)
                        print("FAIL", file=sys.stdout)
                        raise

                warning("DONE - Loaded", variantCount, "missing variants", file=lfh, flush=True)
    database.close()
    indexer.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='load variants missing from VEP output into AnnotatedDB')
    parser.add_argument('-v', '--vcfFile',
                        help="input file (vcf)", required=True)
    parser.add_argument('--commit', action='store_true', help="run in (auto)commit mode", required=False)
    parser.add_argument('--logFileName', help="log file name", required=True)
    parser.add_argument('--gusConfigFile',
                        '--full path to gus config file, else assumes $GUS_HOME/config/gus.config')
  
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()

    INSERT_SQL = """INSERT INTO Variant 
(chromosome, location, is_multi_allelic, bin_index, metaseq_id, vep_output, row_algorithm_id)
VALUES (%s, %s, %s, %s, %s, %s, %s)"""

    algInvocation = AlgorithmInvocation('load_non_vep_annotated_variants_from_vcf.py', json.dumps(vars(args)), args.commit, args.gusConfigFile)
    algInvocId = xstr(algInvocation.getAlgorithmInvocationId())
    algInvocation.close()

    warning("Algorithm Invocation ID", algInvocId)
   
    load_annotation()

    print(algInvocId, file=sys.stdout)
         
