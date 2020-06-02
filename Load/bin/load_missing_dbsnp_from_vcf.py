#!/usr/bin/env python
# pylint: disable=invalid-name,not-an-iterable,unused-import,too-many-locals
"""
Iterates over VCF files and loads any missing variants that were not annotated by VEP
Run after load_dbsnp_from_vep_result.py
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

from CBILDataCommon.Util.utils import warning, xstr, die, qw, execute_cmd, xstrN
from CBILDataCommon.Util.postgres_dbi import Database
from CBILDataCommon.Util.exceptions import print_exception

from AnnotatedVDB.BinIndex.bin_index import BinIndex
from AnnotatedVDB.Util.algorithm_invocation import AlgorithmInvocation
from AnnotatedVDB.Util.vcf_parser import VcfEntryParser


def variant_is_missing(chromosome, refSnpId):
    ''' lookup refsnp id '''
    sql = "SELECT ref_snp_id FROM Variant WHERE chromosome = ? and ref_snp_id = ?"
    with database.cursor() as cursor:
        cursor.execute(sql, (chromosome, refSnpId))
        result = cursor.fetchone()
        if result is None:
            return True
        else: 
            return False


def json2str(value):
    ''' convert JSON to str; accounting for None values '''
    if value is None:
        return 'NULL'
    else:
        return json.dumps(value)


def load_annotation(chromosome):
    ''' extract basic SNP information from VCF line; check against AnnotatedDB; if missing
    load '''

    fname = '00-All.' + xstr(chromosome) + '.vcf.gz'
    fname = path.join(args.dir, fname)
    lineCount = 0
    variantCount = 0
    skipCount = 0

    with database.cursor() as cursor:
        with open(fname, 'r') as fhandle:
            mappedFile = mmap.mmap(fhandle.fileno(), 0, prot=mmap.PROT_READ) # put file in swap
            with gzip.GzipFile(mode='r', fileobj=mappedFile) as gfh:
                for line in gfh:
                    lineCount = lineCount + 1
                    vepResult = {} # will just be the input string, so it matches the annotated variants
                    entry = VcfEntryParser(line.rstrip())

                    # there are json formatting issues w/the input str
                    # so replace w/the parsed entry; which is now a dict
                    vepResult['input'] = entry.get_entry()

                    refSnpId = entry.get_refsnp()

                    refAllele = entry.get('ref')
                    altAllele = entry.get('alt').split(',')

                    chrom = xstr(entry.get('chrom'))
                    position = int(entry.get('pos'))

                    isMultiAllelic = len(altAllele) > 1

                    try:
                        if variant_is_missing('chr' + chrom, refSnpId):
                            warning(refSnpId, "- MISSING -- LOADING")
                            for alt in altAllele:
                                variantCount = variantCount + 1
                                metaseqId = ':'.join(('chr' + chrom, xstr(position), refAllele, alt))
                                positionEnd = entry.infer_variant_end_location(alt)

                                binIndex = indexer.find_bin_index(chrom, position, positionEnd)
                                cursor.execute(INSERT_SQL, 
                                                   ('chr' + xstr(chrom),
                                                    xstr(position),
                                                    isMultiAllelic,
                                                    binIndex,
                                                    refSnpId,
                                                    metaseqId,
                                                    json.dumps(vepResult),
                                                    algInvocId))
                            die("DEBUG: " + refSnpId + " - LOADED")
                        else:
                            warning("DEBUG:", refSnpId + " - NOT MISSING")

                        if lineCount % 25000 == 0:
                            warning("Parsed", lineCount, "lines.")
                            warning("Loaded", variantCount, "missing variants")

                    except Exception:
                        warning("ERROR parsing variant", refSnpId)
                        warning(lineCount, ":", line)
                        raise

         
            if not args.commit:
                database.rollback()
                warning("DONE -- rolling back")

            mappedFile.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='load dbSNP variants missin from VEP output int oAnnotatedDB')
    parser.add_argument('-d', '--dir',
                        help="directory containing VCF files", required=True)
    parser.add_argument('--commit', action='store_true', help="run in (auto)commit mode", required=False)
    parser.add_argument('--gusConfigFile',
                        '--full path to gus config file, else assumes $GUS_HOME/config/gus.config')
    parser.add_argument('--test', action='store_true',
                        help="load 'commitAfter' rows as test")
    parser.add_argument('-c', '--chr', required=True,
                        help="chromosome to load, e.g., 1, 2, M, X")
  
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()

    INSERT_SQL = """INSERT INTO Variant 
(chromosome, location, is_multi_allelic, bin_index, ref_snp_id, metaseq_id, vep_output, row_algorithm_invocation_id)
VALUES (%s, %s, %s, %s, %s, %s, %s, %s)"""

    algInvocation = AlgorithmInvocation('load_missing_dbsnp_from_vcf.py', json.dumps(vars(args)), args.commit)
    algInvocId = xstr(algInvocation.getAlgorithmInvocationId())
    algInvocation.close()

    warning("Algorithm Invocation ID", algInvocId)

    indexer = BinIndex(args.gusConfigFile, verbose=False)

    database = Database(args.gusConfigFile)
    database.connect()
    if args.commit: # b/c adding a few hundred out of 10s of millions of variants
        database.autocommit()

    load_annotation(args.chr)

    database.close()
    indexer.close()
