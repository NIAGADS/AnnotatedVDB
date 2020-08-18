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

import multiprocessing as mp

import threading
from concurrent.futures import ThreadPoolExecutor

from os import path

from CBILDataCommon.Util.utils import warning, xstr, die, qw, execute_cmd, xstrN
from CBILDataCommon.Util.postgres_dbi import Database
from CBILDataCommon.Util.exceptions import print_exception

from AnnotatedVDB.BinIndex.bin_index import BinIndex
from AnnotatedVDB.Util.algorithm_invocation import AlgorithmInvocation
from AnnotatedVDB.Util.vcf_parser import VcfEntryParser

from AnnotatedVDB.Util.chromosomes import Human


def worker(chunk):
 # `chunk` will be a list of CSV rows all with the same name column
    # replace this with your real computation
    # print(chunk)
    return len(chunk)  


def variant_is_missing(database, chromosome, refSnpId):
    ''' lookup refsnp id '''
    sql = "SELECT ref_snp_id FROM Variant WHERE chromosome = %s AND ref_snp_id = %s"
    with database.cursor() as cursor:
        cursor.execute(sql, (chromosome, refSnpId))
        result = cursor.fetchone()

        return result is None



def json2str(value):
    ''' convert JSON to str; accounting for None values '''
    if value is None:
        return 'NULL'
    else:
        return json.dumps(value)


def load_annotation(chromosome):
    ''' extract basic SNP information from VCF line; check against AnnotatedDB; if missing
    load '''

    pool = mp.Pool()
    indexer = BinIndex(args.gusConfigFile, verbose=False)

    database = Database(args.gusConfigFile)
    database.connect()


    fname = '00-All.' + xstr(chromosome) + '.vcf.gz'
    logFname = path.join(args.logDir, fname + '.log')
    fname = path.join(args.dir, fname)

    lineCount = 0
    variantCount = 0

    warning("Parsing", fname)
    warning("Logging:", logFname)
    with database.cursor() as cursor:
        with open(fname, 'r') as fhandle, open(logFname, 'w') as lfh:
            warning("Parsing", fname, file=lfh, flush=True)
            mappedFile = mmap.mmap(fhandle.fileno(), 0, prot=mmap.PROT_READ) # put file in swap
            with gzip.GzipFile(mode='r', fileobj=mappedFile) as gfh:
                for line in gfh:
                    lineCount = lineCount + 1
                    vepResult = {} # will just be the input string, so it matches the annotated variants

                    if line.startswith('#'):
                        continue

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
                        if variant_is_missing(database, 'chr' + chrom, refSnpId):
                            warning(refSnpId, "- MISSING -- LOADING", file=lfh, flush=True)
                            for alt in altAllele:
                                if alt == '.': # most frequently skipped by VEP for this reason
                                    alt = '?'
                                variantCount = variantCount + 1
                                metaseqId = ':'.join(('chr' + chrom, xstr(position), refAllele, alt))
                                positionEnd = entry.infer_variant_end_location(alt)

                                binIndex = indexer.find_bin_index(chrom, position, positionEnd)
                                cursor.execute(INSERT_SQL, 
                                                   ('chr' + chrom,
                                                    xstr(position),
                                                    isMultiAllelic,
                                                    binIndex,
                                                    refSnpId,
                                                    metaseqId,
                                                    json.dumps(vepResult),
                                                    algInvocId))

                                if args.commit:
                                    database.commit()
                                else:
                                    database.rollback()

                        if lineCount % 25000 == 0:
                            warning("Parsed", lineCount, "lines.", file=lfh, flush=True)
                            warning("Loaded", variantCount, "missing variants", file=lfh, flush=True)

                    except Exception:
                        warning("ERROR parsing variant", refSnpId, file=lfh, flush=True)
                        warning(lineCount, ":", line, file=lfh, flush=True)
                        raise

         
            if not args.commit:
                database.rollback()
                warning("DONE -- rolling back")

            mappedFile.close()

    database.close()
    indexer.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='load dbSNP variants missin from VEP output int oAnnotatedDB')
    parser.add_argument('-d', '--dir',
                        help="directory containing VCF files", required=True)
    parser.add_argument('--commit', action='store_true', help="run in (auto)commit mode", required=False)
    parser.add_argument('--logDir', help="log directory", required=True)
    parser.add_argument('--gusConfigFile',
                        '--full path to gus config file, else assumes $GUS_HOME/config/gus.config')
    parser.add_argument('--test', action='store_true',
                        help="load 'commitAfter' rows as test")
    parser.add_argument('--maxWorkers', default=5, type=int)
    parser.add_argument('-c', '--chr', required=True,
                        help="comma separated list of chromosomes to load, or 'all'")
  
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()


    N_CHUNKS = 10
    INSERT_SQL = """INSERT INTO Variant 
(chromosome, location, is_multi_allelic, bin_index, ref_snp_id, metaseq_id, vep_output, row_algorithm_id)
VALUES (%s, %s, %s, %s, %s, %s, %s, %s)"""

    algInvocation = AlgorithmInvocation('load_missing_dbsnp_from_vcf.py', json.dumps(vars(args)), args.commit)
    algInvocId = xstr(algInvocation.getAlgorithmInvocationId())
    algInvocation.close()

    warning("Algorithm Invocation ID", algInvocId)
   
    chrList = args.chr.split(',') if args.chr != 'all' \
      else [c.value for c in Human]

    with ThreadPoolExecutor(max_workers=args.maxWorkers) as executor:
        for c in chrList:
            warning("Create and start thread for chromosome:", xstr(c))
            executor.submit(load_annotation, c)
