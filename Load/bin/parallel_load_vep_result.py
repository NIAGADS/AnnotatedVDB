#!/usr/bin/env python
# pylint: disable=invalid-name,not-an-iterable,unused-import,too-many-locals
"""
Loads AnnotatedVDB from JSON output from running VEP on dbSNP
Loads each chromosome in parallel
"""

from __future__ import print_function

import argparse
import gzip
import mmap
import json
import sys
import csv
import traceback

from datetime import datetime
from os import path
import multiprocessing
from concurrent.futures import ProcessPoolExecutor

from GenomicsDBData.Util.utils import xstr, warning, print_dict, print_args
from GenomicsDBData.Util.postgres_dbi import Database, raise_pg_exception

from AnnotatedVDB.Util.variant_loader import VEPLoader as VariantLoader
from AnnotatedVDB.Util.chromosomes import Human


def parse_chromosome_map():
    """! parse chromosome map
    @returns           dict representation of chromosome map if mapping file was provided
    """

    if args.verbose:
        warning("Loading chromosome map from:", args.chromosomeMap)

    if args.chromosomeMap is None:
        return None

    cMap = {}
    with open(args.chromosomeMap, 'r') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            # source_id	chromosome	chromosome_order_num	length
            key = row['source_id']
            value = row['chromosome'].replace('chr', '')
            cMap[key] = value

    return cMap


def load_annotation(chromosome):
    """! parse over a JSON file, extract position, frequencies,
    ids, and ADSP-ranked most severe consequence; bulk load using COPY """

    lfn = "annotatedvdb_dbsnp_chr_" + xstr(chromosome) + '.log';
    warning("Logging to", lfn)
    
    loader = VariantLoader('dbSNP', logFileName=lfn, 
                           verbose=args.verbose, debug=args.debug)
    if args.verbose:
        loader.log(("Parameters:", print_dict(vars(args), pretty=True)), prefix="INFO")

    database = Database(args.gusConfigFile)
    database.connect()
    
    loader.set_algorithm_invocation('parallel_load_vep_result', 'chr' + xstr(chromosome) + '|' + print_args(args, False))
    loader.log('Algorithm Invocation Id = ' + xstr(loader.alg_invocation_id), prefix="INFO")
    
    loader.initialize_pk_generator(args.genomeBuild, args.seqrepoProxypath)
    loader.initialize_vep_parser(args.rankingFie, True, args.verbose)
    loader.initialize_bin_indexer(args.gusConfigFile)
    loader.initialize_copy_sql() # use default copy fields
    loader.set_chromosome_map(chrmMap)

    fileName = path.join(args.dir, 
                         'chr' + xstr(chromosome) + "." + args.extension)
    loader.log('Parsing ' + fileName, prefix="INFO")

    resume = args.resumeAfter is None # false if need to skip lines
    if not resume:
        loader.log(("--resumeAfter flag specified; Finding skip until point", args.resumeAfter), prefix="INFO")
        loader.set_resume_after_variant(args.resumeAfter)
    with database.cursor() as cursor:
        loader.set_cursor(cursor)
        lineCount = 1
        with open(fileName, 'r') as fhandle:
            mappedFile = mmap.mmap(fhandle.fileno(), 0, prot=mmap.PROT_READ) # put file in swap
            with gzip.GzipFile(mode='r', fileobj=mappedFile) as gfh:
                for line in gfh:
                    try: 
                        if lineCount == 1 \
                            or lineCount % args.commitAfter == 0: 
                            if args.debug:
                                loader.log('Processing new copy object')
                            tstart = datetime.now()
                            
                        loader.parse_result(line.rstrip())

                        if lineCount % args.logAfter == 0 \
                            and lineCount % args.commitAfter != 0:
                            loader.log((lineCount, "lines (", 
                                        loader.get_count('variant'), ") variants"), 
                                       prefix="PARSED")
                      
                        if lineCount % args.commitAfter == 0:
                            if args.debug:
                                tendw = datetime.now()
                                message = 'Copy object prepared in ' + str(tendw - tstart) + '; ' + \
                                    str(loader.copy_buffer(sizeOnly=True)) + ' bytes; transfering to database'
                                loader.log(message, prefix="DEBUG") 
                                
                            loader.insert_variants()

                            message = '{:,}'.format(loader.get_count('variant')) + " variants"
                            messagePrefix = "COMMITTED"
                            if args.commit:
                                database.commit()

                            else:
                                database.rollback()
                                messagePrefix = "PARSED"
                                message += " -- rolling back"

                            if lineCount % args.logAfter == 0:
                                message += "; up to = " + loader.get_current_variant_id()
                                loader.log(message, prefix=messagePrefix)

                                if args.debug:
                                    tend = datetime.now()
                                    loader.log('Database copy time: ' + str(tend - tendw), prefix="DEBUG")
                                    loader.log('        Total time: ' + str(tend - tstart), prefix="DEBUG")

                            if args.test:
                                loader.log("DONE - TEST COMPLETE" , prefix="WARNING")
                                sys.exit("EXITING - TEST COMPLETE")

                    except Exception as err:
                        loader.log("Problem parsing variant: " + loader.get_current_variant_id(), prefix="ERROR")
                        loader.log((lineCount, ":", print_dict(json.loads(line))), prefix="LINE")
                        loader.log(str(err), prefix="ERROR")
                        raise(err)

            # final commit / leftovers
            loader.insert_variants()

            message = '{:,}'.format(loader.get_count('variant')) + " variants"
            messagePrefix = "COMMITTED"
            if args.commit:
                database.commit()
            else:
                database.rollback()
                messagePrefix = "PARSED"
                message += " -- rolling back"
            message += "; up to = " + loader.get_current_variant_id()
            loader.log(message, prefix=messagePrefix)
            
            # summarize new consequences
            loader.log(("Counsequences added during load:", 
                      loader.vep_parser().get_added_conseq_summary()), prefix="INFO")

            mappedFile.close()
            
    database.close()
    loader.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='load AnnotatedDB from  JSON output of VEP against dbSNP')
    parser.add_argument('-d', '--dir',
                        help="directory containing VEP results", required=True)
    parser.add_argument('-e', '--extension', required=True,
                        help="file extension (e.g., json.gz)")
    parser.add_argument('-g', '--genomeBuild', default='GRCh38', help="genome build: GRCh37 or GRCh38")
    parser.add_argument('-s', '--seqrepoProxyPath', required=True,
                        help="full path to local SeqRepo file repository")
    parser.add_argument('-m', '--chromosomeMap', required=False,
                        help="chromosome map")
    parser.add_argument('-r', '--rankingFile', required=True,
                        help="full path to ADSP VEP consequence ranking file")
    parser.add_argument('--commit', action='store_true', help="run in commit mode", required=False)
    parser.add_argument('--gusConfigFile',
                        '--full path to gus config file, else assumes $GUS_HOME/config/gus.config')
    parser.add_argument('--test', action='store_true',
                        help="load 'commitAfter' rows as test")
    parser.add_argument('--resumeAfter',
                        help="refSnpId after which to resume load (log lists lasts committed refSNP)")
    parser.add_argument('-c', '--chr', required=True,
                        help="comma separated list of one or more chromosomes to load, e.g., 1, 2, M, X or `all`")
    parser.add_argument('--commitAfter', type=int, default=500,
                        help="commit after specified inserts")
    parser.add_argument('--maxWorkers', default=10, type=int)
    parser.add_argument('--logAfter', type=int,
                        help="number of inserts to log after completion; will work best if factor/multiple of commitAfter")
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--debug',  action='store_true',
                        help="log database copy time / may print development debug statements")
    args = parser.parse_args()

    if not args.logAfter:
        args.logAfter = args.commitAfter

    chrmMap = parse_chromosome_map()

    chrList = args.chr.split(',') if args.chr != 'all' \
      else [c.value for c in Human]

    numChrs = len(chrList)
    with ProcessPoolExecutor(args.maxWorkers) as executor:
        if numChrs > 1:
            for c in chrList:
                warning("Create and start thread for chromosome:", xstr(c))
                executor.submit(load_annotation, c)
        else: # for debugging
            load_annotation(chrList[0])
