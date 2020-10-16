#!/usr/bin/env python
# pylint: disable=invalid-name,not-an-iterable,unused-import,too-many-locals
"""
Loads AnnotatedVDB from JSON output from running VEP on a non-dbSNP VCF
NOTE: assumes variants were not mapped to a refSNP
NOTE: to skip check against DB to avoid duplicates and to speed up (NOT RECOMMENDED)
use the --skipVerification flag
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

from AnnotatedVDB.BinIndex.bin_index import BinIndex
from AnnotatedVDB.Util.algorithm_invocation import AlgorithmInvocation
from AnnotatedVDB.Util.vcf_parser import VcfEntryParser
from AnnotatedVDB.Util.vep_parser import VepJsonParser, CONSEQUENCE_TYPES

def duplicate(metaseqId, dcursor):
    ''' test against db '''
    if args.skipVerification:
        return False
    sql = "SELECT * FROM find_variant_by_metaseq_id(%s)"
    dcursor.execute(sql, (metaseqId,))
    for record in dcursor:
        return True # if there is a hit, there is a duplicate
    return False # if not, no duplicate


def get_frequencies():
    ''' extract frequencies'''

    frequencies = vepParser.get_frequencies()
    if frequencies is None:
        return None
    return frequencies['values']


def json2str(value):
    ''' convert JSON to str; accounting for None values '''
    if value is None:
        return 'NULL'
    else:
        return json.dumps(value)


def load_annotation():
    ''' parse over a JSON file, extract position, frequencies,
    ids, and ADSP-ranked most severe consequence; bulk load using COPY '''

    fname = args.inputFile
    lineCount = 0
    variantCount = 0
    skipCount = 0

    warning("Parsing variants from:", fname) # should print to plugin log
    with database.cursor() as cursor, database.cursor("RealDictCursor") as dcursor:
        copyObj = io.StringIO()
        with open(fname, 'r') as fh:
            with open(args.logFile, 'w') as lfh:
                warning("Parsing variants from:", fname, file=lfh, flush=True)
                for line in fh:
                    lineCount = lineCount + 1
                    vepResult = json.loads(line.rstrip())

                    vepParser.set_annotation(copy.deepcopy(vepResult))
                    vepInputStr = vepParser.get('input')
                    entry = VcfEntryParser(vepInputStr)

                    # there may be json formatting issues w/the input str
                    # so replace w/the parsed entry; which is now a dict
                    vepResult['input'] = entry.get_entry()

                    if lineCount == 1 or variantCount % 5000 == 0:
                        warning('Processing new copy object', file=lfh, flush=True)
                        tstart = datetime.now()

                    refAllele = entry.get('ref')
                    altAllele = entry.get('alt')  # assumes no multialleic variants

                    if refAllele == '0': # happens sometimes
                        refAllele = '?'
                    if altAllele == '0':
                        altAllele = '?'

                    truncatedRef = truncate(refAllele, 20)
                    truncatedAlt = truncate(altAllele, 20)
                  
                    chrom = xstr(entry.get('chrom'))
                    if chrom == 'MT':
                        chrom = 'M'

                    position = int(entry.get('pos'))

                    try:
                        metaseqId = ':'.join((chrom, xstr(position), truncatedRef, truncatedAlt))

                        if duplicate(metaseqId, dcursor): 
                            warning("SKIPPING:",metaseqId, "- already loaded.", file=lfh, flush=True)
                            continue

                        vepParser.set('ref_allele', refAllele)
                        vepParser.set('alt_allele', altAllele)
                     
                        frequencies = get_frequencies()
                        vepParser.adsp_rank_and_sort_consequences()

                        # for each allele
                        variantCount = variantCount + 1

                        positionEnd = entry.infer_variant_end_location(altAllele)
                        binIndex = indexer.find_bin_index(chrom, position, positionEnd)

                        # NOTE: VEP uses normalized alleles to indicate variant_allele
                        nRef, nAlt = entry.normalize_alleles(refAllele, altAllele, snvDivMinus=True)
                        alleleFreq = None if frequencies is None \
                          else get_allele_frequencies(nAlt, frequencies)

                        msConseq = get_most_severe_consequence(nAlt)
                        valueStr = '#'.join((
                          'chr' + xstr(chrom),
                          xstr(position),
                          binIndex,
                          metaseqId,
                          json2str(alleleFreq),
                          json2str(msConseq),
                          json2str(get_adsp_ranked_allele_consequences(nAlt)),
                          json.dumps(vepResult),
                          algInvocId
                          ))


                        copyObj.write(valueStr + '\n')

                        if variantCount % 5000 == 0:
                            warning("FOUND", variantCount, " new variants", file=lfh, flush=True)
                            tendw = datetime.now()
                            warning('Copy object prepared in ' + str(tendw - tstart) + '; ' +
                                    str(copyObj.tell()) + ' bytes; transfering to database',
                                    file=lfh, flush=True)
                            copyObj.seek(0)
                            cursor.copy_from(copyObj, 'variant', sep='#', null="NULL",
                                             columns=VARIANT_COLUMNS)

                            message = '{:,}'.format(variantCount)
                            if args.commit:
                                database.commit()
                                message = "COMMITTED " + message

                            else:
                                database.rollback()
                                message = "PARSED " + message + " -- rolling back"


                            warning(message, "; up to = ", metaseqId, file=lfh, flush=True)

                            tend = datetime.now()
                            warning('Database copy time: ' + str(tend - tendw), file=lfh, flush=True)
                            warning('        Total time: ' + str(tend - tstart), file=lfh, flush=True)

                            copyObj = io.StringIO() # reset io string

                    except Exception:
                        warning("ERROR parsing variant on line", line, ' - ', metaseqId,
                                file=lfh, flush=True)
                        print("FAIL", file=sys.stdout)
                        raise

                # final commit / leftovers
                copyObj.seek(0)
                cursor.copy_from(copyObj, 'variant', sep='#', null="NULL", columns=VARIANT_COLUMNS)
                message = '{:,}'.format(variantCount)

                if args.commit:
                    database.commit()
                    message = "DONE - COMMITTED " + message
                else:
                    database.rollback()
                    message = "DONE - PARSED " + message + " -- rolling back"
                    
                warning(message, file=lfh, flush=True)



def get_most_severe_consequence(allele):
    ''' retrieve most severe consequence from the VEP JSON Parser,
    for the specified allele; returns None if no consequences are found
    return first hit among transcript, then regulatory feature, then intergenic
    consequences
    '''

    for ct in CONSEQUENCE_TYPES:
        msConseq = get_allele_consequences(allele, vepParser.get(ct + '_consequences'))
        if msConseq is not None:
            return msConseq[0]

    return None


def get_adsp_ranked_allele_consequences(allele):
    ''' get dict of all ADSP ranked allele consequences '''
    adspConseq = {}
    for ctype in CONSEQUENCE_TYPES:
        ctypeKey = ctype + '_consequences'
        cq = get_allele_consequences(allele, vepParser.get(ctypeKey))
        if cq is not None:
            adspConseq[ctypeKey] = cq

    return adspConseq


def get_allele_consequences(allele, conseq):
    ''' check to see if conseq dict contains values for the specified allele,
    if so return first value associated w/the allele as the most
    serious; assumes ADSP ranking and re-sort has occurred '''

    if conseq is None:
        return None

    if allele in conseq:
        return conseq[allele]

    return None


def get_allele_frequencies(allele, frequencies):
    ''' given an allele and frequency dict,
    retrieve and return allele-specific frequencies '''

    if frequencies is None:
        return None

    if 'values' in frequencies:
        if allele in frequencies['values']:
            return frequencies['values'][allele]

    if allele in frequencies:
        return frequencies[allele]

    return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='load AnnotatedDB from  JSON output of VEP against a non-dbSNP VCF; assumes no multiallelic variants')
    parser.add_argument('-i', '--inputFile',
                        help="full path to file containing VEP results", required=True)
    parser.add_argument('--logFile', required=True)
    parser.add_argument('-r', '--rankingFile', required=True,
                        help="full path to ADSP VEP consequence ranking file")
    parser.add_argument('--commit', action='store_true', help="run in commit mode", required=False)
    parser.add_argument('--gusConfigFile',
                        help="full path to gus config file, else assumes $GUS_HOME/config/gus.config")
    parser.add_argument('--skipVerification', help="skip check against DB (for avoiding duplications")
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()


    VARIANT_COLUMNS = qw('chromosome location bin_index metaseq_id allele_frequencies adsp_most_severe_consequence adsp_ranked_consequences vep_output row_algorithm_id', returnTuple=True)

    algInvocation = AlgorithmInvocation('load_vep_result.py', json.dumps(vars(args)), args.commit, args.gusConfigFile)
    algInvocId = xstr(algInvocation.getAlgorithmInvocationId())
    algInvocation.close()

    warning("Algorithm Invocation ID", algInvocId)

    vepParser = VepJsonParser(args.rankingFile, verbose=True)

    indexer = BinIndex(args.gusConfigFile, verbose=False)

    database = Database(args.gusConfigFile)
    database.connect()

    load_annotation()

    database.close()
    indexer.close()
    print(algInvocId, file=sys.stdout)
