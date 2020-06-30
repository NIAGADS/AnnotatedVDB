#!/usr/bin/env python
# pylint: disable=invalid-name,not-an-iterable,unused-import,too-many-locals
"""
Loads AnnotatedVDB from JSON output from running VEP on a non-dbSNP VCF
NOTE: assumes variants were not mapped to a refSNP
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
from AnnotatedVDB.Util.vep_parser import VepJsonParser, CONSEQUENCE_TYPES


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

    resume = True if args.resumeAfter is None else False

    previousSnp = None
    with database.cursor() as cursor:
        copyObj = io.StringIO()
        with open(fname, 'r') as fh:
            with open(args.logFile) as lfh:
                for line in fh:
                    lineCount = lineCount + 1
                    vepResult = json.loads(line.rstrip())

                    vepParser.set_annotation(copy.deepcopy(vepResult))
                    vepInputStr = vepParser.get('input')
                    entry = VcfEntryParser(vepInputStr)

                    # there may be json formatting issues w/the input str
                    # so replace w/the parsed entry; which is now a dict
                    vepResult['input'] = entry.get_entry()

                    if lineCount == 1 or variantCount % 500 == 0:
                        warning('Processing new copy object', file=lfh, flush=True)
                        tstart = datetime.now()

                    refAllele = entry.get('ref')
                    altAllele = entry.get('alt')  # assumes no multialleic variants

                    chrom = xstr(entry.get('chrom'))
                    if chrom == 'MT':
                        chrom = 'M'

                    position = int(entry.get('pos'))

                    metaseqId = ':'.join((chrom, xstr(position), refAllele, altAllele))

                    if not resume:
                        if previousSnp == args.resumeAfter:
                            warning(previousSnp, metaseqId, file=lfh, flush=True)
                            warning("Resuming after:", args.resumeAfter, "- SKIPPED",
                                    skipCount, "lines.", file=lfh, flush=True)
                            resume = True
                        else:
                            previousSnp = metaseqId
                            skipCount = skipCount + 1
                            continue

                    try:

                        vepParser.set('ref_allele', refAllele)
                        vepParser.set('alt_allele', altAllele)

                        frequencies = get_frequencies()
                        warning(frequencies, file=lfh, flush=True)
                        die("DEBUG - done")
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

                        if variantCount % 500 == 0:
                            warning("PARSED", variantCount, file=lfh, flush=True)
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

                            if variantCount % 500 == 0:
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
    parser.add_argument('-r', '--rankingFile', required=True,
                        help="full path to ADSP VEP consequence ranking file")
    parser.add_argument('--commit', action='store_true', help="run in commit mode", required=False)
    parser.add_argument('--gusConfigFile',
                        '--full path to gus config file, else assumes $GUS_HOME/config/gus.config')
    parser.add_argument('--resumeAfter',
                        help="refSnpId after which to resume load (log lists lasts committed refSNP)")
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()


    VARIANT_COLUMNS = qw('chromosome location bin_index metaseq_id allele_frequencies adsp_most_severe_consequence adsp_ranked_consequences vep_output row_algorithm_id', returnTuple=True)

    algInvocation = AlgorithmInvocation('load_vep_result.py', json.dumps(vars(args)), args.commit)
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
