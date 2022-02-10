#!/usr/bin/env python
# pylint: disable=invalid-name,not-an-iterable,unused-import,too-many-locals
"""
Loads AnnotatedVDB from JSON output from running VEP on dbSNP
USE ONLY IF CERTAIN OF NO UNKNOWN CONSEQUENCES
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

from GenomicsDBData.Util.utils import xstr, die, qw, execute_cmd, xstrN, warning
from GenomicsDBData.Util.postgres_dbi import Database
from GenomicsDBData.Util.exceptions import print_exception

from AnnotatedVDB.BinIndex.bin_index import BinIndex
from AnnotatedVDB.Util.algorithm_invocation import AlgorithmInvocation
from AnnotatedVDB.Util.vcf_parser import VcfEntryParser
from AnnotatedVDB.Util.variant_annotator import VariantAnnotator, VARIANT_COLUMNS
from AnnotatedVDB.Util.vep_parser import VepJsonParser, CONSEQUENCE_TYPES

import multiprocessing
from concurrent.futures import ProcessPoolExecutor

from AnnotatedVDB.Util.chromosomes import Human


def evaluate_minor_allele(vepParser, lfh):
    ''' some dbSNP variants are reported in VCF w/single mutations
    but are actually multiallelic (have additional minor allele)
    and reported as such on the dbSNP website.  To catch these we
    need to evaluate the 1000Genomes minor allele frequency '''

    refSnpId = vepParser.get('ref_snp_id')
    refAllele = vepParser.get('ref_allele')
    altAlleles = vepParser.get('alt_allele')
    frequencies = vepParser.get_frequencies(refSnpId)

    if frequencies is None:
        return None

    if 'minor_allele_freq' not in frequencies: # nothing to change
        return frequencies

    minorAlleleFreq = frequencies['minor_allele_freq']
    minorAllele = frequencies['minor_allele']

    match = False

    aFreq = get_allele_frequencies(minorAllele, frequencies)

    # minor allele already in frequency list; update
    if aFreq is not None:
        match = True
        if '1000Genomes' in aFreq:
            frequencies['values'][minorAllele]['1000Genomes'].append({'gmaf': minorAlleleFreq})
        else:
            frequencies['values'][minorAllele]['1000Genomes'] = [{'gmaf': minorAlleleFreq}]

    # the minor allele is one of the alt alleles in the VCF, but not in frequency list; add
    elif minorAllele in altAlleles:
        match = True
        frequencies['values'][minorAllele]['1000Genomes'] = [{'gmaf': minorAlleleFreq}]

    # then check to see if it is a normalized allele (INDEL/DIV/INS/etc)
    else:
        for alt in altAlleles:
            nRef, nAlt = VcfEntryParser(None).normalize_alleles(refAllele, alt, snvDivMinus=True)
            if minorAllele == nAlt:
                match = True
                aFreq = get_allele_frequencies(nAlt, frequencies)
                if aFreq is not None:
                    # this should never happen, b/c it should have been caught above, but...
                    if '1000Genomes' in aFreq:
                        frequencies['values'][minorAllele]['1000Genomes'].append({'gmaf': minorAlleleFreq})
                    else:
                        frequencies['values'][minorAllele]['1000Genomes'] = [{'gmaf': minorAlleleFreq}]
                else:
                    frequencies['values'][minorAllele]['1000Genomes'] = [{'gmaf': minorAlleleFreq}]
                break # if match found, no need to iterate over the rest

            if minorAllele == nRef:
                # don't add it in; can/should be inferred from other allele frequencies
                match = True
                warning("Matched minor allele", minorAllele,
                        "to normalized reference for variant -",
                        refSnpId, refAllele + '/' + '/'.join(altAlleles), "- IGNORED", file=lfh, flush=True)
                break

    if not match:
        warning("Unable to match minor allele", minorAllele,
                "to variant -", refSnpId, "- alleles:",
                refAllele + '/' + '/'.join(altAlleles),
                "Assuming VCF/dbSNP record mismatch.", file=lfh, flush=True)
        frequencies['values'][minorAllele]['1000Genomes'] = [{'gmaf': minorAlleleFreq}]

        # also add it to the alternative alleles, so it will be loaded as its own variant
        vepParser.set('alt_allele', altAlleles.append(minorAllele)) 

    return frequencies


def get_frequencies(vepParser, lfh):
    ''' extract frequencies'''

    frequencies = evaluate_minor_allele(vepParser, lfh)
    if frequencies is None:
        return None
    return frequencies['values']


def json2str(value):
    ''' convert JSON to str; accounting for None values '''
    if value is None:
        return 'NULL'
    else:
        return json.dumps(value)


def load_annotation(chromosome):
    ''' parse over a JSON file, extract position, frequencies,
    ids, and ADSP-ranked most severe consequence; bulk load using COPY '''

    lfn = xstr(chromosome) + '.log';
    lfh = open(lfn, 'w')

    algInvocation = AlgorithmInvocation('parallel_load_dbsnp_vep_result.py', json.dumps({'chromosome': chromosome}), args.commit)
    algInvocId = xstr(algInvocation.getAlgorithmInvocationId())
    algInvocation.close()
    warning("Algorithm Invocation ID", algInvocId, file=lfh, flush=True)

    vepParser = VepJsonParser(args.rankingFile, verbose=True)
    indexer = BinIndex(args.gusConfigFile, verbose=False)

    database = Database(args.gusConfigFile)
    database.connect()

    fname = 'chr' + xstr(chromosome) + '.json.gz'
    fname = path.join(args.dir, fname)
    lineCount = 0
    variantCount = 0
    skipCount = 0

    resume = True if args.resumeAfter is None else False
    if not resume:
        warning("--resumeAfter flag specified; Finding skip until point", args.resumeAfter, file=lfh, flush=True)

    previousSnp = None
    with database.cursor() as cursor:
        copyObj = io.StringIO()
        with open(fname, 'r') as fhandle:
            mappedFile = mmap.mmap(fhandle.fileno(), 0, prot=mmap.PROT_READ) # put file in swap
            with gzip.GzipFile(mode='r', fileobj=mappedFile) as gfh:
                for line in gfh:
                    lineCount = lineCount + 1
                    vepResult = json.loads(line.rstrip())

                    vepParser.set_annotation(copy.deepcopy(vepResult))
                    vepInputStr = vepParser.get('input')
                    entry = VcfEntryParser(vepInputStr)

                    # there are json formatting issues w/the input str
                    # so replace w/the parsed entry; which is now a dict
                    vepResult['input'] = entry.get_entry()

                    refSnpId = entry.get_refsnp()

                    if not resume:
                        if previousSnp == args.resumeAfter:
                            warning(previousSnp, refSnpId, file=lfh, flush=True)
                            warning("Resuming after:", args.resumeAfter, "- SKIPPED",
                                    skipCount, "lines.", file=lfh, flush=True)
                            resume = True
                        else:
                            previousSnp = refSnpId
                            skipCount = skipCount + 1
                            continue

                    if lineCount == 1 or variantCount % args.commitAfter == 0:
                        warning('Processing new copy object', file=lfh, flush=True)
                        tstart = datetime.now()

                    vepParser.set('ref_snp_id', refSnpId)

                    refAllele = entry.get('ref')
                    altAllele = entry.get('alt').split(',')

                    chrom = xstr(entry.get('chrom'))
                    if chrom == 'MT':
                        chrom = 'M'
                    position = int(entry.get('pos'))

                    isMultiAllelic = len(altAllele) > 1

                    try:

                        vepParser.set('is_multi_allelic', isMultiAllelic)
                        vepParser.set('ref_allele', refAllele)
                        vepParser.set('alt_allele', altAllele)

                        frequencies = get_frequencies(vepParser, lfh)
                        vepParser.adsp_rank_and_sort_consequences()

                        # for each allele
                        for alt in altAllele:
                            variantCount = variantCount + 1
                            
                            annotator = VariantAnnotator(refAllele, altAllele, chrom, position);
                            normRef, normAlt = annotator.get_normalized_alleles()
                            annotator.set_location_end(entry.infer_variant_end_location(alt, normRef));
                            
                            binIndex = indexer.find_bin_index(chrom, position, annotator.get_location_end())

                            # NOTE: VEP uses normalized alleles to indicate variant_allele

                            alleleFreq = None if frequencies is None \
                              else get_allele_frequencies(normAlt, frequencies)

                            msConseq = get_most_severe_consequence(normAlt, vepParser)

                            valueStr = '#'.join((
                                'chr' + xstr(chrom),
                                xstr(position),
                                xstr(vepParser.get('is_multi_allelic')),
                                binIndex,
                                refSnpId,
                                annotator.get_metaseq_id(),
                                json2str(annotator.get_display_attributes()),
                                json2str(alleleFreq),
                                json2str(msConseq),
                                json2str(get_adsp_ranked_allele_consequences(normAlt, vepParser)),
                                json.dumps(vepResult),
                                algInvocId
                                ))

                            copyObj.write(valueStr + '\n')

                            if variantCount % args.logAfter == 0 \
                              and variantCount % args.commitAfter != 0:
                                warning("PARSED", variantCount, file=lfh, flush=True)

                            if variantCount % args.commitAfter == 0:
                                tendw = datetime.now()
                                warning('Copy object prepared in ' + str(tendw - tstart) + '; ' +
                                        str(copyObj.tell()) + ' bytes; transfering to database', file=lfh, flush=True)
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

                                if variantCount % args.logAfter == 0:
                                    warning(message, "; up to = ", refSnpId, file=lfh, flush=True)

                                tend = datetime.now()
                                warning('Database copy time: ' + str(tend - tendw), file=lfh, flush=True)
                                warning('        Total time: ' + str(tend - tstart), file=lfh, flush=True)

                                if args.test:
                                    die("Test complete")

                                copyObj = io.StringIO() # reset io string

                    except Exception as e:
                        warning("ERROR parsing variant", refSnpId, file=lfh, flush=True)
                        warning(lineCount, ":", line, file=lfh, flush=True)
                        warning(str(e), file=lfh, flush=True)
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

            mappedFile.close()
        
    database.close()
    indexer.close()
    lfh.close()



def get_most_severe_consequence(allele, vepParser):
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


def get_adsp_ranked_allele_consequences(allele, vepParser):
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
    parser = argparse.ArgumentParser(description='load AnnotatedDB from  JSON output of VEP against dbSNP')
    parser.add_argument('-d', '--dir',
                        help="directory containing VEP results", required=True)
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
                        help="chromosome to load, e.g., 1, 2, M, X")
    parser.add_argument('--commitAfter', type=int, default=500,
                        help="commit after specified inserts")
    parser.add_argument('--maxWorkers', default=10, type=int)
    parser.add_argument('--logAfter', type=int,
                        help="number of inserts to log after completion; will work best if factor/multiple of commitAfter")
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()

    if not args.logAfter:
        args.logAfter = args.commitAfter

    chrList = args.chr.split(',') if args.chr != 'all' \
      else [c.value for c in Human]
 
    with ProcessPoolExecutor(args.maxWorkers) as executor:
        for c in chrList:
            warning("Create and start thread for chromosome:", xstr(c))
            executor.submit(load_annotation, c)
