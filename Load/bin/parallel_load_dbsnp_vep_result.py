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
import io
import copy
import csv
import traceback

from datetime import datetime
from os import path
import multiprocessing
from concurrent.futures import ProcessPoolExecutor

from GenomicsDBData.Util.utils import xstr, die, execute_cmd, xstrN, warning, print_dict, to_numeric
from GenomicsDBData.Util.postgres_dbi import Database, raise_pg_exception
from GenomicsDBData.Util.exceptions import print_exception

from AnnotatedVDB.BinIndex.bin_index import BinIndex
from AnnotatedVDB.Util.algorithm_invocation import AlgorithmInvocation
from AnnotatedVDB.Util.vcf_parser import VcfEntryParser
from AnnotatedVDB.Util.variant_annotator import VariantAnnotator, BASE_LOAD_FIELDS, reverse_complement
from AnnotatedVDB.Util.vep_parser import VepJsonParser, CONSEQUENCE_TYPES
from AnnotatedVDB.Util.primary_key_generator import VariantPKGenerator
from AnnotatedVDB.Util.chromosomes import Human

COPY_SQL = "COPY AnnotatedVDB.Variant("  \
  + ','.join(BASE_LOAD_FIELDS) \
  + ") FROM STDIN WITH (NULL 'NULL', DELIMITER '#')"


def insert_variants(copyObj, cursor):  
    """! perform copy operation to insert variants into the DB
        @param copyObj      io.StringIO stream containing values to be inserted
        @param cursor       database cursor
    """ 
    try:
        copyObj.seek(0)
        cursor.copy_expert(COPY_SQL, copyObj, 2**10)
    except Exception as e:
        raise_pg_exception(e)

  
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


def get_vcf_frequencies(allele, vcfEntry):
    """! retrieve allele frequencies reported in INFO FREQ field
    @param allele          allele to match (not normalized)
    @param vcfEntry        VCF entry (input to VEP)
    @param lfh             log file handle
    @returns               dict of source:frequency for the allele
    """

    zeroValues = ['.', '0']
    
    entryFrequencies = vcfEntry.get_info('FREQ')
    if entryFrequencies is None:
        return None

    # FREQ=GnomAD:0.9986,0.001353|Korea1K:0.9814,0.01861|dbGaP_PopFreq:0.9994,0.0005901
    # altIndex needs to be incremented as first value is for the ref allele)
    altIndex = vcfEntry.get('alt').split(',').index(allele) + 1
    populationFrequencies = {pop.split(':')[0]:pop.split(':')[1] for pop in entryFrequencies.split('|')}
    vcfFreqs = {pop: {'gmaf': to_numeric(freq.split(',')[altIndex])} \
                for pop, freq in populationFrequencies.items() \
                if freq.split(',')[altIndex] not in zeroValues}
  
    return None if len(vcfFreqs) == 0 else vcfFreqs


def update_frequencies(dict1, dict2):
    """! update a frequency dict by combining dict1 & dict2, making checks for NoneTypes
    @returns updated dict frequency dict
    """
    if dict1 is None and dict2 is None: # really not necessary, but there for clarity
        return None
    if dict1 is None:
        return dict2
    if dict2 is None:
        return dict1

    # TODO - update nested dicts
    
    dict1.update(dict2)
    return dict1


def get_frequencies(vepParser, lfh):
    """! extract frequencies"""

    refSnpId = vepParser.get('ref_snp_id')
    frequencies = vepParser.get_frequencies(refSnpId)
    
    if frequencies is None:
        return None
    return frequencies['values']


def json2str(value):
    """! convert JSON to str; accounting for None values """
    if value is None:
        return 'NULL'
    else:
        return json.dumps(value)


def load_annotation(chromosome):
    """! parse over a JSON file, extract position, frequencies,
    ids, and ADSP-ranked most severe consequence; bulk load using COPY """

    lfn = "annotatedvdb_dbsnp_chr_" + xstr(chromosome) + '.log';
    warning("Logging to", lfn)
    lfh = open(lfn, 'w')
    if args.verbose:
        warning("Parameters:", print_dict(vars(args), pretty=True), file=lfh, flush=True)

    algInvocation = AlgorithmInvocation('parallel_load_dbsnp_vep_result.py',
                                        json.dumps({'chromosome': chromosome}), args.commit)
    algInvocId = xstr(algInvocation.getAlgorithmInvocationId())
    algInvocation.close()
    warning("Algorithm Invocation ID", algInvocId, file=lfh, flush=True)

    pkGenerator = VariantPKGenerator(args.genomeBuild, args.seqrepoProxyPath)
    vepParser = VepJsonParser(args.rankingFile, rankConsequencesOnLoad=True, verbose=args.verbose)
    indexer = BinIndex(args.gusConfigFile, verbose=False)

    database = Database(args.gusConfigFile)
    database.connect()

    fname = 'chr' + xstr(chromosome) + "." + args.extension
    fname = path.join(args.dir, fname)
    lineCount = 0
    variantCount = 0
    skipCount = 0
    debugCount = 0

    resume = args.resumeAfter is None
    if not resume:
        warning("--resumeAfter flag specified; Finding skip until point", args.resumeAfter, file=lfh, flush=True)

    previousSnp = None
    with database.cursor() as cursor:
        copyObj = io.StringIO()
        with open(fname, 'r') as fhandle:
            mappedFile = mmap.mmap(fhandle.fileno(), 0, prot=mmap.PROT_READ) # put file in swap
            with gzip.GzipFile(mode='r', fileobj=mappedFile) as gfh:
                for line in gfh:
                    try: 
                        lineCount = lineCount + 1
                        vepResult = json.loads(line.rstrip())

                        vepParser.set_annotation(copy.deepcopy(vepResult))
                        vepInputStr = vepParser.get('input')
                        entry = VcfEntryParser(vepInputStr)

                        if chrmMap is not None:
                            entry.update_chromosome(chrmMap)

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

                        if lineCount == 1 or lineCount % args.commitAfter == 0:
                            if args.debug:
                                warning('Processing new copy object', file=lfh, flush=True)
                            tstart = datetime.now()

                        vepParser.set('ref_snp_id', refSnpId)

                        refAllele = entry.get('ref')
                        altAllele = entry.get('alt').split(',')
    
                        chrom = xstr(entry.get('chrom'))
                        if chrom == 'MT':
                            chrom = 'M'

                        position = int(entry.get('pos'))
                        rsPosition = entry.get_info('RSPOS')
                
                        isMultiAllelic = len(altAllele) > 1

                        vepParser.set('is_multi_allelic', isMultiAllelic)
                        vepParser.set('ref_allele', refAllele)
                        vepParser.set('alt_allele', altAllele)

                        frequencies = get_frequencies(vepParser, lfh)

                        vepParser.adsp_rank_and_sort_consequences()

                    # for each allele
                        for alt in altAllele:
                            variantCount = variantCount + 1

                            annotator = VariantAnnotator(refAllele, alt, chrom, position)
                            normRef, normAlt = annotator.get_normalized_alleles()
                            binIndex = indexer.find_bin_index(chrom, position,
                                                                entry.infer_variant_end_location(alt, normRef))

                            # NOTE: VEP uses left normalized alleles to indicate freq allele
                            # so need to use normalized alleles to match freq allele
                            # to the correct variant alt allele
                            
                            alleleFreq = None if frequencies is None \
                                else get_allele_frequencies(normAlt, frequencies)
        
                            alleleFreq = update_frequencies(alleleFreq, get_vcf_frequencies(alt, entry))
                        
                            msConseq = get_most_severe_consequence(normAlt, vepParser)

                            recordPK = pkGenerator.generate_primary_key(annotator.get_metaseq_id(), refSnpId)
                            
                            # fields: ('chromosome record_primary_key position is_multi_allelic bin_index ref_snp_id metaseq_id display_attributes allele_frequencies adsp_most_severe_consequence adsp_ranked_consequences vep_output row_algorithm_id')
                            valueStr = '#'.join((
                                'chr' + xstr(chrom),
                                recordPK,
                                xstr(position),
                                xstr(vepParser.get('is_multi_allelic')),
                                binIndex,
                                refSnpId,
                                annotator.get_metaseq_id(),
                                json2str(annotator.get_display_attributes(rsPosition)),
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
                                if args.debug:
                                    tendw = datetime.now()
                                    warning('Copy object prepared in ' + str(tendw - tstart) + '; ' +
                                            str(copyObj.tell()) + ' bytes; transfering to database',
                                                file=lfh, flush=True)

                                insert_variants(copyObj, cursor)

                                message = '{:,}'.format(variantCount)
                                if args.commit:
                                    database.commit()
                                    message = "COMMITTED " + message

                                else:
                                    database.rollback()
                                    message = "PARSED " + message + " -- rolling back"

                                if variantCount % args.logAfter == 0:
                                    warning(message, "; up to = ", refSnpId, file=lfh, flush=True)

                                if args.debug:
                                    tend = datetime.now()
                                    warning('Database copy time: ' + str(tend - tendw), file=lfh, flush=True)
                                    warning('        Total time: ' + str(tend - tstart), file=lfh, flush=True)

                                copyObj.close() # free the memory
                                copyObj = io.StringIO() # reset io string

                                if args.test:
                                    warning("DONE - TEST COMPLETE", file=lfh, flush=True)
                                    sys.exit("EXITING - TEST COMPLETE")

                    except Exception as e:
                        warning("ERROR parsing variant", refSnpId, file=lfh, flush=True)
                        warning(lineCount, ":", print_dict(json.loads(line)), file=lfh, flush=True)
                        warning(str(e), file=lfh, flush=True)
                        warning(traceback.format_exc(), file=lfh, flush=True)
                        raise

            # final commit / leftovers
            insert_variants(copyObj, cursor)

            message = '{:,}'.format(variantCount)

            if args.commit:
                database.commit()
                message = "DONE - COMMITTED " + message
            else:
                database.rollback()
                message = "DONE - PARSED " + message + " -- rolling back"

            warning(message, file=lfh, flush=True)
            warning("NEW CONSEQUENCES -", vepParser.get_added_conseq_summary())

            mappedFile.close()
            
    database.close()
    indexer.close()
    lfh.close()


def get_most_severe_consequence(allele, vepParser):
    """ retrieve most severe consequence from the VEP JSON Parser,
    for the specified allele; returns None if no consequences are found
    return first hit among transcript, then regulatory feature, then intergenic
    consequences
    """

    for ct in CONSEQUENCE_TYPES:
        msConseq = get_allele_consequences(allele, vepParser.get(ct + '_consequences'))
        if msConseq is not None:
            return msConseq[0]

    return None


def get_adsp_ranked_allele_consequences(allele, vepParser):
    """ get dict of all ADSP ranked allele consequences """
    adspConseq = {}
    for ctype in CONSEQUENCE_TYPES:
        ctypeKey = ctype + '_consequences'
        cq = get_allele_consequences(allele, vepParser.get(ctypeKey))
        if cq is not None:
            adspConseq[ctypeKey] = cq

    return None if len(adspConseq) == 0 else adspConseq


def get_allele_consequences(allele, conseq):
    """ check to see if conseq dict contains values for the specified allele,
    if so return first value associated w/the allele as the most
    serious; assumes ADSP ranking and re-sort has occurred """

    if conseq is None:
        return None

    if allele in conseq:
        return conseq[allele]

    return None


def get_allele_frequencies(allele, frequencies):
    """ given an allele and frequency dict,
    retrieve and return allele-specific frequencies """

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
