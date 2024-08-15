#!/usr/bin/env python3
# pylint: disable=invalid-name,not-an-iterable,unused-import,too-many-locals
"""
Loads/Updates AnnotatedVDB from ADSP QC pVCF.  Flags novel variants for CADD update
Loads each chromosome in parallel
"""

from __future__ import print_function

import argparse
import glob

from copy import deepcopy
from datetime import datetime
from os import path
from sys import stdout

from concurrent.futures import ProcessPoolExecutor
from psycopg2 import DatabaseError

from GenomicsDBData.Util.utils import xstr, warning, print_dict, print_args, die, get_opener, chunker, disallowed2str
from GenomicsDBData.Util.postgres_dbi import Database, raise_pg_exception

from AnnotatedVDB.Util.loaders import VCFVariantLoader
from AnnotatedVDB.Util.parsers import ChromosomeMap, VcfEntryParser
from AnnotatedVDB.Util.enums import HumanChromosome as Human

NUM_BULK_LOOKUPS = 1000


def load(loader, database, lookups, response):
    variantCount = 0
    for variant in lookups:
        variantCount += 1
        if variant in response: 
            if args.debug:
                loader.log(("Lookup result for", variant + ': ', response[variant]), prefix="DEBUG")
                
            if response[variant] is not None:
                # variants may have multiple matches if a metaseq 
                # matches multiple refsnps, so an array is returned by the lookup                
                for hit in response[variant]:   
                    
                    updateFlags = {'record_primary_key': hit['record_primary_key'], 
                                    'is_adsp_variant': hit['is_adsp_variant'],
                                    'adsp_qc': hit['annotation']['ADSP_QC'] is not None and args.version.lower() in hit['annotation']['ADSP_QC']}
                    if args.debug:
                        loader.log(("Update flags:", updateFlags), prefix="DEBUG")
                        
                    loader.parse_variant(lookups[variant], updateFlags) 
            else: 
                loader.parse_variant(lookups[variant])
                    
        else: # TODO shouldn't happen / probably should raise an error
            loader.parse_variant(lookups[variant])

        if variantCount % args.commitAfter == 0:
            if loader.copy_buffer(sizeOnly=True) > 0:
                loader.load_variants()
            if loader.update_buffer(sizeOnly=True) > 0 :
                loader.update_variants()
            log_load(loader, database)
                
    # residuals
    if loader.copy_buffer(sizeOnly=True) > 0:
            loader.load_variants()
    if loader.update_buffer(sizeOnly=True) > 0 :
        loader.update_variants()
    
    

def log_load(loader, database):
    message = "INSERTED " + '{:,}'.format(loader.get_count('variant')) + " variants"
    messagePrefix = "COMMITTED"
    if args.commit:
        database.commit()
    else:
        database.rollback()
        messagePrefix = "ROLLING BACK"

    message += "; up to = " + loader.get_current_variant_id()
                    
    if loader.get_count('update') > 0:
        message += "; UPDATED " + '{:,}'.format(loader.get_count('update')) + " variants"
    
    if loader.get_count('skipped') > 0:
        message += "; SKIPPED " + '{:,}'.format(loader.get_count('skipped')) + " variants"    
        
    loader.log(message, prefix=messagePrefix)
                                


def bulk_lookup(loader, lookups):
    """ do the bulk lookup in batches of 200 or less to be faster """
    numLookups = len(lookups.keys())
    chunks = chunker(list(lookups.keys()), NUM_BULK_LOOKUPS) if numLookups > NUM_BULK_LOOKUPS else None
    if chunks is None:
        return loader.variant_validator().bulk_lookup(list(lookups.keys()))
    else:
        finalResponse = {}
        chunkNum = 0
        for group in chunks:
            if args.debug:
                chunkNum += 1
                loader.log(("Check lookups - performing chunk", xstr(chunkNum), "; size =", xstr(len(group))), prefix="DEBUG")
            response = loader.variant_validator().bulk_lookup(group, firstHitOnly=False)
            finalResponse.update(response)
            
        if args.debug:
            loader.log(("Check lookups:", len(lookups.keys())==len(finalResponse.keys())), prefix="DEBUG")
        return finalResponse
            

def generate_update_values(loader, entry, flags):
    info = entry.get('info')
    filter = entry.get('filter')
    qual = entry.get('qual')
    format = entry.get('format', raiseError=False) # set value to None if empty
    release = loader.get_datasource() # ADSP release
    
    if loader.debug():
        loader.log(("Generate Update Values:", flags), prefix="DEBUG")

    recordPK = flags['record_primary_key'] if flags is not None else None
    isAdspVariant = flags['is_adsp_variant'] if flags is not None else False
    hasAdspQC = flags['adsp_qc'] if flags is not None else False
    adspFlag = True if filter == 'PASS' else 'NULL'
            
    qcValues = {
        release : {
            "info": info,
            "filter": filter,
            "qual": qual,
            "format": format
        }
    }
    
    qcStr = xstr(qcValues)
    if 'Infinity' in qcStr:
        loader.log(("Infinity Found - String Version: " +  qcStr), prefix="ERROR")
        loader.log(print_dict(qcValues, pretty=True), prefix="ERROR")
        raise ValueError("Infinity found among QC scores")
        
    # returns recordPK, flags, update values dict
    return [recordPK, {'is_adsp_variant': isAdspVariant, 'update': args.updateExistingValues or not hasAdspQC}, 
            {'is_adsp_variant': adspFlag, 'adsp_qc': qcValues}]



def initialize_loader(logFilePrefix, chrom=None):
    """! initialize loader """

    lfn = xstr(logFilePrefix)
    if args.resumeAfter:
        lfn += '_resume_' + args.resumeAfter.replace(':','-') 
    if args.failAt:
        lfn += '_fail_' + args.failAt.replace(':', '-')
    lfn += '.log'
    warning("Logging to", lfn)
    
    try:
        loader = VCFVariantLoader(args.version, logFileName=lfn, verbose=args.verbose, debug=args.debug) # don't want to update is_adsp_flags based on datasource but PASS status

        loader.log(("Parameters:", print_dict(vars(args), pretty=True)), prefix="INFO")

        loader.set_algorithm_invocation('load_vcf_result', xstr(logFilePrefix) + '|' + print_args(args, False))
        loader.log('Algorithm Invocation Id = ' + xstr(loader.alg_invocation_id()), prefix="INFO")
        
        loader.initialize_pk_generator(args.genomeBuild, args.seqrepoProxyPath)
        loader.initialize_bin_indexer(args.gusConfigFile)
        
        if chrom is not None:
            loader.set_chromosome(chrom)

        updateFields = ["is_adsp_variant", "adsp_qc"]
        loader.set_update_fields(updateFields)
        loader.build_update_sql()
        loader.initialize_copy_sql(copyFields=updateFields) # use default copy fields
        
        loader.set_chromosome_map(chrmMap)
        
        loader.initialize_variant_validator(args.gusConfigFile)
        
        loader.set_update_value_generator(generate_update_values)
        loader.set_update_existing(True)
        
        if args.failAt:
            loader.set_fail_at_variant(args.failAt)
            loader.log("Fail at variant: " + loader.fail_at_variant(), prefix="INFO")
        
    except Exception as err:
        warning("ERROR", "Problem initializing Loader")
        raise(err)
        
    return loader


def load_annotation(fileName, logFilePrefix, chrom=None):
    """! parse over a VCF file; bulk load using COPY """
 
    loader = initialize_loader(logFilePrefix, chrom)
    loader.log('Parsing ' + fileName, prefix="INFO")
    
    resume = args.resumeAfter is None # false if need to skip lines
    if not resume:
        loader.log(("--resumeAfter flag specified; Finding skip until point", args.resumeAfter), prefix="INFO")
        # loader.set_resume_after_variant(args.resumeAfter) --> handle resume during lookup before load

    previousLine = None
    vcfHeaderFields = args.vcfHeaderFields.split(',') if args.vcfHeaderFields else None
    lookups = {}
    try: 
        database = Database(args.gusConfigFile)
        database.connect()
        opener = get_opener(args.fileName, compressed=True, binary=True)
        if args.debug:
            loader.log('Preparing to do updates/inserts', prefix="DEBUG")
        with opener(fileName, 'r') as fhandle, database.cursor() as cursor:
            loader.set_cursor(cursor)
            # mappedFile = mmap.mmap(fhandle.fileno(), 0, prot=mmap.PROT_READ) # put file in swap
            lineCount = 0
            variantCount = 0
            for line in fhandle: # iter(mappedFile.readline, b""):
                line = line.decode("utf-8")  # in python 3, mapped file is a  binary IO object
                line = line.rstrip()
                if line.startswith("#"): # skip comments
                    previousLine = line
                    continue                    
                
                if previousLine is not None and previousLine.startswith("#") and vcfHeaderFields is None: 
                    # header, b/c this line does not start with "#"
                    vcfHeaderFields = previousLine.split('\t')
                    previousLine = None
                
                lineCount += 1 
                if args.debug and lineCount % 10000 == 0:
                    loader.log(("Read", lineCount, "lines"), prefix="DEBUG")
                    
                entry = VcfEntryParser(line, vcfHeaderFields, debug=True, loader=loader)
                currentVariant = entry.get_variant()
                for alt in currentVariant['alt_alleles']:                  
                    if not resume:
                        if currentVariant['id'] == args.resumeAfter:
                            loader.log("Resume point found", prefix="INFO")
                            resume = True
                        if lineCount % 50000 == 0:
                            loader.log(("SKIPPED", lineCount, "lines"), prefix="INFO")
                        break
                    
                    variantCount += 1
                    id = ':'.join((xstr(currentVariant['chromosome']),
                        xstr(currentVariant['position']), 
                        currentVariant['ref_allele'],
                        alt))
                    lookups[id] = entry
                
                if lineCount % args.numLookups == 0 and resume:                
                    response = bulk_lookup(loader, lookups)
                    load(loader, database, lookups, response)  
                    del lookups
                    del response
                    lookups = {}
                                        
                    log_load(loader, database)

                if args.test:
                    break

            
            # ============== end mapped file ===================
            
            # commit anything left over
            if len(lookups.keys()) > 0:
                loader.log("Processing residuals", prefix="INFO")
                response = bulk_lookup(loader, lookups)
                load(loader, database, lookups, response)   
                log_load(loader, database)
                
            if args.test:
                loader.log("DONE - TEST COMPLETE" , prefix="WARNING")
                
        # ============== end with open, cursor ===================
        
    except DatabaseError as err:
        problematicVariant = loader.get_current_variant_id()
        if problematicVariant is None:
            problematicVariant = currentVariant['id']
        loader.log("Problem submitting insert or update, error involves any variant from last COMMIT until "  \
            + problematicVariant, prefix="ERROR")
        raise(err)
    except IOError as err:
        raise(err)
    except Exception as err:
        problematicVariant = loader.get_current_variant_id()
        if problematicVariant is None:
            problematicVariant = currentVariant['id']
        loader.log("Problem parsing variant: " + problematicVariant, prefix="ERROR")
        loader.log((lineCount, ":", line), prefix="LINE")
        loader.log(str(err), prefix="ERROR")
        raise(err)
    finally:
        # mappedFile.close()
        database.close()
        loader.close()
        print(loader.get_algorithm_invocation_id(), file=stdout)


def validate_args():
    """! validate the parameters, print warnings, and update some values as necessary """
    if args.failAt:
        warning("--failAt option provided / running in NON-COMMIT mode")
        args.commit = False

    if args.fileName and args.chr:
        die("Both --fileName and --chr supplied, please specify only one, see --usage")

    if not args.fileName:
        if not args.chr:
            warning("Neither --fileName of --chr arguments supplied, assuming parallel load of all chromosomes")
            args.chr = 'all'
        else:        
            if not args.dir:
                die("Loading by chromosome, please supply path to directory containing the VEP results")
            if not args.extension:
                die("Loading by chromosome, please provide file extension. Assume file names are like chr1.extension / TODO: pattern matching")
   
         
def get_input_file_name(chrm):
    """ find the file that matches the chromosome & specified extension """  
    pattern = path.join(args.dir, '*chr' + xstr(chrm) + args.extension)
    fileName = glob.glob(pattern)   # *chr b/c there may be a prefix
    return path.join(args.dir, fileName[0])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(allow_abbrev=False, # otherwise it can substitute --chr for --chromosomeMap
                                    description='load AnnotatedDB from a VCF file, specify either a file or one or more chromosomes')
    parser.add_argument('-d', '--dir',
                        help="directory containing VCF files / only necessary for parallel load")
    parser.add_argument('-e', '--extension', 
                        help="file extension (e.g., vcf.gz) / required for parallel load")
    parser.add_argument('-g', '--genomeBuild', default='GRCh38', help="genome build: GRCh37 or GRCh38")
    parser.add_argument('-s', '--seqrepoProxyPath', required=True,
                        help="full path to local SeqRepo file repository")
    parser.add_argument('-m', '--chromosomeMap', required=False,
                        help="chromosome map")
    parser.add_argument('--commit', action='store_true', help="run in commit mode", required=False)
    parser.add_argument('--gusConfigFile',
                        '--full path to gus config file, else assumes $GUS_HOME/config/gus.config')
    parser.add_argument('--test', action='store_true',
                        help="load 'numLookups' rows as test")
    parser.add_argument('--resumeAfter',
                        help="variantId after which to resume load (log lists lasts committed variantId)")
    parser.add_argument('-c', '--chr', 
                        help="comma separated list of one or more chromosomes to load, e.g., 1, 2, M, X, `all`, `allNoM`, `autosome` / required for parallel load"),
    parser.add_argument('--fileName',
                        help="full path of file to load, if --dir option is provided, will be ignored")
    parser.add_argument('--commitAfter', type=int, default=500,
                        help="commit after specified inserts")
    parser.add_argument('--maxWorkers', default=10, type=int)
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--debug',  action='store_true',
                        help="log database copy time / may print development debug statements")
    parser.add_argument('--failAt', 
                        help="fail on specific variant and log output; if COMMIT = True, COMMIT will be set to False")
    parser.add_argument('--version', required=True, help='release version, e.g., 17K')
    parser.add_argument('--updateExistingValues', action="store_true",
                        help="update fields even if value exists; if flag not specificed will skip records that already have values")
    parser.add_argument('--numLookups', default=50000, type=int,
                        help="number of bulk lookups to do at a time")
    parser.add_argument('--vcfHeaderFields', help="comma separated list of fields to include", default='#CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT')
    args = parser.parse_args()


    validate_args()

    chrmMap = ChromosomeMap(args.chromosomeMap) if args.chromosomeMap else None

    if args.fileName:
        load_annotation(args.fileName, args.fileName + "-vcf-variant-loader")
        
    else:
        chrList = args.chr.split(',') if not args.chr.startswith('all') \
            else [c.value for c in Human]

        numChrs = len(chrList)
        with ProcessPoolExecutor(args.maxWorkers) as executor:
            if numChrs > 1:
                for c in chrList:
                    if args.chr == 'allNoM' and c == 'M':
                        continue 
                    if args.chr == 'autosome' and c in ['X', 'Y', 'M', 'MT']:
                        continue
                    warning("Create and start thread for chromosome:", xstr(c))
                    inputFile = get_input_file_name(c)
                    executor.submit(load_annotation, fileName=inputFile, logFilePrefix=inputFile + "-vcf-variant-loader", chrom=xstr(c))
            else: # debugs better w/out thread overhead, so single file -- no threading
                inputFile = get_input_file_name(chrList[0])
                load_annotation(inputFile, logFilePrefix=inputFile + "-vcf-variant-loader", chrom=xstr(chrList[0]))
