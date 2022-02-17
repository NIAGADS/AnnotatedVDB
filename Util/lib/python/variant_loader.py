"""! @brief Annotated VDB Loader """
#!pylint: disable=invalid-name
##
# @file variant_loader.py
#
# @brief  Annotated VDB Loader
#
# @section variant_loader Description
# provides functions & class for managing variant loads
#
# @section todo_variant_loader TODO
# - create __verify_current_variant() to make sure __current_variant is set before trying to access
# - extract common elements to parent class and make VEPLoader a child
#
# @section libraries_variant_loader Libraries/Modules
#
# - TBD
#
# @section author_variant_loader Author(s)
# - Created by Emily Greenfest-Allen (fossilfriend) 2022

import json 
import copy
import traceback

from io import StringIO

from GenomicsDBData.Util.utils import xstr, warning, print_dict, print_args, to_numeric, deep_update
from GenomicsDBData.Util.list_utils import qw, is_subset, is_equivalent_list
from GenomicsDBData.Util.postgres_dbi import Database, raise_pg_exception
from AnnotatedVDB.BinIndex.bin_index import BinIndex
from AnnotatedVDB.Util.algorithm_invocation import AlgorithmInvocation
from AnnotatedVDB.Util.vcf_parser import VcfEntryParser
from AnnotatedVDB.Util.variant_annotator import VariantAnnotator, BASE_LOAD_FIELDS
from AnnotatedVDB.Util.vep_parser import VepJsonParser
from AnnotatedVDB.Util.primary_key_generator import VariantPKGenerator
from Load.bin.load_dbsnp_vep_result import get_allele_consequences

ALLOWABLE_COPY_FIELDS = ["chromosome", "record_primary_key", "position", 
                         "is_multi_allelic", "is_adsp_variant", "ref_snp_id", 
                         "metaseq_id", "bin_index", "display_attributes", 
                         "allele_frequencies", "cadd_scores", 
                         "adsp_most_severe_consequence", "adsp_ranked_consequences", 
                         "loss_of_function", "vep_output", "adsp_qc", 
                         "gwas_flags", "other_annotation", "row_algorithm_id"]

REQUIRED_COPY_FIELDS = ["chromosome", "record_primary_key", "position", "metaseq_id", "bin_index", "row_algorithm_id"]

DEFAULT_COPY_FIELDS = qw('chromosome record_primary_key position is_multi_allelic bin_index ref_snp_id metaseq_id display_attributes allele_frequencies adsp_most_severe_consequence adsp_ranked_consequences vep_output row_algorithm_id', returnTuple=True)

class VEPLoader(object):
    """! functions for loading variants from VEP result """
    
    def __init__(self, datasource, logFileName=None, verbose=False, debug=False):
        """! VariantLoader base class initializer

            @param copySql             SQL for copy statement
            @param chromosomeMap       chromosome map
            @param datasource          datasource description
            @param logFileName         full path to logging file
            @param verbose             flag for verbose output
            @param debug               flag for debug output
            
            @returns                   An instance of the VariantLoader class with initialized counters and copy buffer
        """
        self.__log_file_handle = None if logFileName is None else open(logFileName, 'w') 
        self.__verbose = verbose
        self.__debug = debug
        self.__commit = self.__args.commit
        self.__alg_invocation_id = None
        self.__pk_generator = None
        self.__vep_parser = None
        self.__bin_indexer = None
        self.__cursor = None # database curspr
        self.__counters = {}
        self.__chromosome_map = None
        self.__resume_after_variant = None
        self.__resume = True
        self.__datasource = datasource.lower() 
        self.__current_variant = {}
        self.__copy_buffer = None
        self.__copy_fields = None
        self.__copy_sql = None
 
        self.__initialize_counters()
        self.initialize_copy_buffer()

    def close(self):
        self.close_copy_buffer()
        self.close_log()
        
        
    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
    

    def set_chromosome_map(self, chrmMap):
        self.__chromosome_map = chrmMap

    def get_current_variant_id(self):
        return self.__current_variant.id


    def __verify_copy_fields(self, copyFields):
        if copyFields is None: # default of DEFAULT_COPY_FIELDS will be used
            return True
        
        if is_subset(copyFields, ALLOWABLE_COPY_FIELDS):
            if is_subset(REQUIRED_COPY_FIELDS, copyFields):
                return True
            else:
                err = ValueError("Copy fields must include the following required fields: " + xstr(REQUIRED_COPY_FIELDS))
                self.log(str(err), prefix="ERROR")
                raise err
        else:
            err = ValueError("Copy fields include invalid columns from AnnotatedVDB.Variant: se variant_loader.py for list of allowable fields")
            self.log(str(err), prefix="ERROR")
            raise err

    def initialize_copy_sql(self, copyFields=None):
        self.__verify_copy_fields(copyFields)   
        self.__copy_fields = DEFAULT_COPY_FIELDS if copyFields is None else copyFields
        
        if not is_equivalent_list(self.__copyFields, DEFAULT_COPY_FIELDS):
            # TODO - make VEPLoader a child of a more generic loader to handle user
            # specified copy fields
            err = NotImplementedError('Currently the VEP Variant Loader can only handle the default copy fields')
            self.log(str(err), prefix="ERROR")
            raise err
     
        self.__copy_sql = "COPY AnnotatedVDB.Variant("  \
         + ','.join(self._copy_fields) \
         + ") FROM STDIN WITH (NULL 'NULL', DELIMITER '#')"


    def close_copy_buffer(self):
        """! close the copy buffer """
        self.__copy_buffer.close()
        
        
    def reset_copy_buffer(self):
        """! reset copy buffer to clean up memory """    
        self.close_copy_buffer()
        self.initialize_copy_buffer()
    
    
    def initialize_copy_buffer(self):
        """! initialize the copy buffer (io.StringIO) """
        self.__copy_buffer = StringIO()
        
        
    def copy_buffer(self, sizeOnly=False):
        """! returns copy buffer
        @param sizeOnly               return size only 
        @returns copy buffer or buffer size """
        
        return self.__copy_buffer.tell() if sizeOnly else self.__copy_buffer
    
    def add_copy_str(self, copyStr):
        """! write a copy string to the copy buffer """
        self.__copy_buffer.write(copyStr + '\n')


    def is_dbsnp(self):
        """!  @returns  flag indicating if datasource is dbsnp
        """
        return self.__datasource == 'dbsnp'
    

    def resume_load(self):
        """! check if resume load flag has been toggled
        @returns    boolean stored in __resume
        """
        return self.__resume
    

    def set_resume_after_variant(self, variantId):
        """! resume after specified variant
        @param variantId           id (metaseq or refsnp) of the variant
        """
        self.__resume_after_variant = variantId
        self.__resume = False # skipping until variant is found
        

    def initialize_vep_parser(self, rankingFile, rankConsequencesOnLoad, verbose):
        """! initialize VEP Parser

            @param rankingFile                file containing consequence ranks
            @param rankConsequencesOnLoad     rank consequences on load?
            @param verbose                  
        """
        self.__vep_parser = VepJsonParser(rankingFile, rankConsequencesOnLoad=rankConsequencesOnLoad, verbose=verbose)
        
    
    def initialize_bin_indexer(self, gusConfigFile):
        """! initialize Bin Indexer 

            @param gusConfigFile             gus.config file w/database connection info
        """
        self.__bin_indexer = BinIndex(gusConfigFile, verbose=False)


    def initialize_pk_generator(self, genomeBuild, seqrepoProxyPath):
        """! initialize primary key generator

            @param genomeBuild             genome build (GRCh38 or GRCh37)
            @param seqrepoProxyPath        full path to the seqrepo file proxy
        """
        self.__pk_generator = VariantPKGenerator(genomeBuild, seqrepoProxyPath)
 
    
    def get_count(self, counter):
        """! get value of a counter
            @param counter           counter whose value should be returned
            @returns                 value of counter
        """
        return self.__counters[counter]
    
    
    def increment_counter(self, counter):
        """! increment counter by 1 """
        self.__counters[counter] = self.__counters[counter] + 1   
        
    
    def __initialize_counters(self):
        """! initialize counters"""
        self.__counters = { 'line' : 0, 'total_variants': 0, 'skipped_variants':0}
            
                
    def cursor(self):
        """! @returns database cursor """    
        return self.__cursor
    
    
    def set_cursor(self, cursor):
        """! set database handle
        @params  cursor              database cursor
        """
        self.__cursor = cursor
            
        
    def current_variant(self):
        """! @returns current variant """
        return self.__current_variant
        
        
    def vep_parser(self):
        """! @returns VEP Parser """
        return self.__vep_parser
    
    
    def bin_indexer(self):
        """! @returns Bin indexer """
        return self.__bin_indexer
        
        
    def alg_invocation_id(self):
        """! @returns algorithm invocation id """
        return self.__alg_invocation_id
    
    
    def pk_generators(self):
        """! @returns primary key generatory"""
        return self.__pk_generator
        
    
    def set_algorithm_invocation(self, callingScript, comment):
        """! create entry in AnnotatedVDB.AlgorithmInvocation table for the data load
          @returns                 algorithm invocation id
        """
        algInvocation = AlgorithmInvocation(callingScript, comment, self.__commit)
        algInvocId = xstr(algInvocation.getAlgorithmInvocationId())
        algInvocation.close()
        return algInvocId
    
        
    def log(self, message, prefix=None):
        """! print to log
            @param message             message to print to log -- string or tuple/list
            @param prefix              e.g., DEBUG, WARNING, ERROR
        """
        if not isinstance(message, str):
            message = ' '.join([xstr(x) for x in message])
        
        if prefix:
            message = prefix + ": " + message          
        warning(message, file=self.__log_file_handle, flush=True)
        
        
    def log_fh(self):
        """! @returns log file handle """
        return self.__log_file_handle
    
    
    def close_log(self):
        """! closes log file handle """
        self.__log_file_handle.close()
    
    
    def __update_resume_status(self, variantId):
        """! check to see if the row should be skipped 
            update __resume flag accordingly
            
            @param variantId            id field in vep output  
        """
        if not self.resume_load():
            self.increment_counter('skip')
            self.__resume = variantId == self.__resume_after_variant
            
            if self.resume_load() is True: # resume load when next line is parsed
                message = ("Resuming after", self.__resume_after_variant)
                self.log(message, prefix="WARNING")
                message = ("Skipped", xstr(self.get_count('skip'), "variants"))
                self.log(message, prefix="INFO")
    
    
    def __get_result_frequencies(self):
        """! get allele frequencies from the VEP result, catching null objects
            @returns                frequency JSON
        """
        matchingVariantId = self.__current_variant.ref_snp_id if self.is_dbnp() else None
        frequencies = self.__vep_parser.get_frequencies(matchingVariantId)
        if frequencies is None:
            return None
        return frequencies['values']     
    
    
    def __get_allele_frequencies(allele, frequencies):
        """! given an allele and frequency dict, retrieve and return allele-specific frequencies
            @param allele          allele to be matched
            @param frequencies     result frequencies as JSON/dict
            
            @returns 
        """
        if frequencies is None:
            return None

        if allele in frequencies:
            return frequencies[allele]

        return None
        
    
    def __parse_vcf_frequencies(self, allele, vcfGMAFs):
        """! retrieve allele frequencies reported in INFO FREQ field
        @param allele          allele to match (not normalized)
        @param vcfGMAF         global minor allele frequencies from the VCF FREQ info field
        @returns               dict of source:frequency for the allele
        """
        
        if vcfGMAFs is None:
            return None
        
        zeroValues = ['.', '0']

        # FREQ=GnomAD:0.9986,0.001353|Korea1K:0.9814,0.01861|dbGaP_PopFreq:0.9994,0.0005901
        # altIndex needs to be incremented as first value is for the ref allele)
        altIndex = self.__current_variant.alt_alleles.split(',').index(allele) + 1
        populationFrequencies = {pop.split(':')[0]:pop.split(':')[1] for pop in vcfGMAFs.split('|')}
        vcfFreqs = {pop: {'gmaf': to_numeric(freq.split(',')[altIndex])} \
                    for pop, freq in populationFrequencies.items() \
                    if freq.split(',')[altIndex] not in zeroValues}
    
        return None if len(vcfFreqs) == 0 else vcfFreqs


    def __add_vcf_frequencies(self, alleleFreqs, vcfFreqs, allele):
        """! add VCF GMAF from the VCF FREQ INFO field to allele frequencies 

            @param alleleFreqs         allele frequencies from VEP result
            @param vcfFreqs            FREQ string from VCF INFO field
            @param allele             _description_
            @returns                 _description_
        """
        
        alleleGMAFs = self.__parse_vcf_frequencies(allele, vcfFreqs)
        if alleleFreqs is None:
            return alleleGMAFs
        if alleleGMAFs is None:
            return alleleFreqs
        
        # update the alleleFreq dict, taking into account nestedness
        return deep_update(alleleFreqs, alleleGMAFs)
    
        
    def __parse_alt_alleles(self, vcfEntry, resultJson):
        """! iterate over alleles and calculate values to load
             add to copy buffer
            @param vcfEntry             the VCF entry for the current variant
            @param resultJSON           the VEP result / to be added to copyStr
        """
        # extract result frequencies
        frequencies = self.__get_result_frequencies()
        variant = self.__current_variant # just for ease/simplicity
        variantExternalId = variant.ref_snp_id if 'ref_snp_id' in variant else None
        
        for alt in variant.alt_alleles:
            self.increment_counter('variant')
            annotator = VariantAnnotator(variant.ref_allele, alt, variant.chromosome, variant.position)       
            normRef, normAlt = annotator.get_normalized_alleles()
            binIndex = self.__bin_indexer.find_bin_index(variant.chromosome, variant.position,
                                                         vcfEntry.infer_variant_end_location(alt, normRef))
            
            # NOTE: VEP uses left normalized alleles to indicate freq allele
            # so need to use normalized alleles to match freq allele / consequence allele
            # to the correct variant alt allele
            alleleFreq = None if frequencies is None \
                else self.__get_allele_frequencies(normAlt, frequencies)    
            alleleFreq = self.__add_vcf_frequencies(alleleFreq, vcfEntry.get_info('FREQ'), alt)
                        
            msConseq = self.__vep_parser.get_most_severe_consequence(normAlt)

            recordPK = self.__pk_generator.generate_primary_key(annotator.get_metaseq_id(), variantExternalId)

            copyValues = ['chr' + xstr(variant.chromosome),
                          recordPK,
                          xstr(variant.position),
                          xstr(variant.is_multi_allelic, falseAsNull=True, nullStr='NULL'),
                          binIndex,
                          xstr(variant.ref_snp_id, nullStr='NULL'),
                          annotator.get_metaseq_id(),
                          xstr(annotator.get_display_attributes(variant.rs_position), nulLStr='NULL'),
                          xstr(alleleFreq, nullStr='NULL'),
                          xstr(msConseq, nullStr='NULL'),
                          xstr(self.__vep_parser.get_allele_consequences(normAlt), nullStr='NULL'),
                          xstr(resultJson),
                          xstr(self.__alg_invocation_id)
                          ]           
       
            self.add_copy_str('#'.join(copyValues))

    
    def parse_result(self, result):
        """! parse & load single line from file 
        @params result             the VEP result in string form
        @params checkSkip          make checks to see if line should be skipped (after resume)
        @returns copy string for db load 
        """
        
        if self.resume_load() is False and self.__resume_after_variant is None:
            err = ValueError('Must set VariantLoader result_afer_variant if resuming load')
            self.log(str(err), prefix="ERROR")
            raise err
        
        try:
            self.increment_counter('line')
            resultJson = json.loads(result)
            self.__vep_parser.set_annotation(copy.deepcopy(resultJson))
            
            entry = VcfEntryParser(self.__vep_parser.get('input'))
            
            if not self.resume_load():
                self.__update_resume_status(entry.get('id'))
                return None
            
            # otherwise proceed with the parsing            
            entry.update_chromosome(self.__chromosome_map)
            
            # there are json formatting issues w/the input str
            # so replace w/the parsed entry; which is now a dict
            resultJson['input'] = entry.get_entry()
            
            # extract identifying variant info for frequent reference
            self.__current_variant = entry.get_variant(dbSNP=self.is_dbsnp(), namespace=True)
            
            # rank consequences
            self.__vep_parser.adsp_rank_and_sort_consequences()
    
            # iterate over alleles
            self.__parse_alt_alleles(entry, resultJson)
            
        except:
            1 # print error
    
    
    def load_variants(self):
        """! perform copy operation to insert variants into the DB """ 
        try:
            self.__copy_buffer.seek(0)
            self.__cursor.copy_expert(self.__copy_sql, self.__copy_buffer, 2**10)
            self.reset_copy_buffer()
        except Exception as e:
            err = raise_pg_exception(e, returnError=True)
            self.log(str(err), prefix="ERROR")
            raise err
