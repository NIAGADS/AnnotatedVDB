"""! @brief Annotated VDB Loader """

#!pylint: disable=invalid-name

##
# @package loaders
# @file variant_loader.py
#
# @brief  Annotated VDB Loader
#
# @section variant_loader Description
# provides functions & class for managing variant loads
# - types
#   + VariantLoader (parent class)
#   + VEPVariantLoader (load from VEP result)
#   + VCFVariantLoader (load from VCF)
#
# @section todo_variant_loader TODO
# - create __verify_current_variant() to make sure __current_variant is set before trying to access
# - extract common elements to parent class and make VEPLoader a child
#
# @section libraries_variant_loader Libraries/Modules
# - json - parse/serialize JSON
# - copy - for deep copy functionality 
# - [GenomicsDBData.Util.utils](https://github.com/NIAGADS/GenomicsDBData/blob/master/Util/lib/python/utils.py)
#   + provides variety of wrappers for standard file, string, and logging operations
# - [GenomicsDBData.Util.list_utils](https://github.com/NIAGADS/GenomicsDBData/blob/master/Util/lib/python/list_utils.py)
#   + provides variety list and set operations
# - [GenomicsDBData.Util.postgres_dbi](https://github.com/NIAGADS/GenomicsDBData/blob/master/Util/lib/python/postgres_dbi.py)
#   + wrapper for connecting to database & handling PG errors
# - [AnnotatedVDB.BinIndex](https://github.com/NIAGADS/AnnotatedVDB/blob/master/BinIndex/lib/python/bin_index.py)
#   + functions to calculating bin index of a sequence feature
# - [AnnotatedVDB.Util.algorithm_invocation](https://github.com/NIAGADS/AnnotatedVDB/blob/master/Util/lib/python/algorithm_invocation.py)
#   + create entry in AnnotatedVDB.AlgorithmInvocation to track load & allow for UNDO
# - [AnnotatedVDB.Util.parsers](https://github.com/NIAGADS/AnnotatedVDB/blob/master/Util/lib/python/parsers)
#   + parsers for VCF and VEP JSON
#
# @section author_variant_loader Author(s)
# - Created by Emily Greenfest-Allen (fossilfriend) 2022

import json 

from copy import deepcopy
from io import StringIO

from GenomicsDBData.Util.utils import xstr, warning, print_dict, print_args, to_numeric, deep_update
from GenomicsDBData.Util.list_utils import qw, is_subset, is_equivalent_list
from GenomicsDBData.Util.postgres_dbi import raise_pg_exception

from AnnotatedVDB.Util.variant_annotator import VariantAnnotator
from AnnotatedVDB.Util.primary_key_generator import VariantPKGenerator
from AnnotatedVDB.Util.parsers import VcfEntryParser, VepJsonParser

from AnnotatedVDB.BinIndex.bin_index import BinIndex
from AnnotatedVDB.Util.algorithm_invocation import AlgorithmInvocation
from AnnotatedVDB.Util.database import VariantRecord

ALLOWABLE_COPY_FIELDS = ["chromosome", "record_primary_key", "position", 
                         "is_multi_allelic", "is_adsp_variant", "ref_snp_id", 
                         "metaseq_id", "bin_index", "display_attributes", 
                         "allele_frequencies", "cadd_scores", 
                         "adsp_most_severe_consequence", "adsp_ranked_consequences", 
                         "loss_of_function", "vep_output", "adsp_qc", 
                         "gwas_flags", "other_annotation", "row_algorithm_id"]

REQUIRED_COPY_FIELDS = ["chromosome", "record_primary_key", "position", "metaseq_id", "bin_index", "row_algorithm_id"]

DEFAULT_COPY_FIELDS = qw('chromosome record_primary_key position is_multi_allelic bin_index ref_snp_id metaseq_id display_attributes allele_frequencies adsp_most_severe_consequence adsp_ranked_consequences vep_output row_algorithm_id', returnTuple=True)

class VariantLoader(object):
    """! functions for loading variants -- use child classes for specific datasources / result types """
    
    def __init__(self, datasource, logFileName=None, verbose=False, debug=False):
        """! VariantLoader base class initializer

            @param datasource          datasource description
            @param logFileName         full path to logging file
            @param verbose             flag for verbose output
            @param debug               flag for debug output
            
            @returns                   An instance of the VariantLoader class with initialized log, counters and copy buffer
        """
        # single underscore so protected / child classes can access
        self._log_file_handle = None if logFileName is None else open(logFileName, 'w') 
        self._verbose = verbose
        self._debug = debug
        self._datasource = datasource.lower() 
        
        self._alg_invocation_id = None
        self._pk_generator = None
        self._bin_indexer = None
        
         ## database cursor
        self._cursor = None
        
        self._counters = {}
        self._chromosome_map = None
        self._resume_after_variant = None
        self._resume = True
        self._current_variant = {}
        
        self._copy_buffer = None
        self._copy_fields = None
        self._copy_sql = None
        
        self._fail_at_variant = None
        self._variant_validator = None
        self._skip_duplicates = False

        self._initialize_counters()
        self.initialize_copy_buffer()
        

    def close(self):
        self.close_copy_buffer()
        self.close_log()
        if self._variant_validator is not None:
            self._variant_validator.close()
    
        
    def _exit__(self, exc_type, exc_value, traceback):
        self.close()
        
    
    def initialize_variant_validator(self, gusConfigFile):
        self._variant_validator=VariantRecord(gusConfigFile, self._verbose, self._debug)
    
    
    def set_skip_duplicates(self, skipDuplicates, gusConfigFile=None):
        self._skip_duplicates = skipDuplicates
        if skipDuplicates:
            self.initialize_variant_validator(gusConfigFile)
            
            
    def is_duplicate(self, variantId, idType, chromosome=None):
        return self._variant_validator.exists(variantId, idType, chromosome)
    
    
    def skip_duplicates(self):
        return self._skip_duplicates
    
    
    def variant_validator(self):
        return self._variant_validator()
    
    
    def is_fail_at_variant(self):
        """! check if current variant is fail at variant
            @returns                 True if match
        """

        return self._fail_at_variant == self._current_variant.id
    
    
    def fail_at_variant(self):
        """! @returns fail at variant """
        return self._fail_at_variant
    
    
    def set_fail_at_variant(self, variant):
        """! for debugging problematic variants that fail on load
            @param variant             id of variant that should be logged before exit
        """
        self._fail_at_variant = variant


    def set_chromosome_map(self, chrmMap):
        self._chromosome_map = chrmMap


    def get_current_variant_id(self):
        if self._debug:
            self.log(print_dict(self._current_variant, pretty=True), prefix="DEBUG")
        return self._current_variant.id


    def _verify_copy_fields(self, copyFields):
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
        self._verify_copy_fields(copyFields)   
        self._copy_fields = DEFAULT_COPY_FIELDS if copyFields is None else copyFields
        
        self._copy_sql = "COPY AnnotatedVDB.Variant("  \
         + ','.join(self._copy_fields) \
         + ") FROM STDIN WITH (NULL 'NULL', DELIMITER '#')"


    def close_copy_buffer(self):
        """! close the copy buffer """
        self._copy_buffer.close()
        
        
    def reset_copy_buffer(self):
        """! reset copy buffer to clean up memory """    
        self.close_copy_buffer()
        self.initialize_copy_buffer()
    
    
    def initialize_copy_buffer(self):
        """! initialize the copy buffer (io.StringIO) """
        self._copy_buffer = StringIO()
        
        
    def copy_buffer(self, sizeOnly=False):
        """! returns copy buffer
        @param sizeOnly               return size only 
        @returns copy buffer or buffer size """
        
        return self._copy_buffer.tell() if sizeOnly else self._copy_buffer
    
    def add_copy_str(self, copyStr):
        """! write a copy string to the copy buffer """
        self._copy_buffer.write(copyStr + '\n')


    def is_dbsnp(self):
        """!  @returns  flag indicating if datasource is dbsnp
        """
        return self._datasource == 'dbsnp'
    

    def resume_load(self):
        """! check if resume load flag has been toggled
        @returns    boolean stored in __resume
        """
        return self._resume
    

    def set_resume_after_variant(self, variantId):
        """! resume after specified variant
        @param variantId           id (metaseq or refsnp) of the variant
        """
        self._resume_after_variant = variantId
        self._resume = False # skipping until variant is found
        
    
    def initialize_bin_indexer(self, gusConfigFile):
        """! initialize Bin Indexer 

            @param gusConfigFile             gus.config file w/database connection info
        """
        self._bin_indexer = BinIndex(gusConfigFile, verbose=False)


    def initialize_pk_generator(self, genomeBuild, seqrepoProxyPath):
        """! initialize primary key generator

            @param genomeBuild             genome build (GRCh38 or GRCh37)
            @param seqrepoProxyPath        full path to the seqrepo file proxy
        """
        self._pk_generator = VariantPKGenerator(genomeBuild, seqrepoProxyPath)
 
    
    def get_count(self, counter):
        """! get value of a counter
            @param counter           counter whose value should be returned
            @returns                 value of counter
        """
        return self._counters[counter]
    
    
    def increment_counter(self, counter):
        """! increment counter by 1 """
        self._counters[counter] = self._counters[counter] + 1   
        
    
    def _initialize_counters(self):
        """! initialize counters"""
        self._counters = { 'line' : 0, 'variant': 0, 'skipped': 0, 'duplicates': 0}
            
                
    def cursor(self):
        """! @returns database cursor """    
        return self._cursor
    
    
    def set_cursor(self, cursor):
        """! set database handle
        @param  cursor              database cursor
        """
        self._cursor = cursor
            
        
    def current_variant(self):
        """! @returns current variant """
        return self._current_variant
    
    
    def bin_indexer(self):
        """! @returns Bin indexer """
        return self._bin_indexer
        
        
    def alg_invocation_id(self):
        """! @returns algorithm invocation id """
        return self._alg_invocation_id
    
    
    def pk_generators(self):
        """! @returns primary key generatory"""
        return self._pk_generator
        
    
    def set_algorithm_invocation(self, callingScript, comment, commit=True):
        """! create entry in AnnotatedVDB.AlgorithmInvocation table for the data load
         save alg_invocation_id
        """
        algInvocation = AlgorithmInvocation(callingScript, comment, commit)
        self._alg_invocation_id = xstr(algInvocation.getAlgorithmInvocationId())
        algInvocation.close()

        
    def log(self, message, prefix=None):
        """! print to log
            @param message             message to print to log -- string or tuple/list
            @param prefix              e.g., DEBUG, WARNING, ERROR
        """
        if not isinstance(message, str):
            message = ' '.join([xstr(x) for x in message])
        
        if prefix:
            message = prefix + ": " + message          
        warning(message, file=self._log_file_handle, flush=True)
        
        
    def log_fh(self):
        """! @returns log file handle """
        return self._log_file_handle
    
    
    def close_log(self):
        """! closes log file handle """
        self._log_file_handle.close()
    
    
    def _update_resume_status(self, variantId):
        """! check to see if the row should be skipped 
            update __resume flag accordingly
            
            @param variantId            id field in vep output  
        """
        if not self.resume_load():
            self.increment_counter('skipped')
            self._resume = variantId == self._resume_after_variant
            
            if self.resume_load() is True: # resume load when next line is parsed
                message = ("Resuming after", self._resume_after_variant)
                self.log(message, prefix="WARNING")
                message = ("Skipped", xstr(self.get_count('skipped'), "variants"))
                self.log(message, prefix="INFO")


    def load_variants(self):
        """! perform copy operation to insert variants into the DB """ 
        try:
            self._copy_buffer.seek(0)
            self._cursor.copy_expert(self._copy_sql, self._copy_buffer, 2**10)
            self.reset_copy_buffer()
        except Exception as e:
            err = raise_pg_exception(e, returnError=True)
            self.log(str(err), prefix="ERROR")
            raise err
        
        
    def parse_result(self, result):
        """! parse & load single line from file 
        @param result             the line from the result file
        @returns copy string for db load 
        """
                
        err = NotImplementedError('Result parsing is not defined for the VariantLoader parent class, please use result-specific loader')
        self.log(str(err), prefix="ERROR")
        raise err
      
    

