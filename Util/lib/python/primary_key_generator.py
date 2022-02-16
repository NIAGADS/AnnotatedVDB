"""! @brief GenomicsDB Record Primary Key Generators"""

##
# @file primary_key_generator.py
#
# @brief  GenomicsDB Record Primary Key Generators
#
# @section primary_key_generator Description
# Defines primary key generators for GenomicsDB records:
## - VariantPKGenerator: generator for variants
#
# @section todo_primary_key_generator TODO
# - evaluate VRS normalization
# - handle EVA sub SNP ids (ss) ids?
# - validate GenomeBuild?
# - should chrM be 'M' or 'MT' to be mapped correctly?
#
# @section libraries_primary_key_generator Libraries/Modules
# - [vrs-python](https://github.com/ga4gh/vrs-python)
#   + provides Python language support for the GA4GH Variation Representation Specification (VRS)
# - [GenomicsDBData.Util.utils](https://github.com/NIAGADS/GenomicsDBData/blob/master/Util/lib/python/utils.py)
#   + provides variety of wrappers for standard file, string, list, and logging operations
#
# @section author_primary_key_generator Author(s)
# - Created by Emily Greenfest-Allen (fossilfriend) 2022

from ga4gh.core import ga4gh_identify, ga4gh_serialize
from ga4gh.vrs.extras.translator import Translator
from ga4gh.vrs.dataproxy import create_dataproxy
from GenomicsDBData.Util.utils import warning, die, xstr, print_dict

class VariantPKGenerator(object):
    """! Generator for variant record primary key.
    
    Keys are generated using the following specifications:
    - [NCBI dbSNP SPDI or Sequence:Position:Deletion:Insertion](https://www.ncbi.nlm.nih.gov/variation/notation/)
    - [GA4GH VRS Computed Sequence Representations](https://vrs.ga4gh.org/en/stable/impl-guide/computed_identifiers.html#computed-identifiers)

    According to the following rules:
    - for SNV's and short INDEX (total length alleles <= maxSequenceLength):
      + S:P:D:I, with restriction that deletion (reference) length should not substituted for the sequence
    - for large INDELS/SVs (total length alleles > maxSequenceLength)
      + S:P:VRS_CI
    - if an external id is provided (e.g., refSNP or ss) it is provided at the end so that PKs are as follows:
      + S:P:D:I:refSNP
      + S:P:VRS_CI:refSNP
   """
    

    def __init__(self, genomeBuild, seqrepoProxyPath, maxSequenceLength=50, normalize=False, verbose=False, debug=False):
        """! VariantPKGenerator base class initializer
        @param genomeBuild          Assembly name (e.g., GRCh38, GRCh37)
        @param seqrepoProxyPath     full path to the file-based seqrepo data repository / required by GA4GH VRS
        @param maxSequenceLength    max length for ref & alt alleles, above which sequence should be digested
        @param normalize            apply GA4GH normalization
        @param verbose              verbose output flag
        @param debug                debug flag

        @return                     An instance of the VariantPKGenerator class with initialized translator
        """
        
        self._verbose = verbose
        self._debug = debug
        self._genomeBuild = genomeBuild
        self._maxSequenceLength = maxSequenceLength
        self._ga4gh_sequence_map = {}
        self._translator = None
        self._set_seqrepo_translator(seqrepoProxyPath, normalize)


    def _set_seqrepo_translator(self, proxyPath, normalize=False):
        """! set the seqrepo translator

        @param proxyPath         full path to the file-based seqrepo data repository
        @param normalize         apply GA4GH normalization
        """

        dataProxy = create_dataproxy("seqrepo+file://" + proxyPath)
        self._translator = Translator(data_proxy=dataProxy)
        self._translator.normalize = normalize



    def translate_vrs(self, vrsDict, formatSpec="spdi"):
        """! translate VRS allele dict back to HGVS or SPDI format / for validation purposes
        @param vrsDict           the VRS Allele Dict object
        @param formatSpec        target format: hgvs or spdi
        @returns                 translation of the VRS Allele Dict
        """
        
        if formatSpec not in ['hgvs', 'spdi']:
            raise IndexError("invalid format spec: specify 'spdi' or 'hgvs'")

        return self._translator.translate_to(vrsDict, formatSpec)

        
    def generate_primary_key(self, metaseqId, externalId=None):
        """! generate and returns the primary key 
        @param metaseqId         metaseq or SPDI formatted (with deletion sequence) variant representation
        @param externalId        refSnp or subSnp ID
        @returns                 generated primary key
        """

        chrm, position, ref, alt = metaseqId.split(':')

        pk = [chrm, position]
        longSequence = False
        if len(ref) + len(alt) <= self._maxSequenceLength:
            pk.extend([ref,alt])
        else:
            pk.append(self.compute_vrs_identifier(metaseqId))

        if externalId is not None:
          pk.append(externalId)
            
        return ':'.join(pk)
        

    def get_vrs_allele_dict(self, metaseqId, serialize=False, toJson=False):
        """! get the GA4GH VRS Allele dict, given a variant
        @param metaseqId          metaseq or SPDI formatted (with deletion sequence) variant representation
        @param serialize          serialize to binary object
        @param toJson             return as JSON object

        @returns                  VRS Allele dict in format specified by flags
        """
        
        # transform metaseq id to gnomad id (replace ':' with '-')
        # using gnomad b/c it accepts chrNum instead of refseq accession (required by SPDI)
        gnomadExpr = metaseqId.replace(':', '-')
        alleleDict = self._translator._from_gnomad(gnomadExpr, self._genomeBuild)
      
        if serialize:
            return ga4gh_serialize(alleleDict)
        if toJson:
            return alleleDict.for_json() 

        return alleleDict
        
        
    def compute_vrs_identifier(self, metaseqId):
        """! return computed GA4GH identifier for the variant
        @param metaseqId          metaseq or SPDI formatted (with deletion sequence) variant representation
        @returns                  VRS Computed Identifier with ga4gh:VA prefix removed
        """
        alleleDict = self.get_vrs_allele_dict(metaseqId)
        
        if self._debug:
            debugOutput = {
                "Input Variant":  metaseqId,
                "VRS Representation" : alleleDict.for_json(),
                }
            warning(print_dict(debugOutput, pretty=True))

        vrsComputedId = ga4gh_identify(alleleDict)
            
        return vrsComputedId.split('.')[1]


        
        

        
