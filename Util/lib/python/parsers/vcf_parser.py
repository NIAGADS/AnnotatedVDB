"""! @brief VCF Entry Parser"""

##
# @package parsers
# @file vcf_parser.py
#
# @brief  VCF Entry Parser
# 
# @section vcf_parser Description
# utils for parsing VCF entries (a line in a VCF file)
#
# @section todo_vcf_parser TODO
#
# - None
#
# @section libraries_vcf_parser Libraries/Modules
# - types: SimpleNamespace -- allows treatment of a dict as a namespace to access with dot notation
# - [GenomicsDBData.Util.utils](https://github.com/NIAGADS/GenomicsDBData/blob/master/Util/lib/python/utils.py)
#   + provides variety of wrappers for standard file, string, list, and logging operations
# - [GenomicsDBData.Util.list_utils](https://github.com/NIAGADS/GenomicsDBData/blob/master/Util/lib/python/list_utils.py)
#
# @section author_vcf_parser Author(s)
# - Created by Emily Greenfest-Allen (fossilfriend) 2019
# - Modified to remove consequence ranking to the ConsequenceParser class by EGA 2022

# pylint: disable=line-too-long,invalid-name,no-self-use

from types import SimpleNamespace

from GenomicsDBData.Util.utils import xstr, warning, convert_str2numeric_values, to_numeric
from GenomicsDBData.Util.list_utils import qw

class VcfEntryParser(object):
    """! utils for parse a single line of a vcf file """

    def __init__(self, entry, headerFields = None, verbose=False, debug=False, loader=None):
        """! VcfEntryParser base class initializer
        @param entry       VCF Entry (row) in string format
        @param headerFields array of fields if pVCF or non-standard VCF
        @param verbose     flag for verbose output
        @param debug       flag for debug output
        @param loader      if loader is present, can use its logging
        @returns           An instance of the VcfEntryParser class with parsed entry if entry is not None
        """   
        
        self.__debug = debug
        self.__verbose = verbose
        self.__loader = loader
        self._header_fields = qw('chrom pos id ref alt qual filter info') \
            if headerFields is None \
            else [x.lower().replace("#", "") for x in headerFields]
            
        self.__entry = None if entry is None else self.parse_entry(entry)

        
    def debug(self):
        return self.__debug
    
    
    def verbose(self):
        return self.__verbose
    
    
    def log(self, message, prefix="DEBUG"):
        if self.__loader:
            self.__loader.log(message, prefix)
        else:
            if isinstance(message, str):
                warning(prefix, message)
            else:
                warning(prefix, ' '.join((xstr(x) for x in message)))
        

    def parse_entry(self, inputStr):
        """! processes the VCF input string and return map
        
        - Example VCF entry
        
        ~~~~~~~~~~~~~{.txt}
        CHROM POS     ID        REF ALT    QUAL FILTER INFO
        `X\t605409\trs780063150\tC\tA\t.\t.\tRS=780063150;RSPOS=605409;dbSNPBuildID=144;SSR=0;SAO=0;VP=0x05000088000d000026000100;GENEINFO=SHOX:6473;WGT=1;VC=SNV;U3;INT;CFL;ASP;KGPhase3;CAF=0.9996,0.0003994;COMMON=0;TOPMED=0.99999203618756371,0.00000796381243628`
        ~~~~~~~~~~~~~
        
        @param inputStr             the line from the VCF file     
        @returns                    the string parsed into a dict
        @exception                  raises error when run into problems parsing the INFO field
        """
        fields = self._header_fields
        values = inputStr.split('\t')
        result = convert_str2numeric_values(dict(zip(fields, values)))
       
        
        # now unpack the info field and save as its own
        try:
            if 'info' in result:
                infoStr = result['info'].replace('\\x2c', ',') # \x escape causes JSON parsing issues
                infoStr = infoStr.replace('\\x59', '/') # \x59 is a semi-colon, but b/c its a delimiter can't use
                infoStr = infoStr.replace('#', ':') # b/c using pound sign as COPY delimiter
                info = dict(item.split('=',1) if '=' in item else [item, True] for item in infoStr.split(';'))
                result['info'] = convert_str2numeric_values(info)
            else:
                result['info'] = {}
            
        except Exception as err:
            warning("ERROR parsing variant -", result['id'], "- unable to split item in VCF entry INFO field:", xstr(result['info']))
            raise err
        return result


    def update_chromosome(self, chrmMap):
        """! update chromosome in entry based upon provided mapping
        @param chrmMap          ChromosomeMap object mapping sequenceID => chrm number
        @exception TypeError if entry is not set / AttributeError if chromosome number not in the map
        """
        self.__verify_entry()
        if chrmMap is not None:
            self.__entry['chrom'] = chrmMap.get(self.__entry['chrom'])

            
    def get_variant(self, dbSNP=False, namespace=False):
        """! extract basic variant attributes from a VCF entry 
        @param dbSNP                  dbSNP VCF (expect fields that may not be in a generic VCF)
        @param namespace              if True return SimpleNamespace, else return dict
        @returns attributes as a simple namespace so they can be accessed in dot notation"""

        attributes = {}

        chrom = xstr(self.get('chrom'))
        if chrom == 'MT':
            chrom = 'M'
        altAlleles = self.get('alt').split(',')
        id = self.get('id')
        if id == '.':
            id = ':'.join((chrom.replace('chr', ''), xstr(self.get('pos')), self.get('ref'), self.get('alt')))
            
        variant =  {
            'id' : id,
            'ref_snp_id' : self.get_refsnp(),
            'ref_allele' : self.get('ref'),
            'alt_alleles' : altAlleles,
            'is_multi_allelic' : len(altAlleles) > 1,
            'chromosome' : xstr(chrom).replace('chr', ''),
            'position' : int(self.get('pos')),
            'rs_position' : self.get_info('RSPOS') 
        }
        
        return SimpleNamespace(**variant) if namespace else variant
    

    def get_refsnp(self):
        """! extract refsnp id from vcf entry dictionary
        @returns ref_snp_id 
        """
        self.__verify_entry()

        if 'rs' in self.__entry['id']:
            return self.__entry['id']
        if 'RS' in self.__entry['info']:
            return 'rs' + str(self.__entry['info']['RS'])
        else:
            return None


    def get_entry(self):
        """! return the parsed entry 
        @returns the parsed entry """
        return self.__entry

    
    def get(self, key):
        """! get the entry value associated with the key """
        self.__verify_entry()
        return self.__entry[key]


    def get_info(self, key):
        """! get the INFO value associated with the key """
        self.__verify_entry()
        if 'info' not in self.__entry:
            return None
        if key in self.__entry['info']:
            return self.__entry['info'][key]
        else:
            return None
        
        
    def get_frequencies(self, allele):
        """! retrieve allele frequencies reported in INFO FREQ field
        @param allele          allele to match (not normalized)
        @param vcfGMAFs         global minor allele frequencies from the VCF FREQ info field
        @returns               dict of source:frequency for the allele
        """
        
        vcfGMAFs = self.get_info('FREQ')
        if vcfGMAFs is None:
            return None
        
        zeroValues = ['.', '0']

        # FREQ=GnomAD:0.9986,0.001353|Korea1K:0.9814,0.01861|dbGaP_PopFreq:0.9994,0.0005901
        # altIndex needs to be incremented as first value is for the ref allele)
        altAlleles = self.get('alt').split(',')
        altIndex = altAlleles.index(allele) + 1
        populationFrequencies = {pop.split(':')[0]:pop.split(':')[1] for pop in vcfGMAFs.split('|')}
        vcfFreqs = {pop: {'gmaf': to_numeric(freq.split(',')[altIndex])} \
                    for pop, freq in populationFrequencies.items() \
                    if freq.split(',')[altIndex] not in zeroValues}
    
        return None if len(vcfFreqs) == 0 else vcfFreqs



    def infer_variant_end_location(self, alt, normRef):
        """! infer span of indels/deletions for a 
        specific alternative allele, modeled off 
        GUS Perl VariantAnnotator & dbSNP normalization conventions
        @param alt               alternative allele
        @param normRef           left normalized reference allele
        @returns                 end location
        """
        self.__verify_entry()
        ref = self.get('ref')

        rLength = len(ref)
        aLength = len(alt)

        position = int(self.get('pos'))
        rsPosition = self.get_info('RSPOS')
        if rsPosition is None:
            rsPosition = position
        else:
            rsPosition = int(rsPosition)

        if rLength == 1 and aLength == 1: # SNV
            return position

        if rLength == aLength: # MNV
            if ref == alt[::-1]: #inversion
                return position + rLength - 1

            # substitution
            return position + len(normRef) - 1

        if rLength > aLength: # deletions
            if len(alt) > 1: # indel
                if len(normRef) == 0: # was normalized; adjust
                    return rsPosition + len(ref) - 2
                return rsPosition + len(ref) - 1
            else: # straight up deletion
                return rsPosition + len(normRef) -  1

        if rLength < aLength: # insertion
            return rsPosition + 1


    def __verify_entry(self):
        """! check that entry is set 
        @returns boolean if entry has a value 
        @exception   failed assertion """
        assert self.__entry is not None, \
            "DEBUG - must set value of _entry in the VCF parser before attempting to access"

