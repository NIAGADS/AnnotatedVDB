'''
utils for parsing VCF files
'''
#pylint: disable=line-too-long,invalid-name

from types import SimpleNamespace

from GenomicsDBData.Util.utils import xstr, die, warning, print_dict, print_args
from GenomicsDBData.Util.list_utils import qw

class VcfEntryParser(object):
    '''! utils for parse a single line of a vcf file '''

    def __init__(self, entry, verbose=False, debug=False):
        """! VcfEntryParser base class initializer
        @param entry       VCF Entry (row) in string format
        @param verbose     flag for verbose output
        @param debug       flag for debug output
        @returns           An instance of the VcfEntryParser class with parsed entry if entry is not None
        """
        
        self.__entry = None if entry is None else self.__parse_entry(entry)
        

    def __parse_entry(self, inputStr):
        ''' processes the VCF input string and return map
        # VCF
        # CHROM POS     ID        REF ALT    QUAL FILTER INFO
        #X\t605409\trs780063150\tC\tA\t.\t.\tRS=780063150;RSPOS=605409;dbSNPBuildID=144;SSR=0;SAO=0;VP=0x05000088000d000026000100;GENEINFO=SHOX:6473;WGT=1;VC=SNV;U3;INT;CFL;ASP;KGPhase3;CAF=0.9996,0.0003994;COMMON=0;TOPMED=0.99999203618756371,0.00000796381243628
        '''

        fields = qw('chrom pos id ref alt qual filter info')
        values = inputStr.split('\t')
        result = convert_str2numeric_values(dict(zip(fields, values)))

        # now unpack the info field and save as its own
        try:
            info = dict(item.split('=',1) if '=' in item else [item, True] for item in result['info'].split(';'))
        except Exception as err:
            warning("ERROR parsing variant -", result['id'], "- unable to split item in VCF entry INFO field:", xstr(result['info']))
            raise err
        
        result['info'] = convert_str2numeric_values(info)

        return result


    def update_chromosome(self, chrmMap):
        '''! update chromosome in entry based upon provided mapping
        @param chrmMap          dict mapping sequence_id => chrm number
        @raises TypeError if entry is not set
        '''
        self.__verify_entry()
        if chrmMap is not None:
            self.__entry['chrom'] = chrmMap[self.__entry['chrom']]

            
    def get_variant(self, dbSNP = False, namespace=False):
        """! extract basic variant attributes from a VCF entry 
        @param dbSNP                  dbSNP VCF (expect fields that may not be in a generic VCF)
        @param namespace              if True return SimpleNamespace, else return dict
        @returns attributes as a simple namespace so they can be accessed in dot notation"""

        attributes = {}
        if dbSNP:
            chrom = xstr(self.get('chrom'))
            if chrom == 'MT':
                chrom = 'M'
            altAlleles = self.get('alt').split(',')
            variant =  {
                'id' : self.get('id'),
                'ref_snp_id' : self.get_refsnp(),
                'ref_allele' : self.get('ref'),
                'alt_alleles' : altAlleles,
                'is_multi_allelic' : len(altAlleles) > 1,
                'chromosome' : xstr(chrom),
                'position' : int(self.get('pos')),
                'rs_position' : self.get_info('RSPOS'),
            }
        else:
            err = NotImplementedError('VcfEntryParser.get_variant not implemented for non-dbSNP VCF variants')
            raise err
        
        return SimpleNamespace(**variant) if namespace else variant
    
    
    def get_refsnp(self):
        ''' extract refsnp id from vcf entry dictionary
        (output from parse_vcf_entry)
        '''
        self.__verify_entry()

        if 'rs' in self.__entry['id']:
            return self.__entry['id']
        if 'RS' in self.__entry['info']:
            return 'rs' + str(self.__entry['info']['RS'])
        else:
            return None


    def get_entry(self):
        ''' return the entry '''
        return self.__entry

    
    def get(self, key):
        ''' get the entry value associated with the key '''
        self.__verify_entry()
        return self.__entry[key]


    def get_info(self, key):
        ''' get the INFO value associated with the key '''
        self.__verify_entry()
        if key in self.__entry['info']:
            return self.__entry['info'][key]
        else:
            return None
        

    def infer_variant_end_location(self, alt, normRef):
        ''' infer span of indels/deletions for a 
        specific alternative allele, modeled off 
        GUS Perl VariantAnnotator, see for more info'''
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
        ''' check that entry is set '''
        assert self.__entry is not None, \
          "DEBUG - must set value of _entry in the VCF parser before attempting to access"



# --------------
# helpers


def is_integer(value):
    if isinstance(value, (float, bool)):
        return False
    try:
        int(value)
        return True
    except ValueError:
        return False


def is_float(value):
    try:
        int(value)
        return True
    except ValueError:
        return False


def convert_str2numeric_values(cdict):
    '''
    converts numeric values in dictionary stored as strings 
    to numeric
    '''

    for key, value in cdict.items():
        if is_float(value): # must check float first b/c integers are a subset
            cdict[key] = float(value)
        if is_integer(value):
            cdict[key] = int(value)

    return cdict
