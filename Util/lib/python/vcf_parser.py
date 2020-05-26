'''
utils for parsing VCF files
'''
#pylint: disable=line-too-long,invalid-name

from CBILDataCommon.Util.utils import qw, die, warning

class VcfEntryParser(object):
    ''' utils for parse a single line of a vcf file '''

    def __init__(self, entry):
        self._entry = None if entry is None else self._parse_entry(entry)
        

    def _parse_entry(self, inputStr):
        ''' processes the VCF input string and return map
        # VCF
        # CHROM POS     ID        REF ALT    QUAL FILTER INFO
        #X\t605409\trs780063150\tC\tA\t.\t.\tRS=780063150;RSPOS=605409;dbSNPBuildID=144;SSR=0;SAO=0;VP=0x05000088000d000026000100;GENEINFO=SHOX:6473;WGT=1;VC=SNV;U3;INT;CFL;ASP;KGPhase3;CAF=0.9996,0.0003994;COMMON=0;TOPMED=0.99999203618756371,0.00000796381243628
        '''

        fields = qw('chrom pos id ref alt qual filter info')
        values = inputStr.split('\t')
        result = convert_str2numeric_values(dict(zip(fields, values)))

        # now unpack the info field and save as its own
        info = dict(item.split('=') if '=' in item else [item, True] for item in result['info'].split(';'))
        result['info'] = convert_str2numeric_values(info)

        return result


    def get_refsnp(self):
        ''' extract refsnp id from vcf entry dictionary
        (output from parse_vcf_entry)
        '''
        self.__verify_entry()

        if 'rs' in self._entry['id']:
            return self._entry['id']
        if 'RS' in self._entry['info']:
            return 'rs' + str(self._entry['info']['RS'])
        else:
            return None


    def get_entry(self):
        ''' return the entry '''
        return self._entry

    
    def get(self, key):
        ''' get the entry value associated with the key '''
        self.__verify_entry()
        return self._entry[key]


    def get_info(self, key):
        ''' get the INFO value associated with the key '''
        self.__verify_entry()
        return self._entry['info'][key]


    def normalize_alleles(self, ref, alt, snvDivMinus=False):
        ''' remove leftmost alleles that are equivalent between
        the ref & alt alleles; e.g. CAGT/CG <-> AGT/G
        snvDivMinux -- return '-' for SNV deletion when True,
        otherwise return empty string
        '''

        rLength = len(ref)
        aLength = len(alt)

        if rLength == 1 and aLength == 1:
            return (ref, alt)
        
        lastMatchingIndex = - 1
        for i in range(rLength):
            r = ref[i:i + 1]
            a = alt[i:i + 1]
            if r == a:
                lastMatchingIndex = i
            else:
                break

        if lastMatchingIndex >= 0:
            normAlt = alt[lastMatchingIndex + 1:len(alt)]
            if not normAlt and snvDivMinus:
                normAlt = '-'
            normRef = ref[lastMatchingIndex + 1:len(ref)]
            if not normRef and snvDivMinus:
                normRef = '-'
            return (normRef, normAlt)
        

        return (ref, alt)


    def infer_variant_end_location(self, alt):
        ''' infer span of indels/deletions for a 
        specific alternative allele, modeled off 
        GUS Perl VariantAnnotator, see for more info'''
        self.__verify_entry()
        ref = self.get('ref')
        normRef, normAlt = self.normalize_alleles(ref, alt)

        rLength = len(ref)
        aLength = len(alt)

        position = int(self.get('pos'))
        rsPosition = int(self.get_info('RSPOS'))

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
        assert self._entry is not None, \
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
