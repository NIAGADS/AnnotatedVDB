""" variant annotator functions """
#!pylint: disable=invalid-name

from GenomicsDBData.Util.utils import xstr, truncate, reverse, warning
from GenomicsDBData.Util.list_utils import qw


BASE_LOAD_FIELDS = qw('chromosome record_primary_key position is_multi_allelic bin_index ref_snp_id metaseq_id display_attributes allele_frequencies adsp_most_severe_consequence adsp_ranked_consequences vep_output row_algorithm_id', returnTuple=True)

def truncate_allele(value):
    ''' wrapper for trunctate to 5 chars '''
    return truncate(value, 5)


class VariantAnnotator(object):
    ''' functions used in different loading scripts / consolidated '''

    def __init__(self, refAllele, altAllele, chrom, position):
        self._ref = refAllele
        self._alt = altAllele
        self._chrom = chrom
        self._position = position
        self._metaseqId = None
        self.__set_metaseq_id()
        

    def get_normalized_alleles(self, snvDivMinus=False):
        ''' public return normalized alleles '''
        return self.__normalize_alleles(snvDivMinus)
    

    def __normalize_alleles(self, snvDivMinus=False):
        ''' remove leftmost alleles that are equivalent between
        the ref & alt alleles; e.g. CAGT/CG <-> AGT/G
        
        if no normalization is possible, keeps normalized alleles as
        default equal to ref & alt

        snvDivMinux -- return '-' for SNV deletion when True,
        otherwise return empty string
        returns a ref, alt tuple
        '''

        rLength = len(self._ref)
        aLength = len(self._alt)

        if rLength == 1 and aLength == 1: # SNV no normalization needed
            return self._ref, self._alt

        else:
            lastMatchingIndex = - 1
            for i in range(rLength):
                r = self._ref[i:i + 1]
                a = self._alt[i:i + 1]
                if r == a:
                    lastMatchingIndex = i
                else:
                    break

            if lastMatchingIndex >= 0: 
                normAlt = self._alt[lastMatchingIndex + 1:len(self._alt)]
                if not normAlt and snvDivMinus:
                    normAlt = '-'
                normRef = self._ref[lastMatchingIndex + 1:len(self._ref)]
                if not normRef and snvDivMinus:
                    normRef = '-'

                return normRef, normAlt
                    
            else: # MNV no normalization needed
                return self._ref, self._alt

        return self._ref, self._alt # not sure under what conditions this would occur, but covering bases
            
        
    def __set_metaseq_id(self):
        ''' generate metaseq id '''
        warning(self._chrom, self._position, self._ref, self._alt)
        self._metaseqId = ':'.join((self._chrom, xstr(self._position), self._ref, self._alt))

        
    def get_metaseq_id(self):
        ''' return metaseq id '''
        return self._metaseqId


    def get_display_attributes(self, rsPosition = None):
        ''' 
        generate display alleles & dbSNP compatible start-end 
        rsPosition = dbSNP property RSPOS 
        '''

        if rsPosition is None:
            rsPosition = self._position

        normRef, normAlt = self.__normalize_alleles() # accurate length version
        nRefLength = len(normRef)
        nAltLength = len(normAlt)
        normRef, normAlt = self.__normalize_alleles(True) # display version (- for empty string)

        refLength = len(self._ref);
        altLength = len(self._alt)
        
        attributes = {
            'location_start': self._posiion,
            'location_end': self._position
        }
            
        if (refLength == 1 and altLength == 1): # SNV
            attributes.update({
                'variant_class': "single nucleotide variant",
                'variant_class_abbrev': "SNV",
                'display_allele': self._ref + '>' + self._alt,
                'sequence_allele': self._ref + '/' + self._alt
            })

        elif refLength == altLength: # MNV
            #inversion
            if self._ref == reverse(self._alt):
                attributes.update({
                    'variant_class': "inversion",
                    'variant_class_abbrev': "MNV",
                    'display_allele': 'inv' + self._ref,
                    'sequence_allele': truncate_allele(self._ref) + '/' + truncate_allele(self._alt),
                    'location_end': self._position + refLength - 1
                })
            else:
                attributes.update({
                    'variant_class': "substitution",
                    'variant_class_abbrev': "MNV",
                    'display_allele': normRef + ">" + normAlt,
                    'sequence_allele': truncate_allele(normRef) + '/' + truncate_allele(normAlt),
                    'location_start': rsPosition,
                    'location_end': rsPosition + nRefLength - 1
                })
        # end MNV
        
        elif (refLength > altLength): # deletions
            attributes.update({'location_start': rsPosition})

            if nAltLength > 1: # INDEL
                attributes.update({
                    'variant_class': 'indel',
                    'variant_class_abbrev': 'INDEL'
                    })

                if nRefLength == 0:
                    displayRef = self._ref[1:] # strip first character from reference
                    attributes.update({
                        'location_end': rsPosition + refLength - 1,
                        'display_allele': "del" + displayRef + "ins" + normAlt,
                        'sequence_allele': truncate_allele(displayRef) + "/" + truncate_allele(normAlt)
                        })
                else:
                    attributes.update({
                        'location_end': rsPosition + nRefLength - 1,
                        'display_allele': "del" + normRef + "ins" + normAlt,
                        'sequence_allele': truncate_allele(normRef) + "/" + truncate_allele(normAlt)
                        })
                # end INDEL
                
            else: # deletion
                attributes.update({
                    'variant_class': "deletion",
                    'variant_class_abbrev': "DEL",
                    'location_end': rsPosition + nRefLength - 1,
                    'display_allele': "del" + normRef,
                    'sequence_allele': truncate_allele(normRef) + "/-"
                    })
                    
        # end deletions
        
        elif refLength < altLength: # insertions
            attributes.update({'location_start': rsPosition})

            if refLength > 1: # INDEL
                attributes.update({
                    'variant_class': 'indel',
                    'variant_class_abbrev': 'INDEL'
                    })
                
                if nRefLength == 0:
                    displayRef = self._ref[1:] # strip first character from reference
                    attributes.update({
                        'location_end': rsPosition + refLength - 1,
                        'display_allele': "del" + displayRef + "ins" + normAlt,
                        'sequence_allele': truncate_allele(displayRef) + "/" + truncate_allele(normAlt)
                        })
                else:
                    attributes.update({
                        'location_end': rsPosition + nRefLength - 1,
                        'display_allele': "del" + normRef + "ins" + normAlt,
                        'sequence_allele': truncate_allele(normRef) + "/" + truncate_allele(normAlt)
                        })
                # end INDEL
            else: # insertion
                attributes.update({
                    'variant_class': "insertion",
                    'variant_class_abbreve': "INS",
                    'location_end': rsPosition + 1,
                    'display_allele': "ins" + normAlt,
                    'sequence_allele': "-/" + truncate_allele(normAlt)
                })
                
        return attributes
