""" variant annotator functions """
#!pylint: disable=invalid-name
##
# @package util

from GenomicsDBData.Util.utils import xstr, truncate, reverse, warning

def truncate_allele(value):
    """ wrapper for trunctate to 5 chars """
    return truncate(value, 8)

def reverse_complement(seq):
    """! @returns reverse complement of the specified sequence (seq)
    """
    mapping = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(mapping)[::-1]
    
class VariantAnnotator(object):
    """ functions used to generate variant annotations """

    def __init__(self, refAllele, altAllele, chrom, position):
        self.__ref = refAllele
        self.__alt = altAllele
        self.__chrom = chrom
        self.__position = position
        self.__metaseqId = None
        self.__set_metaseq_id()


    def get_normalized_alleles(self, snvDivMinus=False):
        """! public wrapper for __normalize_alleles / LEGACY
        @returns normalized alleles """
        return self.__normalize_alleles(snvDivMinus)
    
    
    def infer_variant_end_location(self, rsPosition=None):
        """! infer span of indels/deletions for a 
        specific alternative allele, modeled off 
        GUS Perl VariantAnnotator & dbSNP normalization conventions
        @returns                 end location
        """

        ref = self.__ref
        alt = self.__alt
        
        normRef, normAlt = self.get_normalized_alleles()

        rLength = len(ref)
        aLength = len(alt)

        position = int(self.__position)

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
                    return position + len(ref)
                return position + len(normRef) 
            else: # straight up deletion
                return position + len(normRef)

        if rLength < aLength: # insertion
            return position + 1


    def __normalize_alleles(self, snvDivMinus=False):
        """! left normalize VCF alleles
        - remove leftmost alleles that are equivalent between the ref & alt alleles;
        e.g. CAGT/CG <-> AGT/G

        - if no normalization is possible, keeps normalized alleles as
        default equal to ref & alt

        @param snvDivMinus       return '-' for SNV deletion when True, otherwise return empty string
        @return                  a tuple containing the normalized alleles (ref, alt)
        """

        rLength = len(self.__ref)
        aLength = len(self.__alt)
        
        if rLength == 1 and aLength == 1: # SNV no normalization needed
            return self.__ref, self.__alt

        lastMatchingIndex = - 1
        for i in range(rLength):
            r = self.__ref[i:i + 1]
            a = self.__alt[i:i + 1]
            if r == a:
                lastMatchingIndex = i
            else:
                break

        if lastMatchingIndex >= 0:
            normAlt = self.__alt[lastMatchingIndex + 1:len(self.__alt)]
            if not normAlt and snvDivMinus:
                normAlt = '-'

            normRef = self.__ref[lastMatchingIndex + 1:len(self.__ref)]
            if not normRef and snvDivMinus:
                normRef = '-'

            return normRef, normAlt

        # MNV no normalization needed
        return self.__ref, self.__alt


    def __set_metaseq_id(self):
        """! generate metaseq id and set value"""
        self.__metaseqId = ':'.join((self.__chrom, xstr(self.__position), self.__ref, self.__alt))


    def get_metaseq_id(self):
        """! @returns metaseq id """
        return self.__metaseqId


    def get_display_attributes(self, rsPosition = None):
        """! generate and return display alleles & dbSNP compatible start-end
        @param rsPosition       dbSNP property RSPOS
        @returns                 dict containing display attributes
        """
        
        position = self.__position
            
        refLength = len(self.__ref)
        altLength = len(self.__alt)

        normRef, normAlt = self.__normalize_alleles() # accurate length version
        nRefLength = len(normRef)
        nAltLength = len(normAlt)
        normRef, normAlt = self.__normalize_alleles(True) # display version (- for empty string)
        
        endLocation = self.infer_variant_end_location()

        attributes = {
            'location_start': position,
            'location_end': position
        }
     
        if (refLength == 1 and altLength == 1): # SNV
            attributes.update({
                'variant_class': "single nucleotide variant",
                'variant_class_abbrev': "SNV",
                'display_allele': self.__ref + '>' + self.__alt,
                'sequence_allele': self.__ref + '/' + self.__alt
            })

        elif refLength == altLength: # MNV
            #inversion
            if self.__ref == reverse(self.__alt):
                attributes.update({
                    'variant_class': "inversion",
                    'variant_class_abbrev': "MNV",
                    'display_allele': 'inv' + self.__ref,
                    'sequence_allele': truncate_allele(self.__ref) + '/' + truncate_allele(self.__alt),
                    'location_end': endLocation
                })
            else:
                attributes.update({
                    'variant_class': "substitution",
                    'variant_class_abbrev': "MNV",
                    'display_allele': normRef + ">" + normAlt,
                    'sequence_allele': truncate_allele(normRef) + '/' + truncate_allele(normAlt),
                    'location_start': position,
                    'location_end': endLocation
                })
        # end MNV

        elif refLength > altLength: # deletions
            attributes.update({'location_start': position})

            if nAltLength > 1: # INDEL
                attributes.update({
                    'variant_class': 'indel',
                    'variant_class_abbrev': 'INDEL'
                    })

                if nRefLength == 0:
                    displayRef = self.__ref[1:] # strip first character from reference
                    if displayRef == normAlt: # duplication
                        attributes.update({
                            'location_start': position + 1,
                            'location_end': position + 1,
                            'display_allele': 'dup' + normAlt,
                            'sequence_allele': 'dup' + truncate_allele(normAlt),
                            'variant_class': 'duplication',
                            'variant_class_abbrev': 'DUP'
                        })
                    else:
                        attributes.update({
                            'location_end': endLocation,
                            'location_start': position + 1,
                            'display_allele': "del" + displayRef + "ins" + normAlt,
                            'sequence_allele': truncate_allele(displayRef) + "/" + truncate_allele(normAlt)
                        })
                else:
                    attributes.update({
                        'location_end': endLocation,
                        'location_start': position + 1,
                        'display_allele': "del" + normRef + "ins" + normAlt,
                        'sequence_allele': truncate_allele(normRef) + "/" + truncate_allele(normAlt)
                        })
                # end INDEL

            else: # deletion
                attributes.update({
                    'variant_class': "deletion",
                    'variant_class_abbrev': "DEL",
                    'location_end': endLocation,
                    'display_allele': "del" + normRef,
                    'sequence_allele': truncate_allele(normRef) + "/-"
                    })

        # end deletions

        elif refLength < altLength: # insertions
            attributes.update({'location_start': position + 1})

            if refLength > 1: # INDEL
                attributes.update({
                    'variant_class': 'indel',
                    'variant_class_abbrev': 'INDEL'
                    })

                if nRefLength == 0:
                    displayRef = self.__ref[1:] # strip first character from reference
                    if displayRef == normAlt: # duplication
                        attributes.update({
                            'location_end': position + 1,
                            'display_allele': 'dup' + normAlt,
                            'sequence_allele': 'dup' + truncate_allele(normAlt),
                            'variant_class': 'duplication',
                            'variant_class_abbrev': 'DUP'
                            })
                    else:
                        attributes.update({
                            'location_end': endLocation,
                            'display_allele': "del" + displayRef + "ins" + normAlt,
                            'sequence_allele': truncate_allele(displayRef) + "/" + truncate_allele(normAlt)
                        })
                else:
                    attributes.update({
                        'location_end': endLocation,
                        'display_allele': "del" + normRef + "ins" + normAlt,
                        'sequence_allele': truncate_allele(normRef) + "/" + truncate_allele(normAlt)
                        })
                # end INDEL
            else: # insertion
                attributes.update({
                    'variant_class': "insertion",
                    'variant_class_abbreve': "INS",
                    'location_end': position + 1,
                    'display_allele': "ins" + normAlt,
                    'sequence_allele': "-/" + truncate_allele(normAlt)
                })

        return attributes
