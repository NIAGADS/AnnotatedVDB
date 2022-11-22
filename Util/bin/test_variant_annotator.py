#!/usr/bin/env python3
from AnnotatedVDB.Util.variant_annotator import VariantAnnotator, reverse_complement
from GenomicsDBData.Util.utils import warning, print_dict

variants = [
    '22:11212877:TAAAATATCAAAGTACACCAAATACATATTATATACTGTACAC:T',
    '22:11212877:TAAAATATCAAAGTACACCAAATACATATTATATACTGTACAC:TAAAATATCAAAGTACACCAAATACATATTATATACTGTACACAAAATATCAAAGTACACCAAATACATATTATATACTGTACAC'
]

for v in variants:
    chrm, pos, ref, alt = v.split(':')
    annotator = VariantAnnotator(ref, alt, chrm, int(pos))
    warning("Normalized Alleles:", annotator.get_normalized_alleles(True) )
    warning("Display Attributes:", print_dict(annotator.get_display_attributes(), pretty=True))
    warning("Reverse Comp:", reverse_complement(alt))
