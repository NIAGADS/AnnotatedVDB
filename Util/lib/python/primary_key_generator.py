""" 
generator for variant record primary key

rules:

for SNV's and short INDELS (<=20 bp per allele):
=================================================
1. use the NCBI dbSNP SPDI - https://www.ncbi.nlm.nih.gov/variation/notation/
or Sequence:Position:Deletion:Insertion 
 --> w/restriction that deletion/reference length should NOT be substituted for the deletion/reference sequence b/c we cannot guarantee that the deletion sequence comes from the reference genome

2. if the variant can be associated with a NCBI refSNP, append the rsId,
so that the full PK is:

SPDI_refSnpId

TODO: 3. if the variant can be associated with an EVA short snp identifier (and not a refSnpId), append the ssId, so that the full PK is:

SPDI_ssId


for large INDELS/SVS (any single allele is >20 bp):
=================================================
keep SP and concatenated external database ID (refSNP or ss)
use GA4GH sequence coding algorithm to generate unique computed sequence representation for any allelic sequence >20bp

see: https://vrs.ga4gh.org/en/stable/impl-guide/computed_identifiers.html#computed-identifiers

e.g., for the following variant, the insertion/alt sequence will be replaced with the unique computed sequence identifier

1:110852777:CCTGCTCCTT:CACCCCCACAGCTGTTACCCAGCGCCACACACAGAGCAGACGCTGAATCACTGCTTATTGACTGAATCAGCAATGGGGTACCTGCTCCTG_rs71575164
-->

e.g., for the following variant, both allelic sequences will be replaced with the computed identifier:
1:110241183:CCTTTCCCGCTTCTCTGTCCTGCAGCCAGCTGATCGTGGGACTTCACTCCAAAATTATGTGAGCCAGTTCCCACAAGAGATAC:CTTTCCTGCTTCTCTGTCCTGCAGCCAGCTGATTGTGGGACTTCACTCCAAAATTGTGTGAGCCAATTCCCATAAGAGATAA_rs386634503
-->

"""



class VariantPKGenerator(object):
    ''' generator for variant record primary key '''
    

    def __init__(self, rankingFileName, verbose=True):
        self._verbose = verbose
