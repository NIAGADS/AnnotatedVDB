#!/usr/bin/env python3
''' test variant validator '''

from AnnotatedVDB.Util.database.variant import VariantRecord

variants = {
    "1:110852777:CCTGCTCCTT:CACCCCCACAGCTGTTACCCAGCGCCACACACAGAGCAGACGCTGAATCACTGCTTATTGACTGAATCAGCAATGGGGTACCTGCTCCTG" : "rs71575164",
    "1:110241183:CCTTTCCCGCTTCTCTGTCCTGCAGCCAGCTGATCGTGGGACTTCACTCCAAAATTATGTGAGCCAGTTCCCACAAGAGATAC:CTTTCCTGCTTCTCTGTCCTGCAGCCAGCTGATTGTGGGACTTCACTCCAAAATTGTGTGAGCCAATTCCCATAAGAGATAA" : "rs386634503",
    "13:32936731:G:C": None,
    "M:11257:C:T": "rs377469212",
    "1:148893911:TGGCCAACA:TAGCCAACG": "rs71261250",
    "1:203342335:T:TTTCCTTCCTTTCTTCCTTCCTTTCTCCTTCCTTCCTTTCTCCCTTCCTTCCTTCC": "rs71142598",
    "1:178495286:CTGGCTCTGGAATCTAGAGTTTAGATTAATTTATTAATGGGTAAAAAGACAAATGAGCTGGAACAAAAAAAGTCCGTGTGACTCCAAAACCCAACCCCTGGTCTTGTTCTTCATGTTCAGGAACACTGACTTCACATTTCTTTCTTAACTGAGTGTTCATAATCATAACTTGATTGATGATCAAATGTCTCTCCCACCTCCTTTCTCATCAGTATTCATATGGCAGCTTGTTATTCCTCCCAGGTAAACTGGTTGAGAGGGAATATAAATTTGTTTAAGAGGTACCAATGTCAACTTTGCTCTAAAGCATCATGAACTTGGCCAGCTGCGGTGGCTCATGACTGTAATCCCAGCACTTTGGGAGGCCAAGGTGGGTGGATCACCTGAGGTCAGGAGTTTGAGACCAGCCTGACCAACATAGAGAAACCCCA:CACTTTGGGAGGCCGAGACGGGCAGATCACGAGGTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCCG": None
    }


validator = VariantRecord(gusConfigFile=None, verbose=False, debug=True)

for metaseqId, refSnpId in variants.items():
    
    if refSnpId is not None:
        print(refSnpId)
        print(validator.exists(refSnpId, 'REFSNP', chromosome=None))