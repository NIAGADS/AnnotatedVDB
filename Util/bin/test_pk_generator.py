#!/usr/bin/env python
''' test primary key generator '''
import argparse
import requests

from os import path

from AnnotatedVDB.Util.primary_key_generator import VariantPKGenerator
from GenomicsDBData.Util.utils import truncate, warning, xstr, die

from ga4gh.core import sha512t24u, ga4gh_digest, ga4gh_serialize, ga4gh_identify, is_pjs_instance
from ga4gh.vrs import models, vrs_deref, vrs_enref

NCBI_SERVICES_URL="https://api.ncbi.nlm.nih.gov/variation/v0/"
GA4GHR_VR_ENDPOINT="ga4gh-vr/"
REFSQ_CHR1="refseq:NC_000001.11"

GENOME_BUILD="GRCh38"




def validate(metaseqId, generator):
    ''' validate against NCBI '''
    ncbiGA4GR_VR_endpoint = path.join(NCBI_SERVICES_URL, GA4GHR_VR_ENDPOINT)
    headers = {'Content-type': 'application/vnd.ga4gh.vr+json'}    

    chrm, position = metaseqId.split(':')[0:2]
    sequenceId = 'GRCh38:' + xstr(chrm)

    alleleDict = generator.get_vrs_allele_dict(metaseqId)
    # gnomad = generator.translate_vrs(alleleDict, 'gnomad')
    # warning("reverse translation", gnomad.replace('-', ':'))

    spdi = generator.translate_vrs(alleleDict, 'spdi')
    warning("INFO", "SPDI - ", spdi[0])
    

        
def run():
    ''' run test '''

    variants = {
        "1:110852777:CCTGCTCCTT:CACCCCCACAGCTGTTACCCAGCGCCACACACAGAGCAGACGCTGAATCACTGCTTATTGACTGAATCAGCAATGGGGTACCTGCTCCTG" : "rs71575164",
        "1:110241183:CCTTTCCCGCTTCTCTGTCCTGCAGCCAGCTGATCGTGGGACTTCACTCCAAAATTATGTGAGCCAGTTCCCACAAGAGATAC:CTTTCCTGCTTCTCTGTCCTGCAGCCAGCTGATTGTGGGACTTCACTCCAAAATTGTGTGAGCCAATTCCCATAAGAGATAA" : "rs386634503",
        "13:32936731:G:C": None,
        "1:148893911:TGGCCAACA:TAGCCAACG": "rs71261250",
        "1:203342335:T:TTTCCTTCCTTTCTTCCTTCCTTTCTCCTTCCTTCCTTTCTCCCTTCCTTCCTTCC": "rs71142598",
        "1:178495286:CTGGCTCTGGAATCTAGAGTTTAGATTAATTTATTAATGGGTAAAAAGACAAATGAGCTGGAACAAAAAAAGTCCGTGTGACTCCAAAACCCAACCCCTGGTCTTGTTCTTCATGTTCAGGAACACTGACTTCACATTTCTTTCTTAACTGAGTGTTCATAATCATAACTTGATTGATGATCAAATGTCTCTCCCACCTCCTTTCTCATCAGTATTCATATGGCAGCTTGTTATTCCTCCCAGGTAAACTGGTTGAGAGGGAATATAAATTTGTTTAAGAGGTACCAATGTCAACTTTGCTCTAAAGCATCATGAACTTGGCCAGCTGCGGTGGCTCATGACTGTAATCCCAGCACTTTGGGAGGCCAAGGTGGGTGGATCACCTGAGGTCAGGAGTTTGAGACCAGCCTGACCAACATAGAGAAACCCCA:CACTTTGGGAGGCCGAGACGGGCAGATCACGAGGTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCCG": None
        }

    
    generator = VariantPKGenerator(GENOME_BUILD, args.proxyPath, debug=args.debug)
    ngenerator = VariantPKGenerator(GENOME_BUILD, args.proxyPath, debug=args.debug, normalize=True) 
    for metaseqId, refSnp in variants.items():
        warning("\n")
        warning("# ============= AS IS =========== ")
        pk = generator.generate_primary_key(metaseqId, refSnp)
        warning("INFO", "PK", pk)
        validate(metaseqId, generator)
        warning("# ============= NORMALIZED =========== ")
        pk = ngenerator.generate_primary_key(metaseqId, refSnp)
        warning("INFO", "PK", pk)
        validate(metaseqId, ngenerator)

  
    
    
    
     


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--proxyPath', help="path to seqrepo", required=True)
    parser.add_argument('--debug', action='store_true')
    args = parser.parse_args()



    run()
