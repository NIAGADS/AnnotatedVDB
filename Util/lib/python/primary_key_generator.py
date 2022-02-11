""" 
generator for variant record primary key

rules:

for SNV's and short INDELS (total length alleles <= 50)
======================================================
1. use the NCBI dbSNP SPDI - https://www.ncbi.nlm.nih.gov/variation/notation/
or Sequence:Position:Deletion:Insertion 
 --> w/restriction that deletion/reference length should NOT be substituted for the deletion/reference sequence b/c we cannot guarantee that the deletion sequence comes from the reference genome

2. if the variant can be associated with a NCBI refSNP, append the rsId,
so that the full PK is:

SPDI:refSnpId

TODO: 3. if the variant can be associated with an EVA short snp identifier (and not a refSnpId), append the ssId, so that the full PK is:

SPDI:ssId


for large INDELS/SVS (total length alleles > 50)
======================================================
S:P:VRS Computed Identifier:externalID

see: https://vrs.ga4gh.org/en/stable/impl-guide/computed_identifiers.html#computed-identifiers

e.g. ga4gh:VA.EgHPXXhULTwoP4-ACfs-YCXaeUQJBjH_

for this use case will strip the prefix code 'ga4gh:VA.' as it will be identical for each serialization & the . and : will interfer w/parsing / i.e. only retain the digest

e.g., 
1:110852777:CCTGCTCCTT:CACCCCCACAGCTGTTACCCAGCGCCACACACAGAGCAGACGCTGAATCACTGCTTATTGACTGAATCAGCAATGGGGTACCTGCTCCTG:rs71575164
-->

e.g., 
1:110241183:CCTTTCCCGCTTCTCTGTCCTGCAGCCAGCTGATCGTGGGACTTCACTCCAAAATTATGTGAGCCAGTTCCCACAAGAGATAC:CTTTCCTGCTTCTCTGTCCTGCAGCCAGCTGATTGTGGGACTTCACTCCAAAATTGTGTGAGCCAATTCCCATAAGAGATAA:rs386634503
-->



"""
import json # for debug
from ga4gh.core import ga4gh_identify, ga4gh_serialize
from ga4gh.vrs.extras.translator import Translator
from ga4gh.vrs.dataproxy import create_dataproxy
from GenomicsDBData.Util.utils import warning, die, xstr


class VariantPKGenerator(object):
    ''' generator for variant record primary key '''
    

    def __init__(self, genomeBuild, seqrepoProxyPath, maxSequenceLength=50, normalize=False, verbose=True, debug=False):
        self._verbose = verbose
        self._debug = debug
        self._genomeBuild = genomeBuild
        self._maxSequenceLength = maxSequenceLength
        self._ga4gh_sequence_map = {}
        self._translator = None
        self._set_seqrepo_translator(seqrepoProxyPath, normalize)


    def _set_seqrepo_translator(self, proxyPath, normalize=False):
        ''' set seqrepo translator '''
        dataProxy = create_dataproxy("seqrepo+file://" + proxyPath)
        self._translator = Translator(data_proxy=dataProxy)
        self._translator.normalize = normalize



    def translate_vrs(self, vrsDict, formatSpec="spdi"):
        ''' translate vrsDict to human readable format '''
        if formatSpec not in ['hgvs', 'spdi']:
            raise IndexError("invalid format spec: spd, hgvs")

        return self._translator.translate_to(vrsDict, formatSpec)

        
    def generate_primary_key(self, metaseqId, externalId=None):
        ''' generate the primary key from the metaseqId & externalId'''
        chrm, position, ref, alt = metaseqId.split(':')

        pk = [chrm, position]
        if len(ref) + len(alt) <= self._maxSequenceLength:
            pk.extend([ref,alt])
        else:
            pk.append(self.get_vrs_identifier(metaseqId))

        if externalId is not None:
          pk.append(externalId)

        return ':'.join(pk)
        

    def _get_ga4gh_vr_sequence_id(self, chrm):
        ''' fetch ga4gh-vr sequence id for the chromosome '''
        if chrm not in self._ga4gh_sequence_map:
            gvrId = translate_sequence_identifier(self._genomeBuild + ":" + xstr(chrm), "ga4gh")
            self._ga4gh_sequence_map[chrm] = gvrId[0]

        return self._ga4gh_sequence_map[chrm]



    def get_vrs_allele_dict(self, metaseqId, serialize=False, toJson=False):
        ''' get a ga4gh VRS allele dict for the specified variant '''
        # transform in to gnomad id (replace ':' with '-')
        gnomadExpr = metaseqId.replace(':', '-')
        alleleDict = self._translator._from_gnomad(gnomadExpr, self._genomeBuild)

        
        if serialize:
            return ga4gh_serialize(alleleDict)
        if toJson:
            return alleleDict.for_json() # json.loads(alleleDict.decode("utf-8"))

        return alleleDict
        
        

    def get_vrs_identifier(self, metaseqId):
        ''' return computed ga4gh identifier for the allelic sequence '''

        alleleDict = self.get_vrs_allele_dict(metaseqId)
        
        if self._debug:
            debugOutput = {
                "Input Variant":  metaseqId,
                "VRS Representation" : alleleDict.for_json(),
                }
            warning(json.dumps(debugOutput, indent=4, sort_keys=True))

        vrsComputedId = ga4gh_identify(alleleDict)
        
      
        return vrsComputedId


        
        

        
