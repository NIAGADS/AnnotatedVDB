'''
Enum for handling sets of VEP consequences &  facilitate ADSP consequence ranking
'''
# pylint: disable=line-too-long,invalid-name,no-self-use
from enum import Enum
import GenomicsDBData.Util.list_utils as lu

class ConseqGroup(Enum):
    __order__ = 'HIGH_IMPACT NMD NON_CODING_TRANSCRIPT LOW_IMPACT'
    HIGH_IMPACT = ['transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant',
                   'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost',
                   'inframe_insertion', 'inframe_deletion', 'missense_variant',
                   'protein_altering_variant', 'splice_region_variant',
                   'incomplete_terminal_codon_variant', 'stop_retained_variant',
                   'start_retained_variant', 'synonymous_variant',
                   'coding_sequence_variant', '5_prime_UTR_variant', '3_prime_UTR_variant',
                   'intron_variant'] # GRP 4

    NMD = ['NMD_transcript_variant'] # GRP 1

    NON_CODING_TRANSCRIPT = ['non_coding_transcript_exon_variant', 'non_coding_transcript_variant'] # GRP 2

    LOW_IMPACT = ['mature_miRNA_variant', 'non_coding_transcript_variant',
                  'non_coding_transcript_exon_variant', 'upstream_gene_variant',
                  'downstream_gene_variant', 'TF_binding_site_variant', 'TFBS_ablation',
                  'TF_binding_site_variant', 'regulatory_region_variant', 'intergenic_variant'] # GRP 3


    @classmethod
    def get_all_terms(cls):
        ''' return all terms
        usage: ConseqGroups.get_all_terms '''
        terms = []
        for x in ConseqGroup:
            if x.name != 'NON_CODING_TRANSCRIPT': # subset of LOW_IMPACT & want to preserve order
                terms += x.value

        return terms


    @classmethod
    def get_complete_indexed_dict(cls):
        ''' returns indexed OrderedDict containing all ConseqGroup terms
        usage: ConseqGroups.get_complete_indexed_dict()
        '''
        values = ConseqGroup.get_all_terms()
        return lu.list_to_indexed_dict(values)


    @classmethod
    def validate_terms(cls, conseqs):
        ''' verify that all terms fall into a conseq group
        raise error if novel term is found so that code can be updated to assign to a specific
        consequence group

        see https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html
        for possible updates

        usage: ConseqGroup.validate_terms(conseqs)
        '''
        validTerms = ConseqGroup.get_all_terms()

        for c in conseqs:
            if not lu.is_subset(c.split(','), validTerms):
                for term in c.split(','): # determine which specific term failed
                    if term not in validTerms:
                        raise IndexError("Consequence combination `" + c \
                                             + "` contains an invalid consequence: `" + term \
                                             + "`. Please update the `ConseqGroup` enum " \
                                             + "(conseq_group_enum.py) after reviewing " \
                                             + "https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html")

        return True


    def __str__(self):
        ''' usage: str(ConseqGroup.HIGH_IMPACT) '''
        return ','.join(self.value)


    def get_group_members(self, conseqs, requireSubset=True):
        ''' given a list of combinations, extracts those that
          include group members following ADSP Annotation rules

         usage: ConseqGroup.LOW_IMPACT.get_members(conseqs)
        '''

        ConseqGroup.validate_terms(conseqs)

        if requireSubset: # all terms must be in this group for a consequence to be a member
            return [x for x in conseqs if lu.is_subset(x.split(','), self.value)]

        else:
            if self.name == 'HIGH_IMPACT': # need to ignore consequence combos that include NMD/non-coding
                return [x for x in conseqs if lu.is_overlapping_list(x.split(','), self.value) \
                            and not lu.is_overlapping_list(x.split(','), ConseqGroup.NON_CODING_TRANSCRIPT.value) \
                            and not lu.is_overlapping_list(x.split(','), ConseqGroup.NMD.value)]
            else:
                return [x for x in conseqs if lu.is_overlapping_list(x.split(','), self.value)]


    def toDict(self):
        ''' returns indexed OrderedDict containing ConseqGroup terms
        usage: ConseqGroups.LOW_IMPACT.toDict()
        '''
        return lu.list_to_indexed_dict(self.value)
