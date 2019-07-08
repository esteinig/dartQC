import re

from dartqc.FilterResult import FilterResult
from dartqc.Dataset import Dataset
from dartqc.PipelineOptions import Filter


class SequenceWhitelistFilter(Filter):
    def get_order(self) -> int:
        return 1

    def get_cmd_type(self):
        return Filter.LIST_OF_LISTS

    def get_name(self) -> str:
        return "Sequence Blacklist"

    def get_cmd_names(self) -> [str]:
        return ["--sequence_blacklist"]

    def get_cmd_help(self) -> str:
        return "Remove all SNPs that have this sequence (Note: this is sequence, not ref_sequence) - Pattern: [sequence, sequence, ...]"

    def get_description(self) -> str:
        return "Remove all SNPs that have this sequence (Note: this is sequence, not ref_sequence)"

    def filter(self, dataset: Dataset, blacklist: [str], unknown_args: [], **kwargs) -> FilterResult:
        silenced = FilterResult()

        for snp_def in dataset.snps:
            if snp_def.sequence in blacklist and snp_def.allele_id not in dataset.filtered.snps:
                silenced.silenced_snp(snp_def.allele_id)

        return silenced

SequenceWhitelistFilter()
