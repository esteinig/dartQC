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
        return "Sequence Whitelist"

    def get_cmd_names(self) -> [str]:
        return ["--sequence_whitelist"]

    def get_cmd_help(self) -> str:
        return "Remove SNPs if their sequence isn't in this list (Note: this is sequence, not ref_sequence) - Pattern: [sequence, sequence, ...]"

    def get_description(self) -> str:
        return "Remove SNPs if their sequence isn't in this list (Note: this is sequence, not ref_sequence)"

    def filter(self, dataset: Dataset, whitelist: [str], unknown_args: [], **kwargs) -> FilterResult:
        silenced = FilterResult()

        sequence_list = {snp_def.allele_id : snp_def.sequence_ref for snp_def in dataset.snps if snp_def.allele_id not in dataset.filtered.snps}

        # Check if the whitelist is a file instead of a list of items.
        if len(whitelist) == 1 and ("/" in whitelist[0] or "\\" in whitelist[0]):
            file_path = whitelist[0]

            with open(file_path, "r") as file:
                whitelist = list(set(value.strip() for line in file.readlines() for value in line.split(",") if value.strip() in sequence_list.values()))
        else:
            whitelist = [ref_seq for ref_seq in whitelist if ref_seq in sequence_list.values()]

        silenced.snps = [snp.allele_id for snp, seq in sequence_list.items() if seq not in whitelist]
        return silenced

SequenceWhitelistFilter()
