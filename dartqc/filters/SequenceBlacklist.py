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

        sequence_list = {snp_def.allele_id : snp_def.sequence_ref for snp_def in dataset.snps if snp_def.allele_id not in dataset.filtered.snps}

        # Check if the blacklist is a file instead of a list of items.
        if len(blacklist) == 1 and ("/" in blacklist[0] or "\\" in blacklist[0]):
            file_path = blacklist[0]

            with open(file_path, "r") as file:
                blacklist = list(set(value.strip() for line in file.readlines() for value in line.split(",") if value.strip() in sequence_list.values()))
        else:
            blacklist = [ref_seq for ref_seq in blacklist if ref_seq in sequence_list.values()]

        silenced.snps = [snp for (snp, seq) in sequence_list.items() if seq in blacklist]
        return silenced

SequenceWhitelistFilter()
