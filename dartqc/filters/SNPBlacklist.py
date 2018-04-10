import re

from dartqc.FilterResult import FilterResult
from dartqc.Dataset import Dataset
from dartqc.PipelineOptions import Filter


class MinSNPDataFilter(Filter):
    def get_order(self) -> int:
        return 1

    def get_cmd_type(self):
        return Filter.LIST_OF_LISTS

    def get_name(self) -> str:
        return "SNP Blacklist"

    def get_cmd_names(self) -> [str]:
        return ["--snp_blacklist"]

    def get_cmd_help(self) -> str:
        return "Remove SNPs - Pattern: [[allele_id, allele_id, ...], ...]"

    def get_description(self) -> str:
        return "Remove list of SNPs"

    def filter(self, dataset: Dataset, threshold: [str], unknown_args: [], **kwargs) -> FilterResult:
        silenced = FilterResult()

        for allele_id in threshold:
            silenced.silenced_snp(allele_id)

        return silenced

MinSNPDataFilter()
