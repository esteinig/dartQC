import re

from dartqc.FilterResult import FilterResult
from dartqc.Dataset import Dataset
from dartqc.PipelineOptions import Filter


class SNPWhitelistFilter(Filter):
    def get_order(self) -> int:
        return 1

    def get_cmd_type(self):
        return Filter.LIST_OF_LISTS

    def get_name(self) -> str:
        return "SNP Whitelist"

    def get_cmd_names(self) -> [str]:
        return ["--snp_whitelist"]

    def get_cmd_help(self) -> str:
        return "Remove SNPs that don't match - Pattern: [allele_id, allele_id, ...]"

    def get_description(self) -> str:
        return "Remove any SNPs not in this list"

    def filter(self, dataset: Dataset, whitelist: [str], unknown_args: [], **kwargs) -> FilterResult:
        silenced = FilterResult()

        for snp_def in dataset.snps:
            if snp_def.allele_id not in whitelist and snp_def.allele_id not in dataset.filtered.snps:
                silenced.silenced_snp(snp_def.allele_id)

        return silenced

SNPWhitelistFilter()
