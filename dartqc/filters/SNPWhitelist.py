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

        snp_list = [snp_def.allele_id for snp_def in dataset.snps if snp_def.allele_id not in dataset.filtered.snps]

        # Check if the whitelist is a file instead of a list of items.
        if len(whitelist) == 1 and ("/" in whitelist[0] or "\\" in whitelist[0]):
            file_path = whitelist[0]

            with open(file_path, "r") as file:
                whitelist = list(set(value.strip() for line in file.readlines() for value in line.split(",") if value.strip() in snp_list))
        else:
            whitelist = [snp for snp in whitelist if snp in snp_list]

        silenced.snps = [snp for snp in snp_list if snp not in whitelist]
        return silenced

SNPWhitelistFilter()
