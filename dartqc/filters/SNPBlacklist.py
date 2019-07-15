import re

from dartqc.FilterResult import FilterResult
from dartqc.Dataset import Dataset
from dartqc.PipelineOptions import Filter


class SNPBlacklistFilter(Filter):
    def get_order(self) -> int:
        return 1

    def get_cmd_type(self):
        return Filter.LIST_OF_LISTS

    def get_name(self) -> str:
        return "SNP Blacklist"

    def get_cmd_names(self) -> [str]:
        return ["--snp_blacklist"]

    def get_cmd_help(self) -> str:
        return "Remove SNPs - Pattern: [allele_id, allele_id, ...]"

    def get_description(self) -> str:
        return "Remove list of SNPs"

    def filter(self, dataset: Dataset, blacklist: [str], unknown_args: [], **kwargs) -> FilterResult:
        silenced = FilterResult()

        snp_list = [snp_def.allele_id for snp_def in dataset.snps if snp_def.allele_id not in dataset.filtered.snps]

        # Check if the blacklist is a file instead of a list of items.
        if len(blacklist) == 1 and ("/" in blacklist[0] or "\\" in blacklist[0]):
            file_path = blacklist[0]

            with open(file_path, "r") as file:
                blacklist = list(set(value.strip() for line in file.readlines() for value in line.split(",") if value.strip() in snp_list))
        else:
            blacklist = [snp for snp in blacklist if snp in snp_list]

        silenced.snps = blacklist
        return silenced

SNPBlacklistFilter()
