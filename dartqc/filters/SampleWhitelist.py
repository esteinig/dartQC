import re

from dartqc.FilterResult import FilterResult
from dartqc.Dataset import Dataset
from dartqc.PipelineOptions import Filter


class SampleWhitelistFilter(Filter):
    def get_order(self) -> int:
        return 0

    def get_cmd_type(self):
        return Filter.LIST_OF_LISTS

    def get_name(self) -> str:
        return "Sample Whitelist"

    def get_cmd_names(self) -> [str]:
        return ["--sample_whitelist"]

    def get_cmd_help(self) -> str:
        return "Remove samples not in this list - Pattern: [sample_id, sample_id, ...]"

    def get_description(self) -> str:
        return "Remove samples not in this list"

    def filter(self, dataset: Dataset, whitelist: [str], unknown_args: [], **kwargs) -> FilterResult:
        silenced = FilterResult()

        sample_list = [sample_def.id for sample_def in dataset.samples if sample_def.id not in dataset.filtered.samples]

        # Check if the whitelist is a file instead of a list of items.
        if len(whitelist) == 1 and ("/" in whitelist[0] or "\\" in whitelist[0]):
            file_path = whitelist[0]

            with open(file_path, "r") as file:
                whitelist = list(set(value.strip() for line in file.readlines() for value in line.split(",") if value.strip() in sample_list))
        else:
            whitelist = [sample for sample in whitelist if sample in sample_list]

        silenced.samples = [sample for sample in sample_list if sample not in whitelist]
        return silenced

SampleWhitelistFilter()
