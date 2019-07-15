import re

from dartqc.FilterResult import FilterResult
from dartqc.Dataset import Dataset
from dartqc.PipelineOptions import Filter


class SampleBlacklistFilter(Filter):
    def get_order(self) -> int:
        return 0

    def get_cmd_type(self):
        return Filter.LIST_OF_LISTS

    def get_name(self) -> str:
        return "Sample Blacklist"

    def get_cmd_names(self) -> [str]:
        return ["--sample_blacklist"]

    def get_cmd_help(self) -> str:
        return "Remove samples - Pattern: [sample_id, sample_id, ...]"

    def get_description(self) -> str:
        return "Remove list of samples"

    def filter(self, dataset: Dataset, blacklist: [str], unknown_args: [], **kwargs) -> FilterResult:
        silenced = FilterResult()

        sample_list = [sample_def.id for sample_def in dataset.samples if sample_def.id not in dataset.filtered.samples]

        # Check if the blacklist is a file instead of a list of items.
        if len(blacklist) == 1 and ("/" in blacklist[0] or "\\" in blacklist[0]):
            file_path = blacklist[0]

            with open(file_path, "r") as file:
                blacklist = list(set(value.strip() for line in file.readlines() for value in line.split(",") if value.strip() in sample_list))
        else:
            blacklist = [sample for sample in blacklist if sample in sample_list]

        silenced.samples = blacklist
        return silenced

SampleBlacklistFilter()
