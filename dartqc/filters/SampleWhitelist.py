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

        for sample_def in dataset.samples:
            if sample_def.id not in whitelist and sample_def.id not in dataset.filtered.samples:
                silenced.silenced_sample(re.sub(r"['\"]", "", sample_def.id))

        return silenced

SampleWhitelistFilter()
