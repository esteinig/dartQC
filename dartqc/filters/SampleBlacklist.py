import re

from dartqc.FilterResult import FilterResult
from dartqc.Dataset import Dataset
from dartqc.PipelineOptions import Filter


class MinSNPDataFilter(Filter):
    def get_order(self) -> int:
        return 0

    def get_cmd_type(self):
        return Filter.LIST_OF_LISTS

    def get_name(self) -> str:
        return "sample_blacklist"

    def get_cmd_help(self) -> str:
        return "Remove samples before filtering - Pattern: [[sample_id, sample_id, ...], ...]"

    def filter(self, dataset: Dataset, threshold: [str], unknown_args: [], **kwargs) -> FilterResult:
        silenced = FilterResult()

        for sample_id in threshold:
            silenced.silenced_sample(sample_id)

        return silenced

MinSNPDataFilter()
