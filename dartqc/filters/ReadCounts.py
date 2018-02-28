import numpy

from Dataset import Dataset
from PipelineOptions import Filter
from FilterResult import FilterResult


class ReadCountsFilter(Filter):
    def get_order(self) -> int:
        return 0

    def get_cmd_type(self):
        return lambda s: [float(item.strip()) if len(item.strip()) > 0 else None for item in s[1:-1].split(',')]

    def get_name(self) -> str:
        return "read_counts"

    def get_cmd_help(self) -> str:
        return "Silence call if read count is < given value"

    def filter(self, dataset: Dataset, threshold: float, unknown_args: [], **kwargs) -> FilterResult:
        silenced = FilterResult()
        filtered_calls = dataset.get_filtered_calls()

        for snp, snp_counts in dataset.read_counts.items():
            # Identify which samples are silenced for this SNP
            filter_vector = [False if sum(sample_read_count) <= threshold else True for sample_read_count in snp_counts]

            # Silence the calls
            for idx, call in enumerate(filtered_calls[snp]):
                if not filter_vector[idx] and filtered_calls[snp][idx] != Dataset.missing:
                    silenced.silenced_call(snp, dataset.samples[idx].id)

        return silenced

ReadCountsFilter()
