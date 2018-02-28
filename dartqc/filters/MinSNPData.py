from FilterResult import FilterResult
from Dataset import Dataset
from PipelineOptions import Filter


class MinSNPDataFilter(Filter):
    def get_order(self) -> int:
        return 11

    def get_cmd_type(self):
        return lambda s: [float(item.strip()) if len(item.strip()) > 0 else None for item in s[1:-1].split(',')]

    def get_name(self) -> str:
        return "min_snp"

    def get_cmd_help(self) -> str:
        return "filter snps <= call rate of snp"

    def filter(self, dataset: Dataset, threshold: float, unknown_args: [], **kwargs) -> FilterResult:
        silenced = FilterResult()

        filtered_calls = dataset.get_filtered_calls()

        for allele_id, snp_calls in filtered_calls.items():
            perc_data = 1 - (snp_calls.count(dataset.missing) / len(dataset.samples))
            if perc_data < threshold:
                silenced.silenced_snp(allele_id)

        return silenced

MinSNPDataFilter()
