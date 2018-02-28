import logging
import pandas

from Dataset import Dataset
from FilterResult import FilterResult
from PipelineOptions import Filter

log = logging.getLogger(__file__)

class MinSampleDataFilter(Filter):
    def get_name(self):
        return "min_sample"

    def get_order(self) -> int:
        return 10

    def get_cmd_type(self):
        return lambda s: [float(item.strip()) if len(item.strip()) > 0 else None for item in s[1:-1].split(',')]

    def get_cmd_help(self) -> str:
        return "filter samples > missingness per sample"

    def filter(self, dataset: Dataset, threshold: float, unknown_args: [], **kwargs) -> FilterResult:
        silenced = FilterResult()

        filtered_calls = dataset.get_filtered_calls()

        # Convert calls from tuple to 0,1,2,- for pandas comparison/sum to work
        simple_calls = Dataset.encode_single(filtered_calls)

        df = pandas.DataFrame(simple_calls)
        mind = (df == dataset.encoded_missing).sum(axis=1)
        mind /= len(dataset.snps)  # Series

        for idx, min_data in enumerate(mind):
            if (1 - min_data) < threshold:
                silenced.silenced_sample(dataset.samples[idx].id)

        return silenced


MinSampleDataFilter()
