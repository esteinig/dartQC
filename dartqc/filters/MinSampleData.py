import logging

import numpy
import pandas
import re

from dartqc.Dataset import Dataset
from dartqc.FilterResult import FilterResult
from dartqc.PipelineOptions import Filter

log = logging.getLogger(__file__)


class MinSampleDataFilter(Filter):
    def get_name(self):
        return "Min Sample Data"

    def get_cmd_names(self) -> [str]:
        return ["--min_sample","--mind"]

    def get_order(self) -> int:
        return 10

    def get_cmd_type(self):
        return Filter.LIST_OF_FLOAT

    def get_cmd_help(self) -> str:
        return "filter samples > missingness per sample. Pattern: [.7,.8,.9]"

    def get_description(self):
        return "filter samples > missingness per sample (0 to 1 decimal -> bigger requires more data)"

    def filter(self, dataset: Dataset, threshold: float, unknown_args: [], **kwargs) -> FilterResult:
        silenced = FilterResult()

        filtered_calls, filtered_snps, filtered_samples = dataset.get_filtered_calls()

        # Get calls as matrix - swap SNP/sample axes & trim back to only first allele's call (much quicker processing)
        numpy_matrix = numpy.asarray([filtered_calls[snp.allele_id] for snp in filtered_snps])
        numpy_matrix = numpy.stack(numpy_matrix, axis=2)  # [SNPs][samples][calls] -> [SNPs][calls][samples]
        numpy_matrix = numpy.stack(numpy_matrix, axis=1)  # [SNPs][calls][samples] -> [calls][samples][SNPs]
        numpy_matrix = numpy_matrix[0]  # Only get first allele calls (if "-" -> missing)

        # ignored_samples = numpy.asarray([True if filtered_samples[idx] in dataset.filtered.samples else False
        #                                  for idx in range(len(dataset.samples))])

        for idx, sample_calls in enumerate(numpy_matrix):
            # if ignored_samples[idx]:
            #     continue

            # first_allele_calls = numpy.dstack(sample_calls)[0][0]  # Get just the first allele (if "-" => missing)
            missing_idxs = numpy.where(sample_calls == "-")[0]
            missing = len(missing_idxs) / len(filtered_snps)

            if (1 - missing) <= threshold:
                silenced.silenced_sample(filtered_samples[idx].id)

            if idx % 1000 == 0:
                log.debug("Completed {} of {}".format(idx, len(numpy_matrix)))


        # # Convert calls from tuple to 0,1,2,- for pandas comparison/sum to work
        # simple_calls = Dataset.encode_single(filtered_calls)
        #
        # df = pandas.DataFrame(simple_calls)
        # mind = (df == dataset.encoded_missing).sum(axis=1)
        # mind /= len(dataset.snps)  # Series
        #
        # for idx, min_data in enumerate(mind):
        #     if (1 - min_data) < threshold:
        #         silenced.silenced_sample(dataset.samples[idx].id)
        #
        #     if idx % 100 == 0:
        #         log.debug("Completed {} of {}".format(idx, len(mind)))

        return silenced


MinSampleDataFilter()
