import logging

import numpy
import pandas

from dartqc.Dataset import Dataset
from dartqc.FilterResult import FilterResult
from dartqc.PipelineOptions import Filter

log = logging.getLogger(__file__)


class MinSampleDataFilter(Filter):
    def get_name(self):
        return "min_sample"

    def get_alt_cmd_names(self) -> [str]:
        return ["--mind"]

    def get_order(self) -> int:
        return 10

    def get_cmd_type(self):
        return lambda s: [float(item.strip()) if len(item.strip()) > 0 else None for item in s[1:-1].split(',')]

    def get_cmd_help(self) -> str:
        return "filter samples > missingness per sample"

    def filter(self, dataset: Dataset, threshold: float, unknown_args: [], **kwargs) -> FilterResult:
        silenced = FilterResult()

        filtered_calls = dataset.get_filtered_calls()

        # Get calls as matrix - swap SNP/sample axes & trim back to only first allele's call (much quicker processing)
        numpy_matrix = numpy.asarray([filtered_calls[snp.allele_id] for snp in dataset.snps])
        numpy_matrix = numpy.stack(numpy_matrix, axis=2)  # [SNPs][samples][calls] -> [SNPs][calls][samples]
        numpy_matrix = numpy.stack(numpy_matrix, axis=1)  # [SNPs][calls][samples] -> [calls][samples][SNPs]
        numpy_matrix = numpy_matrix[0]  # Only get first allele calls (if "-" -> missing)

        ignored_samples = numpy.asarray([True if dataset.samples[idx] in dataset.filtered.samples else False
                                         for idx in range(len(dataset.samples))])

        for idx, sample_calls in enumerate(numpy_matrix):
            if ignored_samples[idx]:
                continue

            # first_allele_calls = numpy.dstack(sample_calls)[0][0]  # Get just the first allele (if "-" => missing)
            missing_idxs = numpy.where(sample_calls == "-")[0]
            missing = len(missing_idxs) / len(dataset.snps)

            if (1 - missing) < threshold:
                silenced.silenced_sample(dataset.samples[idx].id)

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
