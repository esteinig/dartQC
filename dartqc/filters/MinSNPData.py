import logging
import numpy
import re

from dartqc.FilterResult import FilterResult
from dartqc.Dataset import Dataset
from dartqc.PipelineOptions import Filter

log = logging.getLogger(__file__)


class MinSNPDataFilter(Filter):
    def get_order(self) -> int:
        return 11

    def get_cmd_type(self):
        return Filter.LIST_OF_FLOAT

    def get_name(self) -> str:
        return "Min SNP Data"

    def get_cmd_names(self) -> [str]:
        return ["--min_snp","--call_rate"]

    def get_cmd_help(self) -> str:
        return "filter snps <= call rate of snp"

    def get_description(self) -> str:
        return "filter snps <= call rate of snp (0 to 1 decimal -> bigger requires more data)"

    def filter(self, dataset: Dataset, threshold: float, unknown_args: [], **kwargs) -> FilterResult:
        silenced = FilterResult()

        filtered_calls = dataset.get_filtered_calls()

        # Get calls as matrix & trim back to only first allele's call (much quicker processing)
        numpy_matrix = numpy.asarray([filtered_calls[snp.allele_id] for snp in dataset.snps])
        numpy_matrix = numpy.dstack(numpy_matrix)  # [SNPs][samples][calls] -> [samples][calls][SNPs]
        numpy_matrix = numpy.dstack(numpy_matrix)  # [SNPs][calls][samples] -> [calls][SNPs][samples]
        numpy_matrix = numpy_matrix[0]  # Only get first allele calls (if "-" -> missing)

        ignored_snps = numpy.asarray([True if dataset.snps[idx].allele_id in dataset.filtered.snps else False
                                      for idx in range(len(dataset.snps))])

        for snp_idx, snp_def in enumerate(dataset.snps):
            # Ignore filtered SNPs
            if ignored_snps[snp_idx]:
                continue

            # first_allele_calls = numpy.dstack(filtered_calls[snp_def.allele_id])[0][0]
            missing_idxs = numpy.where(numpy_matrix[snp_idx] == "-")[0]
            perc_data = 1 - len(missing_idxs) / len(dataset.samples)
            if perc_data < threshold:
                silenced.silenced_snp(snp_def.allele_id)

            if snp_idx % 5000 == 0:
                log.debug("Completed {} of {}".format(snp_idx, len(dataset.snps)))

        return silenced

MinSNPDataFilter()
