import logging
import numpy

from dartqc.Dataset import Dataset
from dartqc.PipelineOptions import Filter
from dartqc.FilterResult import FilterResult

log = logging.getLogger(__file__)


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

        # numpy_matrix = numpy.asarray([dataset.read_counts[snp.allele_id] for snp in dataset.snps])
        # numpy_matrix = numpy.sum(numpy_matrix, axis=2)
        # numpy_matrix = [numpy.where(row <= threshold)[0] for row in numpy_matrix]

        ignored_snps = numpy.asarray([True if dataset.snps[idx].allele_id in dataset.filtered.snps else False
                                      for idx in range(len(dataset.snps))])

        for snp_idx, snp_def in enumerate(dataset.snps):
            # Ignore filtered SNPs
            if ignored_snps[snp_idx]:
                continue

            # Note: Using the commented numpy_matrix way may be slightly faster at expense of some memory
            fail_idxs = numpy.sum(dataset.read_counts[snp_def.allele_id], axis=1)
            fail_idxs = numpy.where(fail_idxs[snp_idx] <= threshold)[0].tolist()

            # fail_idxs = numpy_matrix[snp_idx]

            # Silence the calls
            for idx in fail_idxs:
                sample_id = dataset.samples[idx].id

                # if sample_id not in dataset.filtered.samples and snp_def.allele_id not in dataset.filtered.snps \
                #     and (snp_def.allele_id not in dataset.filtered.calls or sample_id not in dataset.filtered.calls[snp_def.allele_id]) \
                call = dataset.calls[snp_def.allele_id][idx]
                if call[0] != "-" and call[1] != "-":
                    # silenced.silenced_call(snp_def.allele_id, dataset.samples[idx].id)
                    if snp_def.allele_id not in silenced.calls:
                        silenced.calls[snp_def.allele_id] = []

                    silenced.calls[snp_def.allele_id].append(sample_id)

            if snp_idx % 5000 == 0:
                log.debug("Completed {} of {}".format(snp_idx, len(dataset.snps)))

        # Remove any that were already filtered previously
        for allele_id in dataset.filtered.snps:
            if allele_id in silenced.calls:
                del silenced.calls[allele_id]

        for sample_id in dataset.filtered.samples:
            for allele_id, samples in silenced.calls.items():
                samples.remove(sample_id)

        for allele_id, samples in dataset.filtered.calls:
            if allele_id in silenced.calls:
                silenced.calls[allele_id] = [sample_id for sample_id in silenced.calls[allele_id] if sample_id not in samples]

        return silenced


ReadCountsFilter()
