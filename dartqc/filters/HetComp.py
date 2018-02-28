import logging
import numpy
import re

import PipelineOptions
from Dataset import Dataset
from FilterResult import FilterResult

log = logging.getLogger(__file__)


class HetCompFilter(PipelineOptions.Filter):
    def get_order(self) -> int:
        return 5

    def get_cmd_type(self):
        return lambda s: [re.sub(r'[\(\)]', "", item).split(",")
                          for item in re.split(r"\),\(", re.sub(r'[\[\] ]', "", s))
                          if len(s.strip()) > 0]

    def get_name(self) -> str:
        return "count_comp"

    def get_cmd_help(self) -> str:
        return "Comparison of read counts between alleles.  Pattern: [(min_ratio, max_ratio), (eg. 0.05, 0.1),...]: < min ratio = homozygous, > max ratio = heterozygous, calls between are silenced"

    def filter(self, dataset: Dataset, threshold: (str, str), unknown_args: [], **kwargs) -> FilterResult:
        silenced = FilterResult()

        threshold[0] = float(threshold[0])
        threshold[1] = float(threshold[1])

        filtered_calls = dataset.get_filtered_calls()

        for allele_id, snp_counts in dataset.read_counts.items():
            # Identify which samples are silenced for this SNP
            for idx, sample_read_count in enumerate(snp_counts):
                if sample_read_count[0] == 0 or sample_read_count[1] == 0:
                    continue

                ratio = sample_read_count[0] / sample_read_count[1]
                if ratio > 1:
                    ratio = 1 / ratio

                if ratio < threshold[0] and filtered_calls[allele_id][idx] == Dataset.heterozygous:
                    # Homo
                    homo_call = Dataset.homozygous_minor
                    if sample_read_count[0] > sample_read_count[1]:
                        homo_call = Dataset.homozygous_major

                    silenced.add_call_change(allele_id, dataset.samples[idx].id, homo_call)
                elif ratio > threshold[1] and filtered_calls[allele_id][idx] != Dataset.heterozygous and filtered_calls[allele_id][idx] != Dataset.missing:
                    # Het
                    silenced.add_call_change(allele_id, dataset.samples[idx].id, Dataset.heterozygous)

                    log.warning("Homozygous call converted to heterozygous (1,1) based on read count ratio of {:0.3f}, counts: {:2d} & {:2d} original call: {}"
                                .format(ratio, sample_read_count[0], sample_read_count[1], dataset.calls[allele_id][idx]))
                else:
                    # Uncertain
                    silenced.silenced_call(allele_id, dataset.samples[idx].id)

        return silenced

HetCompFilter()
