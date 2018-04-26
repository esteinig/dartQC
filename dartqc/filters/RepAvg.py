import logging
import numpy
import re

import copy

from dartqc import PipelineOptions
from dartqc.Dataset import Dataset
from dartqc.FilterResult import FilterResult

log = logging.getLogger(__file__)


class RepAvgFilter(PipelineOptions.Filter):
    def get_name(self) -> str:
        return "Replication Average"

    def get_cmd_names(self):
        return ["--rep_average"]

    def get_cmd_help(self) -> str:
        return "Filter based on a average replication of duplicate samples.  " \
               "Pattern: [threshold, threshold, ...]. Where threshold can bo 0 to 1"

    def get_description(self) -> str:
        return "Filter based on a replication average of duplicate samples (calculated based on replicated read counts)"

    def get_cmd_type(self):
        return PipelineOptions.Filter.LIST_OF_FLOAT

    def get_order(self) -> int:
        return 6

    def filter(self, dataset: Dataset, threshold, unknown_args: [], **kwargs) -> FilterResult:
        silenced = FilterResult()

        rep_averages = RepAvgFilter.calc_rep_avg(dataset)

        for allele_id, rep_avg in rep_averages.items():
            if rep_avg < threshold:
                silenced.silenced_snp(allele_id)

        return silenced

    @staticmethod
    def calc_rep_avg(dataset: Dataset) -> {str: float}:
        # filtered_counts, filtered_snps, filtered_samples = dataset.get_filtered_counts()

        silenced_reps = [idx for idx, sample_id in enumerate(dataset.replicates) if sample_id in dataset.filtered.samples]

        rep_avgs = {}
        for allele_id, rep_counts in dataset.replicate_counts.items():
            if allele_id in dataset.filtered.snps:
                continue

            rep_avg = 0
            for idx, counts in enumerate(rep_counts):
                # Ignore replicates for samples that have been silenced!
                if idx in silenced_reps:
                    continue

                # Find the average value for each allele
                allele_1_avg = 0.0
                allele_2_avg = 0.0
                for rep_count in counts:
                    allele_1_avg += 0 if rep_count[0] == 0 else 1
                    allele_2_avg += 0 if rep_count[1] == 0 else 1

                allele_1_avg /= len(rep_count)
                allele_2_avg /= len(rep_count)

                # If there happened to be more than 2 replicates for a single sample, we want the match distance from the most popular call (0 or 1)
                # So the error is the closest distance from 0 or 1
                if allele_1_avg > 0.5:
                    allele_1_avg = 1 - allele_1_avg

                if allele_2_avg > 0.5:
                    allele_2_avg = 1 - allele_2_avg

                # Add the error for each allele
                rep_avg += allele_1_avg + allele_2_avg

            # Divide the tot error by # replicates to get % error
            rep_avg /= len(rep_counts)

            rep_avgs[allele_id] = 1 - rep_avg

        return rep_avgs

RepAvgFilter()
