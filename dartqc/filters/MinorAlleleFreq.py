import logging

import re

from Dataset import Dataset
from FilterResult import FilterResult
from PipelineOptions import Filter

log = logging.getLogger(__file__)


class MAFFilter(Filter):
    def get_order(self) -> int:
        return 200

    def get_cmd_type(self):
        return lambda s: [re.sub(r'[\(\)]', "", item).split(",")
                          for item in re.split(r"\),\(", re.sub(r'[\[\] ]', "", s))
                          if len(s.strip()) > 0]

    def get_name(self) -> str:
        return "maf"

    def get_cmd_help(self) -> str:
        return "filter snps <= minor allele frequency"

    def filter(self, dataset: Dataset, threshold: (str, str), unknown_args: [], **kwargs) -> FilterResult:
        threshold[0] = float(threshold[0])
        threshold[1] = int(threshold[1])

        if threshold[0] < 1 / len(dataset.samples):
            log.warning("Filtering MAF at < (1 / # samples) is = filtering with MAF of 0 -> Does nothing...")

        all_maf = MAFFilter.calculate_maf(dataset, not kwargs["ignore_pops"], kwargs["pop_blacklist"])
        silenced = FilterResult()

        snp_pop_good_cnts = {snp.allele_id: 0 for snp in dataset.snps}

        # For all populations, how many exceed the required MAF threshold
        for pop_name, pop_maf in all_maf.items():
            for allele_id in pop_maf.keys():
                snp_pop_good_cnts[allele_id] += 1 if pop_maf[allele_id] > threshold[0] else 0

        # If less than required number of pops (threshold[1]) exceed the MAF threshold (threshold[0]) - silence the SNP
        for allele_id, count in snp_pop_good_cnts.items():
            if count < threshold[1]:
                silenced.silenced_snp(allele_id)

        return silenced

    @staticmethod
    def calculate_maf(dataset: Dataset, use_pops: bool = True, pops_blackilst: [str] = None) -> {str: {str: float}}:
        """
        Calculates minor allele frequency for all SNPs in the dataset.
        Returns the minimum allele frequency for each SNP {clone_id: MAF}.
        """

        filtered_calls = dataset.get_filtered_calls()

        # If not using populations, just add all samples as the default 'pop' population.
        pops = {"pop": list(range(len(dataset.samples)))}
        if use_pops:
            pops = dataset.get_populations()

        # Remove blacklisted populations (ignore these samples)
        for pop_name in pops:
            if pops_blackilst is not None and pop_name in pops_blackilst:
                del pops[pop_name]

        maf_values = {}
        for allele_id, snp_calls in filtered_calls.items():
            for pop_name, pop_sample_idxs in pops.items():
                if pop_name not in maf_values:
                    maf_values[pop_name] = {}

                pop_calls = [call for idx, call in enumerate(snp_calls) if idx in pop_sample_idxs]

                num_calls = len(pop_calls) - pop_calls.count(dataset.missing)
                het_count = pop_calls.count(dataset.heterozygous)

                if num_calls > 0:
                    freq_allele_one = (pop_calls.count(dataset.homozygous_major) + (het_count / 2)) / num_calls
                    freq_allele_two = (pop_calls.count(dataset.homozygous_minor) + (het_count / 2)) / num_calls

                    maf_values[pop_name][allele_id] = min(freq_allele_one, freq_allele_two)
                else:
                    maf_values[pop_name][allele_id] = 0

        return maf_values


MAFFilter()
