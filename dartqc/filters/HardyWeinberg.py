import re

import logging

import numpy
from scipy.stats import stats

from dartqc.Dataset import Dataset
from dartqc.FilterResult import FilterResult
from dartqc.PipelineOptions import Filter

log = logging.getLogger(__file__)


class HWEFilter(Filter):
    def get_order(self) -> int:
        return 201

    def get_cmd_type(self):
        return Filter.LIST_OF_LISTS

    def get_name(self) -> str:
        return "HWE"

    def get_cmd_names(self) -> [str]:
        return ["--hwe"]

    def get_cmd_help(self) -> str:
        return "filter snps <= p-value of hardy-weinberg test - Pattern: [(HWE thresh, # req. passing pops, ignored pop name(s), ...),...] if # req. passing pops is missing/None then pops are ignored."

    def get_description(self):
        return "filter snps <= p-value of hardy-weinberg test (0 to 1 decimal -> bigger means better HWE match)"

    def filter(self, dataset: Dataset, threshold: [str], unknown_args: [], **kwargs) -> FilterResult:

        # TODO:  Update to work with pops as per MAF
        hwe_thresh = float(threshold[0])

        no_pops = True
        req_success_pops = None
        ignored_pops = []

        if len(threshold) > 1 and threshold[1] is not None and threshold[1] != "None" \
                and threshold[1] != "" and threshold[1] != "null":
            req_success_pops = int(threshold[1])
            no_pops = False

            for pop in threshold[2:]:
                ignored_pops.append(pop)

        ignored_snps = numpy.asarray([True if snp_def.allele_id in dataset.filtered.snps else False
                                      for idx, snp_def in enumerate(dataset.snps)])

        all_hwe_values = HWEFilter.calculate_hwe(dataset, not no_pops)

        log.info("HWE values calculated - filtering")

        silenced = FilterResult()
        snp_pop_good_cnts = {snp.allele_id: 0 for snp in dataset.snps}

        # For all populations, how many exceed the required MAF threshold
        for pop_name, pop_hwe in all_hwe_values.items():
            if pop_name not in ignored_pops:
                for allele_id, hwe_val in pop_hwe.items():
                    snp_pop_good_cnts[allele_id] += 1 if hwe_val > hwe_thresh else 0

        # If less than required number of pops exceed the MAF threshold (threshold[0]) - silence the SNP
        for snp_idx, snp_def in enumerate(dataset.snps):
            # Ignore filtered SNPs
            if ignored_snps[snp_idx]:
                continue

            if req_success_pops is None or snp_pop_good_cnts[snp_def.allele_id] < req_success_pops:
                silenced.silenced_snp(snp_def.allele_id)

            if snp_idx % 5000 == 0:
                log.debug("Completed {} of {}".format(snp_idx, len(dataset.snps)))

        return silenced

    @staticmethod
    def calculate_hwe(dataset: Dataset, use_pops: bool = True, pops_blackilst: [str] = None) -> {str: {str: float}}:
        """
        Calculates p-value for HWE using ChiSquare statistic: remove missing, get observed counts, get observed
        frequencies, get expected counts, calculate test values using (O-E)**2 / E and return ChiSquare probability
        with 1 degree of Freedom (bi-allelic SNP).
        """

        log.info("Calculating HWE values")

        filtered_calls = dataset.get_filtered_calls()

        pops = {"pop": list(range(len(dataset.samples)))}
        if use_pops:
            pops = {k: v for k, v in dataset.get_populations().items()}

        # Remove blacklisted populations (ignore these samples)
        for pop_name in pops:
            if pops_blackilst is not None and pop_name in pops_blackilst:
                del pops[pop_name]

        ignored_snps = numpy.asarray([True if snp_def.allele_id in dataset.filtered.snps else False
                                      for idx, snp_def in enumerate(dataset.snps)])

        hwe_values = {pop: {} for pop in pops}
        for snp_idx, snp_def in enumerate(dataset.snps):
            # Ignore filtered SNPs
            if ignored_snps[snp_idx]:
                continue

            for pop_name, pop_sample_idxs in pops.items():
                pop_calls = filtered_calls[snp_def.allele_id][pop_sample_idxs]

                split_allele_calls = numpy.dstack(pop_calls)[0]


                # num_calls = len(pop_calls) - len(missing_idxs)
                allele_1_idxs = set(numpy.where(split_allele_calls[0] == "1")[0])
                allele_2_idxs = set(numpy.where(split_allele_calls[1] == "1")[0])

                if len(allele_1_idxs) + len(allele_2_idxs) > 0:
                    hetero_obs_idxs = allele_1_idxs & allele_2_idxs
                    num_hetero = len(hetero_obs_idxs)

                    # Missing = total - het - (allele 1 - het) - (allele 2 - het)
                    num_calls = len(allele_1_idxs) + len(allele_2_idxs) - num_hetero
                    num_missing = len(pop_sample_idxs) - num_calls

                    num_homo = len(pop_calls) - num_hetero - num_missing
                    num_major_obs = num_homo - len(allele_1_idxs)
                    num_minor_obs = num_homo - num_major_obs

                    p = (num_major_obs + (num_hetero / 2)) / num_calls
                    q = (num_minor_obs + (num_hetero / 2)) / num_calls

                    if (p + q) != 1:
                        ValueError("Sum of observed allele frequencies (p + q) does not equal one.")

                    # Get expected counts under HWE
                    hetero_exp = num_calls * (2 * p * q)
                    major_exp = num_calls * (p ** 2)
                    minor_exp = num_calls * (q ** 2)

                    if hetero_exp > 0:
                        hetero_test = ((num_hetero - hetero_exp) ** 2) / hetero_exp
                        major_test = ((num_major_obs - major_exp) ** 2) / major_exp
                        minor_test = ((num_minor_obs - minor_exp) ** 2) / minor_exp

                        hwe_values[pop_name][snp_def.allele_id] = stats.chisqprob(
                            sum([hetero_test, major_test, minor_test]), 1)

            if snp_idx % 5000 == 0:
                log.debug("Completed {} of {}".format(snp_idx, len(dataset.snps)))

        return hwe_values


HWEFilter()
