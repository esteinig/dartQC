import logging

import re

import numpy

from dartqc.Dataset import Dataset
from dartqc.FilterResult import FilterResult
from dartqc.PipelineOptions import Filter

log = logging.getLogger(__file__)


class MAFFilter(Filter):
    def get_order(self) -> int:
        return 200

    def get_cmd_type(self):
        return Filter.LIST_OF_LISTS

    def get_name(self) -> str:
        return "MAF"

    def get_cmd_names(self) -> [str]:
        return ["--maf"]

    def get_cmd_help(self) -> str:
        return "filter snps <= minor allele frequency - Pattern: [[MAF thresh, # req. passing pops, ignored pop name(s), ...],...] if # req. passing pops is missing/None then pops are ignored."

    def get_description(self):
        return "filter snps <= minor allele frequency (0 to 0.5 decimal where bigger means more heterozygous).  Make sure MAF value is larger than 1/# samples or it will do nothing"

    def filter(self, dataset: Dataset, threshold: [str], unknown_args: [], **kwargs) -> FilterResult:
        maf_thresh = float(threshold[0])
        no_pops = True
        req_success_pops = None
        ignored_pops = []

        # Convert the CMD line args into appropriate variables.
        if len(threshold) > 1 and threshold[1] is not None and threshold[1] != "None" \
                and threshold[1] != "" and threshold[1] != "null":
            req_success_pops = int(threshold[1])
            no_pops = False

            for pop in threshold[2:]:
                ignored_pops.append(pop)

        if maf_thresh < 1 / len(dataset.samples):
            log.warning("Filtering MAF at < (1 / # samples) is = filtering with MAF of 0 -> Does nothing...")

        # Calculate the MAF value for each SNP (using pops as per CMD args)
        all_maf = MAFFilter.calculate_maf(dataset, not no_pops, ignored_pops)

        filtered_snps = [snp_def for snp_def in dataset.snps if snp_def.allele_id not in dataset.filtered.snps]

        log.info("MAF values calculated - filtering")

        silenced = FilterResult()
        snp_pop_good_cnts = {snp.allele_id: 0 for snp in filtered_snps}

        # ignored_snps = numpy.asarray([True if dataset.snps[idx].allele_id in dataset.filtered.snps else False
        #                               for idx in range(len(dataset.snps))])

        # For all populations, how many exceed the required MAF threshold
        for pop_name, pop_maf in all_maf.items():
            if pop_name not in ignored_pops:
                for allele_id in pop_maf.keys():
                    snp_pop_good_cnts[allele_id] += 1 if pop_maf[allele_id] > maf_thresh else 0


        for snp_idx, snp_def in enumerate(filtered_snps):
            # Ignore filtered SNPs
            # if ignored_snps[snp_idx]:
            #     continue

            # If less than required number of pops exceed the MAF threshold (threshold[0]) - silence the SNP
            if (no_pops and snp_pop_good_cnts[snp_def.allele_id] < 1) or (not no_pops and snp_pop_good_cnts[snp_def.allele_id] < req_success_pops):
                silenced.silenced_snp(snp_def.allele_id)

            if snp_idx % 5000 == 0:
                log.debug("Completed {} of {}".format(snp_idx, len(filtered_snps)))

        return silenced

    @staticmethod
    def calculate_maf(dataset: Dataset, use_pops: bool = True, pops_blackilst: [str] = None, allele_list: [str] = None) -> {str: {str: float}}:
        """
        Calculates minor allele frequency for all SNPs in the dataset.
        Returns the minimum allele frequency for each SNP {clone_id: MAF}.
        """

        log.info("Calculating MAF values")

        filtered_calls, filtered_snps, filtered_samples = dataset.get_filtered_calls()

        # If not using populations, just add all samples as the default 'pop' population.
        pops = {"pop": list(range(len(filtered_samples)))}
        if use_pops:
            pops = {k: v for k, v in dataset.get_populations(filtered_samples).items()}

        # Remove blacklisted populations (ignore these samples)
        if pops_blackilst is not None:
            for pop_name in pops_blackilst:
                if pop_name in pops:
                    del pops[pop_name]

        ignored_snps = None
        if allele_list is not None:
            ignored_snps = numpy.asarray([True if filtered_snps[idx].allele_id not in allele_list else False
                                          for idx in range(len(filtered_snps))])

        maf_values = {pop: {} for pop in pops}
        for snp_idx, snp_def in enumerate(filtered_snps):
            # Ignore filtered SNPs
            if ignored_snps is not None and ignored_snps[snp_idx]:
                continue

            for pop_name, pop_sample_idxs in pops.items():
                pop_calls = filtered_calls[snp_def.allele_id][pop_sample_idxs]
                # pop_calls = numpy.asarray([filtered_calls[snp_def.allele_id][idx] for idx in pop_sample_idxs])
                split_allele_calls = numpy.dstack(pop_calls)[0]

                num_missing = len(numpy.where(split_allele_calls[0] == "-")[0])
                num_allele_1 = len(numpy.where(split_allele_calls[0] == "1")[0])
                num_allele_2 = len(numpy.where(split_allele_calls[1] == "1")[0])

                num_calls = (len(pop_calls) - num_missing) * 2

                if num_calls > 0:
                    freq_allele_one = num_allele_1 / num_calls
                    freq_allele_two = num_allele_2 / num_calls

                    maf_values[pop_name][snp_def.allele_id] = min(freq_allele_one, freq_allele_two)
                else:
                    maf_values[pop_name][snp_def.allele_id] = 0

            if snp_idx % 5000 == 0:
                log.debug("Completed {} of {}".format(snp_idx, len(filtered_snps)))

        return maf_values


MAFFilter()
