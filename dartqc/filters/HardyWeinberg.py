from scipy.stats import stats

from Dataset import Dataset
from FilterResult import FilterResult
from PipelineOptions import Filter


class HWEFilter(Filter):
    def get_order(self) -> int:
        return 201

    def get_cmd_type(self):
        return lambda s: [float(item.strip()) if len(item.strip()) > 0 else None for item in s[1:-1].split(',')]

    def get_name(self) -> str:
        return "hwe"

    def get_cmd_help(self) -> str:
        return "filter snps <= p-value of hardy-weinberg test"

    def filter(self, dataset: Dataset, threshold: float, unknown_args: [], **kwargs) -> FilterResult:
        all_hwe_values = HWEFilter.calculate_hwe(dataset, True)

        # TODO:  Update to work with pops as per MAF

        silenced = FilterResult()
        filtered_calls = dataset.get_filtered_calls()

        pops = dataset.get_populations()

        for pop_name, pop_hwe in all_hwe_values.items():
            for allele_id in pop_hwe.keys():
                for idx, sample_def in enumerate(dataset.samples):
                    # Only silence calls that aren't already missing/silenced
                    # Only silence if the sample/call is from this population
                    if filtered_calls[allele_id][idx] != dataset.missing and idx in pops[pop_name] and pop_hwe[allele_id] < threshold:
                        silenced.silenced_call(allele_id, sample_def.id)

        return silenced

    @staticmethod
    def calculate_hwe(dataset: Dataset, use_pops: bool = True) -> {str: {str: float}}:
        """
        Calculates p-value for HWE using ChiSquare statistic: remove missing, get observed counts, get observed
        frequencies, get expected counts, calculate test values using (O-E)**2 / E and return ChiSquare probability
        with 1 degree of Freedom (bi-allelic SNP).

        """

        filtered_calls = dataset.get_filtered_calls()

        pops = list(range(len(dataset.samples)))
        if use_pops:
            pops = dataset.get_populations()

        hwe_values = {}
        for allele_id, snp_calls in filtered_calls.items():
            for pop_name, pop_sample_idxs in pops.items():
                pop_calls = [call for idx, call in enumerate(snp_calls) if idx in pop_sample_idxs]

                adjusted_samples = len(pop_calls) - pop_calls.count(dataset.missing)
                hetero_obs, major_obs, minor_obs = HWEFilter._get_observed(pop_calls)

            if adjusted_samples > 0:
                p = (major_obs + (hetero_obs / 2)) / adjusted_samples
                q = (minor_obs + (hetero_obs / 2)) / adjusted_samples

                if (p + q) != 1:
                    ValueError("Sum of observed allele frequencies (p + q) does not equal one.")

                hetero_exp, major_exp, minor_exp = HWEFilter._get_expected(p, q, adjusted_samples)

                if hetero_exp > 0:
                    hetero_test = ((hetero_obs - hetero_exp) ** 2) / hetero_exp
                    major_test = ((major_obs - major_exp) ** 2) / major_exp
                    minor_test = ((minor_obs - minor_exp) ** 2) / minor_exp

                    if pop_name not in hwe_values:
                        hwe_values[pop_name] = {}

                    hwe_values[pop_name][allele_id] = stats.chisqprob(sum([hetero_test, major_test, minor_test]), 1)

        return hwe_values

    @staticmethod
    def _get_observed(calls):
        """" Get observed counts in the genotype for each SNP """
        return calls.count(Dataset.heterozygous), calls.count(Dataset.homozygous_major), \
               calls.count(Dataset.homozygous_minor)

    @staticmethod
    def _get_expected(p, q, adjusted_samples):
        """ Get expected counts under HWE """
        return adjusted_samples * (2 * p * q), adjusted_samples * (p ** 2), adjusted_samples * (q ** 2)

HWEFilter()
