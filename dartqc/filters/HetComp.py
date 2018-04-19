import logging
import numpy
import re

from dartqc.PipelineOptions import Filter
from dartqc.Dataset import Dataset
from dartqc.FilterResult import FilterResult
import copy

log = logging.getLogger(__file__)


class HetCompFilter(Filter):
    def get_order(self) -> int:
        return 4

    def get_cmd_type(self):
        return Filter.LIST_OF_LISTS

    def get_name(self) -> str:
        return "Count Comparison"

    def get_cmd_names(self) -> [str]:
        return ["--count_comp"]

    def get_cmd_help(self) -> str:
        return "Comparison of read counts between alleles - Pattern: [[min_ratio, max_ratio],...] where < min ratio = homozygous, > max ratio = heterozygous, calls between are silenced"

    def get_description(self):
        return "Comparison of read counts between alleles (0 to 1 decimal values where < min ratio = homozygous, > max ratio = heterozygous, calls between are silenced"

    def filter(self, dataset: Dataset, threshold: [str], unknown_args: [], **kwargs) -> FilterResult:
        silenced = FilterResult()

        min_ratio = float(threshold[0])
        max_ratio = float(threshold[1])

        filtered_calls, filtered_snps, filtered_samples = dataset.get_filtered_calls()

        # ignored_snps = numpy.asarray([True if dataset.snps[idx].allele_id in dataset.filtered.snps else False
        #                               for idx in range(len(filtered_snps))])

        # TODO: Convert to do maths across whole matrix

        for snp_idx, snp_def in enumerate(filtered_snps):
            # if ignored_snps[snp_idx]:
            #     continue

            allele_id = snp_def.allele_id
            snp_counts = dataset.read_counts[snp_def.allele_id]

            sample_idxs = {sample_def.id: idx for idx, sample_def in enumerate(filtered_samples)}
            del_idxs = []
            for idx, sample_def in enumerate(dataset.samples):
                if sample_def.id not in sample_idxs:
                    del_idxs.append(idx)

            snp_counts = numpy.delete(snp_counts, del_idxs, 0)


            # first_allele_cnts = snp_counts[:, 0]
            # second_allele_cnts = snp_counts[:, 0]
            #
            # ratios = numpy.divide(first_allele_cnts, second_allele_cnts)
            # ratios = numpy.asarray([0 if numpy.isnan(ratio) or numpy.isinf(ratio)
            #                         else ratio if ratio < 1 else 1 / ratio for ratio in ratios])
            #
            # silence_idxs = numpy.where((ratios > min_ratio) & (ratios < max_ratio))[0].tolist()
            # silenced.calls[allele_id] = [filtered_samples[idx].id for idx in silence_idxs]
            #
            # homo_idxs = numpy.where(ratios < min_ratio)[0].tolist()
            # homo_idxs = [idx for idx in homo_idxs if tuple(dataset.calls[allele_id][idx]) == Dataset.heterozygous]
            # silenced.call_changes[allele_id] = {filtered_samples[idx].id: Dataset.homozygous_major
            #                                     if first_allele_cnts[idx] > second_allele_cnts[idx]
            #                                     else Dataset.homozygous_minor for idx in silence_idxs}

            # Identify which samples are silenced for this SNP
            for idx, sample_read_count in enumerate(snp_counts):
                # Ignore anything that is already silenced/missing
                if filtered_calls[allele_id][idx][0] == "-":  # Assume that silenced calls have both set as -
                    continue

                # Shortcut for any with a 0 read count (must be homo) - this is common (+NaN) so prioritise speed.
                if sample_read_count[0] == 0 or sample_read_count[1] == 0:
                    if filtered_calls[allele_id][idx][0] == "1" and filtered_calls[allele_id][idx][1] == "1":
                        if sample_read_count[0] == 0 and sample_read_count[1] == 0:
                            if allele_id not in silenced.calls:
                                silenced.calls[allele_id] = []

                            silenced.calls[allele_id].append(filtered_samples[idx].id)
                        elif sample_read_count[0] > sample_read_count[1]:
                            if allele_id not in silenced.call_changes:
                                silenced.call_changes[allele_id] = {}

                            silenced.call_changes[allele_id][filtered_samples[idx].id] = Dataset.homozygous_minor
                        else:
                            if allele_id not in silenced.call_changes:
                                silenced.call_changes[allele_id] = {}

                            silenced.call_changes[allele_id][filtered_samples[idx].id] = Dataset.homozygous_minor

                    continue

                # Calculate the read count ratios (0's are removed so no NaN probs)
                ratio = sample_read_count[0] / sample_read_count[1]
                if ratio > 1:
                    ratio = 1 / ratio

                if ratio < min_ratio and filtered_calls[allele_id][idx][0] == "1" \
                        and filtered_calls[allele_id][idx][1] == "1":
                    # Homo
                    homo_call = Dataset.homozygous_minor
                    if sample_read_count[0] > sample_read_count[1]:
                        homo_call = Dataset.homozygous_major

                    silenced.add_call_change(allele_id, filtered_samples[idx].id, homo_call)
                elif ratio > max_ratio \
                        and (filtered_calls[allele_id][idx][0] != "1" or filtered_calls[allele_id][idx][1] != "1"):

                    homo_replicate = False
                    rep_counts = sample_read_count
                    if filtered_samples[idx].id in dataset.replicates:
                        rep_counts = dataset.replicate_counts[allele_id][dataset.replicates.index(filtered_samples[idx].id)]

                        for counts in rep_counts:
                            if counts[0] == 0 or counts[1] == 0:
                                homo_replicate = True

                    if not homo_replicate:
                        # Het
                        silenced.add_call_change(allele_id, filtered_samples[idx].id, Dataset.heterozygous)

                        log.warning(
                            "Homozygous call converted to heterozygous (1,1) based on read count ratio of {:0.3f},"
                            " total counts: {:2d} & {:2d} replicate counts: {}  original call: {} - {}:{}"
                            .format(ratio, sample_read_count[0], sample_read_count[1],
                                    rep_counts, dataset.calls[allele_id][idx], allele_id,
                                    filtered_samples[idx].id))
                else:
                    # Uncertain
                    silenced.silenced_call(allele_id, filtered_samples[idx].id)

            if snp_idx % 5000 == 0:
                log.debug("Completed {} of {}".format(snp_idx, len(filtered_snps)))

        return silenced


HetCompFilter()
