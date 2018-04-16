import json
import logging
import textwrap

# from dartqc.Dataset import Dataset
import numpy

log = logging.getLogger(__file__)


# TODO: Update read/write to CSV - better to read + should be smaller

class FilterResult:
    def __init__(self):
        self.samples = []
        self.snps = []
        self.calls = {}

        # Key should be a callable function
        # Valye should be function parameters
        self.call_changes = {}

    def silenced_sample(self, sample_id: str):
        if type(sample_id) != str:
            log.error("Silencing SampleDef instead of sample ID (str)!  This is a code error that needs to be fixed.")

        if sample_id not in self.samples:
            self.samples.append(sample_id)

    def silenced_snp(self, allele_id: str):
        if type(allele_id) != str:
            log.error("Silencing SNPDef instead of SNP allele_id (str)!  This is a code error that needs to be fixed.")

        if allele_id not in self.snps:
            self.snps.append(allele_id)

    def silenced_call(self, allele_id: str, sample_id: str):
        # if type(sample_id) != str:
        #     log.error(
        #         "Silencing call using SampleDef instead of sample ID (str)!  This is a code error that needs to be fixed.")
        #
        # if type(allele_id) != str:
        #     log.error(
        #         "Silencing call using SNPDef instead of SNP allele_id (str)!  This is a code error that needs to be fixed.")

        if allele_id not in self.calls:
            self.calls[allele_id] = []

        if sample_id not in self.calls[allele_id]:
            self.calls[allele_id].append(sample_id)

    def add_call_change(self, allele_id, sample_id, new_value: []):
        if type(sample_id) != str:
            log.error(
                "Changing call using SampleDef instead of sample ID (str)!  This is a code error that needs to be fixed.")

        if type(allele_id) != str:
            log.error(
                "Changing call using SNPDef instead of SNP allele_id (str)!  This is a code error that needs to be fixed.")

        if allele_id not in self.call_changes:
            self.call_changes[allele_id] = {}

        self.call_changes[allele_id][sample_id] = new_value

    def is_call_filtered(self, allele_id, sample_id):
        return allele_id in self.snps or sample_id in self.samples or (allele_id in self.calls and sample_id in self.calls[allele_id])

    def get_tot_silenced_calls(self):
        return sum([len(silenced_calls) for snp, silenced_calls in self.calls.items()])

    def get_tot_changed_calls(self):
        return sum([len(samples) for snp, samples in self.call_changes.items()])

    def log(self, filter: str, threshold: float, dataset, fp, incl_summary: bool = False) -> None:
        if not incl_summary:
            filtered_snps = [snp_def for snp_def in dataset.snps if snp_def.allele_id not in dataset.filtered.snps]
            filtered_samples = [sample_def for sample_def in dataset.samples if sample_def.id not in dataset.filtered.samples]
        else:
            filtered_snps = dataset.snps
            filtered_samples = dataset.samples

        tot_filtered_calls = sum([len(silenced_calls) for snp, silenced_calls in self.calls.items()])
        tot_calls = len(filtered_samples) * len(filtered_snps)

        tot_changed_calls = sum([len(samples) for snp, samples in self.call_changes.items()])

        filter_msg = textwrap.dedent("""                Filter: {} {}
                ----------------------------------------------
                Samples silenced:  {} of {} ({:.1f}%)
                SNPs silenced:     {} of {} ({:.1f}%)
                Calls silenced:    {} Calls of {} ({:.1f}%) across {} SNPs
                Calls changed:     {} Calls of {} ({:.1f}%) across {} SNPs
                """.format(filter, "at " + str(threshold) if threshold is not None else "",
                           len(self.samples), len(filtered_samples), 0 if len(filtered_samples) == 0 else len(self.samples) / len(filtered_samples) * 100.0,
                           len(self.snps), len(filtered_snps), 0 if len(filtered_snps) == 0 else len(self.snps) / len(filtered_snps) * 100.0,
                           tot_filtered_calls, tot_calls, 0 if tot_calls == 0 else tot_filtered_calls / tot_calls * 100.0, len(self.calls),
                           tot_changed_calls, tot_calls,0 if tot_calls == 0 else tot_changed_calls / tot_calls * 100.0, len(self.call_changes)))

        # log.debug("Test 1")
        if incl_summary:
            tot_calls = len(dataset.snps) * len(dataset.samples)

            all_missing_idxs = numpy.asarray([dataset.calls[snp.allele_id] for snp in dataset.snps])
            all_missing_idxs = numpy.dstack(all_missing_idxs)  # [SNPs][samples][calls] -> [samples][calls][SNPs]
            all_missing_idxs = numpy.dstack(all_missing_idxs)  # [SNPs][calls][samples] -> [calls][SNPs][samples]
            all_missing_idxs = all_missing_idxs[0]  # Only get first allele calls (if "-" -> missing)
            all_missing_idxs = [numpy.where(all_missing_idxs[snp_idx] == "-")[0] for snp_idx, snp_def in enumerate(dataset.snps)]

            sample_idxs = {sample_def.id: idx for idx, sample_def in enumerate(dataset.samples)}
            snp_idx_map = {snp_def.allele_id: idx for idx, snp_def in enumerate(dataset.snps)}
            num_samples = len(dataset.samples)

            # log.debug("Test 1.5")

            # Find all missing calls
            tot_missing_calls = 0
            for snp_idx, snp_def in enumerate(dataset.snps):
                tot_missing_calls += len(all_missing_idxs[snp_idx])

            # log.debug("Test 4")

            # Add count of silenced SNPs
            sum_all_silenced_calls = 0
            # all_filtered_calls = {}
            for allele_id in self.snps:
                # all_filtered_calls[allele_id] = [sample_def.id for sample_def in filtered_samples]
                sum_all_silenced_calls += num_samples - len(all_missing_idxs[snp_idx_map[allele_id]])

            # log.debug("Test 2")

            #  Add count of silenced samples
            silenced_sample_idxs = [sample_idxs[sample_id] for sample_id in self.samples]
            if len(silenced_sample_idxs) > 0:
                for snp_def in dataset.snps:
                    if snp_def.allele_id not in self.snps:
                        sum_all_silenced_calls += len(silenced_sample_idxs) - len(numpy.intersect1d(silenced_sample_idxs, all_missing_idxs[snp_idx_map[snp_def.allele_id]]))

                        # silenced_slice = numpy_matrix[snp_def.allele_id][silenced_sample_idxs]
                        # sum_all_silenced_calls += len(silenced_slice) - len(numpy.where(silenced_slice == "-")[0])

            # log.debug("Test 3")

            # Add count of silenced calls
            for allele_id, sample_ids in self.calls.items():
                if allele_id not in self.snps:
                    call_idxs = [sample_idxs[sample_id] for sample_id in sample_ids if sample_id not in self.samples]

                    sum_all_silenced_calls += len(call_idxs) - len(numpy.intersect1d(call_idxs, all_missing_idxs[snp_idx_map[allele_id]]))

                    # silenced_slice = numpy_matrix[allele_id][call_idxs]
                    # sum_all_silenced_calls += len(numpy.where(silenced_slice != "-")[0])


            # log.debug("Test 5")
            filter_msg += textwrap.dedent("""        
                    Total missing calls:  {} ({:.1f}%)
                    Total Calls silenced: {} of {} ({:.1f}%)
                    Remaining calls:      {} ({:.1f}%)
                    ----------------------------------------------
                    """.format(tot_missing_calls, tot_missing_calls / tot_calls * 100.0,
                               sum_all_silenced_calls, tot_calls, sum_all_silenced_calls / tot_calls * 100.0,
                               tot_calls - tot_missing_calls - sum_all_silenced_calls,
                               (tot_calls - tot_missing_calls - sum_all_silenced_calls) / tot_calls * 100.0))
            # log.debug("Test 6")
        else:
            filter_msg += "----------------------------------------------\n"

        log.info(filter_msg)

        if fp is not None:
            fp.write(filter_msg)

    def write_json(self, path):
        with open(path, "w") as outfile:
            jsonstr = json.dumps(self.__dict__, default=lambda o: o.__dict__, indent=4)
            outfile.write(jsonstr)

    @staticmethod
    def read_json(path):
        with open(path) as results_file:
            dict_val = json.load(results_file)

            results = FilterResult()
            # results.snps = dict_val["snps"]
            # results.samples = dict_val["samples"]
            # results.calls = dict_val["calls"]
            # results.call_changes = dict_val["call_changes"]

            for k, v in dict_val.items():
                setattr(results, k, v)

            return results
