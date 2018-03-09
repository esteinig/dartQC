import json
import logging
import textwrap

# from dartqc.Dataset import Dataset

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

    def log(self, filter: str, threshold: float, dataset, fp, incl_summary: bool = False) -> None:
        tot_filtered_calls = sum([len(samples) for snp, samples in self.calls.items()])
        tot_calls = sum([len(samples) for snp, samples in dataset.calls.items()])

        tot_changed_calls = sum([len(samples) for snp, samples in self.call_changes.items()])

        filter_msg = textwrap.dedent("""                Filter: {} {}
                ----------------------------------------------
                Samples silenced:  {} of {} ({:.1f}%)
                SNPs silenced:     {} of {} ({:.1f}%)
                Calls silenced:    {} Calls of {} ({:.1f}%) across {} SNPs
                Calls changed:     {} Calls of {} ({:.1f}%) across {} SNPs
                """.format(filter, "at " + str(threshold) if threshold is not None else "",
                           len(self.samples), len(dataset.samples), len(self.samples) / len(dataset.samples),
                           len(self.snps), len(dataset.snps), len(self.snps) / len(dataset.snps),
                           tot_filtered_calls, tot_calls, tot_filtered_calls / tot_calls, len(self.calls),
                           tot_changed_calls, tot_calls,tot_changed_calls / tot_calls, len(self.call_changes)))

        if incl_summary:
            sum_all_silenced_calls = 0
            all_filtered_calls = {}
            for allele_id in self.snps:
                all_filtered_calls[allele_id] = [sample_def.id for sample_def in dataset.samples]
                sum_all_silenced_calls += len(dataset.samples) - dataset.calls[allele_id].tolist().count(
                    dataset.missing)

            for sample_id in self.samples:
                for snp_def in dataset.snps:
                    if snp_def.allele_id not in all_filtered_calls:
                        all_filtered_calls[snp_def.allele_id] = []

                    sample_idx = 0
                    while sample_id != dataset.samples[sample_idx].id:
                        sample_idx += 1

                    if sample_id not in all_filtered_calls[snp_def.allele_id] and dataset.calls[snp_def.allele_id][
                        sample_idx] != dataset.missing:
                        all_filtered_calls[snp_def.allele_id].append(sample_id)
                        sum_all_silenced_calls += 1

            for allele_id, sample_ids in self.calls.items():
                if allele_id not in all_filtered_calls:
                    all_filtered_calls[allele_id] = []

                for sample_id in sample_ids:
                    sample_idx = 0
                    while sample_id != dataset.samples[sample_idx].id:
                        sample_idx += 1

                    if sample_id not in all_filtered_calls[allele_id] \
                            and tuple(dataset.calls[allele_id][sample_idx]) != dataset.missing:
                        all_filtered_calls[allele_id].append(sample_id)
                        sum_all_silenced_calls += 1

            # if filter == "Final Results":
            #     test3 = 1

            tot_missing_calls = sum([samples.tolist().count(dataset.missing) for snp, samples in dataset.calls.items()])

            filter_msg += textwrap.dedent("""        
                    Total missing calls:  {} ({:.1f}%)
                    Total Calls silenced: {} of {} ({:.1f}%)
                    Remaining calls:      {} ({:.1f}%)
                    ----------------------------------------------
                    """.format(tot_missing_calls, tot_missing_calls / tot_calls * 100.0,
                               sum_all_silenced_calls, tot_calls, sum_all_silenced_calls / tot_calls * 100.0,
                               tot_calls - tot_missing_calls - sum_all_silenced_calls,
                               (tot_calls - tot_missing_calls - sum_all_silenced_calls) / tot_calls * 100.0))
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
