# import json
import re

import logging

try:
    import simplejson as json
except ImportError:
    import json

from FilterResult import FilterResult

log = logging.getLogger(__file__)


class SNPDef:
    def __init__(self, clone_id: str, allele_id: str, sequence: str = None, sequence_ref: str = None,
                 rep_average: float = None, snp: str = None, all_headers: {} = None):

        self.clone_id = clone_id
        self.allele_id = allele_id
        self.sequence_ref = sequence_ref
        self.sequence = sequence
        self.snp = snp
        self.rep_average = rep_average
        self.all_headers = all_headers if all_headers is not None else {}

    def __repr__(self):
        return self.allele_id


class SampleDef:
    def __init__(self, id: str, population: str = "pop"):
        self.id = id
        self.population = population

    def __repr__(self):
        return self.id


class Dataset:
    homozygous_major = ("1", "0")
    homozygous_minor = ("0", "1")
    heterozygous = ("1", "1")
    missing = ("-", "-")

    encoded_heterozygous = "0"
    encoded_homozygous_minor = "1"
    encoded_homozygous_major = "2"
    encoded_missing = "-"

    def __init__(self, type: str, working_dir: str, batch_id: str, snps: [SNPDef], samples: [SampleDef],
                 calls: {} = None,
                 read_counts: {} = None):

        self.type = type

        # Validate that defs, calls and read counts all have matching SNP clone IDs & matching # samples
        val_errs = []
        if len(calls) != len(read_counts):
            val_errs.append("Calls and read counts have a different number of SNPs!  (input data error - continuable)")

        if len(snps) != len(read_counts):
            val_errs.append("SNP def mismatch - data has different number of SNPs than identified (input data error - continuable)")

        for snp in snps:
            if snp.allele_id not in calls:
                val_errs.append("Calls data missing SNP  (miss-named?): " + snp.allele_id)

            if snp.allele_id not in read_counts:
                val_errs.append("Read counts missing SNP (miss-named?): " + snp.allele_id)

            if len(calls[snp.allele_id]) != len(samples):
                val_errs.append("Call samples length incorrect for SNP: " + snp.allele_id)

            if len(read_counts[snp.allele_id]) != len(samples):
                val_errs.append("Read count samples length incorrect for SNP: " + snp.allele_id)

        if len(val_errs) > 0:
            log.error("Dataset read errors: \n{}\n\n".format("\n".join(val_errs)))

        self.batch_id = batch_id
        self.working_dir = working_dir

        self.snps = snps
        self.samples = samples

        self.calls = {}
        for allele_id, snp_calls in calls.items():
            self.calls[allele_id] = [
                tuple(call) if (call[0] == "1" or call[0] == "1") and (call[1] == "0" or call[1] == "1")
                else ("-", "-") for call in snp_calls]

        # Make sure counts are tuples (JSON dump/read converts to list)
        self.read_counts = read_counts
        for allele_id, snp_calls in self.read_counts.items():
            self.read_counts[allele_id] = [tuple(counts) for counts in snp_calls]

        # These are the results of the currently running filter - reset at the start of each filter run
        self.filtered = FilterResult()

    def clear_filters(self):
        self.filtered = FilterResult()

    def filter(self, filtered: FilterResult):
        self.filtered.samples.extend(filtered.samples)
        self.filtered.samples = list(set(self.filtered.samples))

        self.filtered.snps.extend(filtered.snps)
        self.filtered.snps = list(set(self.filtered.snps))

        for snp, calls in filtered.calls.items():
            if snp in self.filtered.calls:
                self.filtered.calls[snp].extend(calls)
                self.filtered.calls[snp] = list(set(self.filtered.calls[snp]))
            else:
                self.filtered.calls[snp] = calls

        for snp, changes in filtered.call_changes.items():
            if snp in self.filtered.call_changes:
                for sample_id, new_value in changes.items():
                    self.filtered.call_changes[snp][sample_id] = new_value
            else:
                self.filtered.call_changes[snp] = changes

    def get_filtered_calls(self) -> {str: []}:
        """
        Return a copy of self.calls which have all set calls silenced
        :param set_name: Name/idx of the parameter set to silence for
        :return: Calls dict {snp: [calls]}
        """
        filtered_calls = {}
        for snp, calls in self.calls.items():
            if snp not in self.filtered.snps:
                filtered_calls[snp] = [call if self.samples[idx].id not in self.filtered.samples and (
                    snp not in self.filtered.calls or self.samples[idx].id not in
                    self.filtered.calls[snp]) else ("-", "-") for idx, call in enumerate(calls)]

        return filtered_calls

    def is_snp_filtered(self, allele_id):
        return allele_id in self.filtered.snps

    def get_snp_def(self, allele_id):

        idx = 0
        snp_def = self.snps[0]
        while snp_def.allele_id != allele_id:
            idx += 1
            snp_def = self.snps[idx]

        if snp_def.allele_id != allele_id:
            return None

        return snp_def

    def get_sample_def(self, sample_id):
        idx = 0
        sample_def = self.samples[0]
        while sample_def.id != sample_id:
            idx += 1
            sample_def = self.samples[idx]

        if sample_def.id != sample_id:
            return None

        return sample_def

    def get_populations(self):
        pops = {}

        for idx, sample_def in enumerate(self.samples):
            if sample_def.population not in pops:
                pops[sample_def.population] = []

            pops[sample_def.population].append(idx)

        return pops

    @staticmethod
    def encode_single(calls: {str: []}):
        single_calls = {}

        for allele_id, tuple_calls in calls.items():
            single_calls[allele_id] = [Dataset.encoded_heterozygous if call == Dataset.heterozygous
                                       else Dataset.encoded_homozygous_minor if call == Dataset.homozygous_minor
            else Dataset.encoded_homozygous_major if call == Dataset.homozygous_major
            else Dataset.encoded_missing for call in tuple_calls]

        return single_calls

    def write_json(self, path):
        with open(path, "w") as outfile:
            jsonstr = json.dumps(self.__dict__, default=lambda o: o.__dict__)
            outfile.write(jsonstr)
            # json.dump(self.__dict__, default=lambda o: o.__dict__, fp=outfile)

    @staticmethod
    def read_json(path):
        with open(path) as dataset_file:
            dict_val = json.load(dataset_file)

            snps = []
            for snp_dict in dict_val["snps"]:
                snp = SNPDef(**snp_dict)
                snps.append(snp)

            samples = []
            for sample_dict in dict_val["samples"]:
                sample = SampleDef(*sample_dict.values())
                samples.append(sample)

            dataset = Dataset(dict_val["type"], dict_val["working_dir"], dict_val["batch_id"], snps, samples,
                              dict_val["calls"],
                              dict_val["read_counts"])

            del dict_val["snps"]
            del dict_val["samples"]
            del dict_val["calls"]
            del dict_val["read_counts"]

            for k, v in dict_val.items():
                setattr(dataset, k, v)

            return dataset

    def __repr__(self):
        return self.batch_id + " @ " + self.working_dir
