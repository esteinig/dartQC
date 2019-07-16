# import json
import re

import logging

import sys

import numpy
import time

try:
    import simplejson as json
except ImportError:
    import json

from dartqc.FilterResult import FilterResult

import copy

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
    def __init__(self, id: str, population: str = "pop", all_headers: {} = None):
        self.id = id
        self.population = population
        self.all_headers = all_headers if all_headers is not None else {}

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

    data_cols = ["calls", "read_counts", "replicate_counts"]

    def __init__(self, type: str, working_dir: str, batch_id: str, snps: [SNPDef], samples: [SampleDef],
                 calls: {}, read_counts: {}, replicates: [str], replicate_counts: {str: [[]]}):

        self.type = type

        self.batch_id = batch_id
        self.working_dir = working_dir

        self.snps = snps
        self.samples = samples

        self.calls = calls
        self.read_counts = read_counts

        self.replicates = replicates
        self.replicate_counts = replicate_counts

        # These are the results of the currently running filter - reset at the start of each filter run
        self.filtered = FilterResult()

        # Validate this dataset if it is newly created (don't validate yet if loading from file)
        if calls is not None or read_counts is not None or replicate_counts is not None:
            self._validate_input_data()

    def _validate_input_data(self):
        # Validate that defs, calls and read counts all have matching SNP clone IDs & matching # samples
        val_errs = []
        if len(self.calls) != len(self.read_counts):
            val_errs.append("Calls and read counts have a different number of SNPs!  (input data error - continuable)")

        if len(self.snps) != len(self.read_counts):
            val_errs.append(
                "SNP def mismatch - data has different number of SNPs than identified (input data error - continuable)")

        for snp in self.snps:
            if snp.allele_id not in self.calls:
                val_errs.append("Calls data missing SNP  (miss-named?): " + snp.allele_id)

            if snp.allele_id not in self.read_counts:
                val_errs.append("Read counts missing SNP (miss-named?): " + snp.allele_id)

            if len(self.calls[snp.allele_id]) != len(self.samples):
                val_errs.append("Call samples length incorrect for SNP: " + snp.allele_id + " - " + str(len(self.calls[snp.allele_id])) + " vs " + str(len(self.samples)))

            if len(self.read_counts[snp.allele_id]) != len(self.samples):
                val_errs.append("Read count samples length incorrect for SNP: " + snp.allele_id + " - " + str(len(self.read_counts[snp.allele_id])) + " vs " + str(len(self.samples)))

        if len(val_errs) > 0:
            log.error("Dataset read errors: \n{}\n\n".format("\n".join(val_errs)))

        if len(self.replicates) != len(set(self.replicates)):
            log.error("Duplicate replicate IDs - dataset reading error (each replicate name should only be given once)")

    def clear_filters(self):
        self.filtered = FilterResult()

    def filter(self, filtered: FilterResult):
        start = time.time()

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

                # log.info("Time to add filters to dataset: {:.1f}".format(time.time() - start))

    def get_filtered_calls(self) -> {str: []}:
        """
        Return a copy of self.calls which have all set calls silenced
        :param set_name: Name/idx of the parameter set to silence for
        :return: Calls dict {snp: [calls]}
        """

        start = time.time()

        # Get a copy of orig calls.
        # filtered_calls = numpy.copy(self.calls).tolist()
        filtered_calls = copy.deepcopy(self.calls)
        filtered_snps = copy.copy(self.snps)
        filtered_samples = copy.copy(self.samples)

        # Silence whole SNPs
        for allele_id in self.filtered.snps:
            snp_idx = None

            for idx, snp_def in enumerate(filtered_snps):
                if snp_def.allele_id == allele_id:
                    snp_idx = idx
                    break

            if snp_idx is None:
                log.error("Silenced SNP not found!")
                continue

            del filtered_calls[allele_id]
            del filtered_snps[snp_idx]

            # Need to actually remove from dataset - otherwise it creates invalid missingness!
            # filtered_calls[allele_id] = numpy.asarray([("-", "-") for idx in range(len(self.samples))])

        # Silence whole samples
        sample_idxs = {sample.id: idx for idx, sample in enumerate(self.samples)}
        del_sample_idxs = []
        for sample_id in self.filtered.samples:
            sample_idx = sample_idxs[sample_id]

            if sample_idx is None:
                log.error("Silenced Sample not found!")
                continue

            del_sample_idxs.append(sample_idx)

        for idx in sorted(del_sample_idxs, reverse=True):
            del filtered_samples[idx]

        # Delete samples from numpy (do at end with array for performance)
        for allele_id, calls in filtered_calls.items():
            filtered_calls[allele_id] = numpy.delete(calls, del_sample_idxs, 0)

            # Delete the sample column entirely - otherwise it creates invalid missingness
            # for allele_id in filtered_calls:
            #     filtered_calls[allele_id][sample_idxs[sample_id]] = Dataset.missing

        # Silence calls
        sample_idxs = {sample.id: idx for idx, sample in enumerate(filtered_samples)}
        snp_idxs = {snp_def.allele_id: idx for idx, snp_def in enumerate(filtered_snps)}
        for allele_id, samples in self.filtered.calls.items():
            # Check that this SNP hasn't already been silenced
            if allele_id in snp_idxs:
                for sample_id in samples:
                    # Check that this sample hasn't been silenced
                    if sample_id in sample_idxs:
                        filtered_calls[allele_id][sample_idxs[sample_id]] = Dataset.missing

        # Change calls
        for allele_id, changes in self.filtered.call_changes.items():
            # Check that this SNP hasn't already been silenced
            if allele_id in snp_idxs:
                for sample_id, new_val in changes.items():
                    # Check that this sample hasn't been silenced
                    if sample_id in sample_idxs:
                        filtered_calls[allele_id][sample_idxs[sample_id]] = new_val

        # filtered_calls = {}
        # for snp, calls in self.calls.items():
        #     if snp not in self.filtered.snps:
        #         filtered_calls[snp] = numpy.asarray([call if self.samples[idx].id not in self.filtered.samples and (
        #             snp not in self.filtered.calls or self.samples[idx].id not in
        #             self.filtered.calls[snp]) else ("-", "-") for idx, call in enumerate(calls)])
        #     else:
        #         filtered_calls[snp] = numpy.asarray([["-", "-"] for sample in self.samples])

        # log.info("Time to get filtered dataset: {:.1f}".format(time.time() - start))

        return filtered_calls, filtered_snps, filtered_samples

    def get_filtered_counts(self) -> ():
        filtered_read_counts = copy.deepcopy(self.read_counts)
        filtered_snps = copy.copy(self.snps)
        filtered_samples = copy.copy(self.samples)

        for allele_id in self.filtered.snps:
            snp_idx = None

            for idx, snp_def in enumerate(filtered_snps):
                if snp_def.allele_id == allele_id:
                    snp_idx = idx
                    break

            if snp_idx is None:
                log.error("Silenced SNP not found!")
                continue

            del filtered_read_counts[allele_id]
            del filtered_snps[snp_idx]

        # Silenece whole samples
        sample_idxs = {sample.id: idx for idx, sample in enumerate(self.samples)}
        del_sample_idxs = []
        for sample_id in self.filtered.samples:
            sample_idx = sample_idxs[sample_id]

            if sample_idx is None:
                log.error("Silenced Sample not found!")
                continue

            del_sample_idxs.append(sample_idx)

        for idx in sorted(del_sample_idxs, reverse=True):
            del filtered_samples[idx]

        # Delete samples from numpy (do at end with array for performance)
        for allele_id, counts in filtered_read_counts.items():
            filtered_read_counts[allele_id] = numpy.delete(counts, del_sample_idxs, 0)

        sample_idxs = {sample.id: idx for idx, sample in enumerate(filtered_samples)}
        for allele_id, samples in self.filtered.calls.items():
            if allele_id in filtered_read_counts:
                for sample_id in samples:
                    if sample_id in sample_idxs:
                        filtered_read_counts[allele_id][sample_idxs[sample_id]] = (0, 0)

        return filtered_read_counts, filtered_snps, filtered_samples

    # def get_filtered_calls(self) -> {str: []}:
    #     """
    #     Return a copy of self.calls which have all set calls silenced
    #     :param set_name: Name/idx of the parameter set to silence for
    #     :return: Calls dict {snp: [calls]}
    #     """
    #
    #     start = time.time()
    #
    #     # Get a copy of orig calls.
    #     filtered_calls = numpy.copy(self.calls).tolist()
    #
    #     # Silence whole SNPs
    #     for allele_id in self.filtered.snps:
    #         filtered_calls[allele_id] = numpy.asarray([("-", "-") for idx in range(len(self.samples))])
    #
    #     # Silenece whole samples
    #     sample_idxs = {sample.id: idx for idx, sample in enumerate(self.samples)}
    #     for sample_id in self.filtered.samples:
    #         # Delete the sample column entirely - otherwise it creates invalid missingness
    #         for allele_id in filtered_calls:
    #             filtered_calls[allele_id][sample_idxs[sample_id]] = Dataset.missing
    #
    #     # Silence calls
    #     for allele_id, samples in self.filtered.calls.items():
    #         for sample_id in samples:
    #             filtered_calls[allele_id][sample_idxs[sample_id]] = Dataset.missing
    #
    #     # Change calls
    #     for allele_id, changes in self.filtered.call_changes.items():
    #         for sample_id, new_val in changes.items():
    #             filtered_calls[allele_id][sample_idxs[sample_id]] = new_val
    #
    #     log.info("Time to get filtered dataset: {:.1f}".format(time.time() - start))
    #
    #     return filtered_calls #, filtered_snps, filtered_samples

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

    def get_populations(self, samples=None):
        pops = {}

        if samples is None:
            samples = self.samples

        for idx, sample_def in enumerate(samples):
            if sample_def.population not in pops:
                pops[sample_def.population] = []

            pops[sample_def.population].append(idx)

        return pops

    @staticmethod
    def encode_single(calls: {str: []}):
        single_calls = {}

        for allele_id, tuple_calls in calls.items():
            single_calls[allele_id] = [Dataset.encoded_heterozygous if tuple(call) == Dataset.heterozygous
                                       else Dataset.encoded_homozygous_minor if tuple(call) == Dataset.homozygous_minor
            else Dataset.encoded_homozygous_major if tuple(call) == Dataset.homozygous_major
            else Dataset.encoded_missing for call in tuple_calls]

        return single_calls

    def write_json(self, path):
        start_time = time.time()

        dict_data = self.__dict__
        filtered = self.filtered
        del dict_data["filtered"]
        numpy.save(path, dict_data)

        self.filtered = filtered

        info_path = path[:path.rfind(".")] + "_info.json"
        dataset_info = {
            "samples": [sample_def.id for sample_def in self.samples],
            "snps": [snp_def.allele_id for snp_def in self.snps],
            "sequences": [snp_def.sequence_ref for snp_def in self.snps],
            "pops": {sample_def.population: 0 for sample_def in self.samples}
        }

        for sample_def in self.samples:
            dataset_info["pops"][sample_def.population] += 1

        with open(info_path, "w") as info_out:
            json.dump(dataset_info, indent=4, fp=info_out)

        # test = numpy.asarray(self.__dict__)
        # numpy.savez(path, working_dir=self.working_dir, batch_id=self.batch_id, snps=self.snps, samples=self.samples,
        #             type=self.type, calls=self.calls, read_counts=self.read_counts, replicates=self.replicates,
        #             replicate_counts=self.replicate_counts)
        # with open(path, "w") as outfile:
        #     # jsonstr = json.dumps(self.__dict__, default=lambda o: o.__dict__)
        #     # outfile.write(jsonstr)
        #
        #     # json.dump(self.__dict__, default=lambda o: o.__dict__, fp=outfile)
        #
        #     outfile.write("{\n")
        #
        #     # Need custom JSON writing as dumping it all to a string then writing is a massive amount of memory.
        #     # json.dump directly to file is incredibly slow - not feasible to use.
        #
        #     first_row = True
        #     for key, val in self.__dict__.items():
        #         if key not in Dataset.data_cols and key != "filtered":
        #             outfile.write('{}\t"{}": {}'.format((",\n" if not first_row else ""), key, json.dumps(val, default=lambda o: o.__dict__)))
        #         elif key in Dataset.data_cols:
        #             # Need to custom write out the numpy as its incompatible with json.dumps
        #             outfile.write(',\n\t"' + key + '": {')
        #
        #             first_row = True
        #             for allele_id, snp_data in val.items():
        #                 list_data = snp_data.tolist()
        #                 if len(list_data) > 0 and type(list_data[0]) == numpy.ndarray:
        #                     list_data = [item.tolist() for item in list_data]
        #
        #                 outfile.write(
        #                     '{}"{}":{}'.format("," if not first_row else "", allele_id, json.dumps(list_data)))
        #                 first_row = False
        #
        #             outfile.write("}")
        #
        #         first_row = False
        #
        #     outfile.write("\n}")
        #     outfile.flush()

        log.info("Time to write dataset to file: {}".format(time.time() - start_time))

    @staticmethod
    def read_json(path):
        start_time = time.time()

        dataset = numpy.load(path).tolist()
        dataset = Dataset(**dataset)

        log.info("Dataset loaded - {} Samples, {} SNPs".format(len(dataset.samples), len(dataset.snps)))

        # file_data = numpy.load(path)
        # for key in file_data:
        #     # if key not in Dataset.data_cols:
        #     dict_val[key] = file_data[key].tolist()
        #

        # snps = []
        # for snp_dict in file_data["snps"]:
        #     snp = SNPDef(**snp_dict)
        #     snps.append(snp)
        # dict_val["snps"] = snps
        #
        # samples = []
        # for sample_dict in file_data["samples"]:
        #     if "pop" in sample_dict:
        #         del sample_dict["pop"]  # temp fix
        #
        #     sample = SampleDef(*sample_dict.values())
        #     samples.append(sample)
        # dict_val["samples"] = samples



        # with open(path) as dataset_file:
        #     line = dataset_file.readline()
        #
        #     # dict_val = json.load(dataset_file)
        #     dict_val = {}
        #
        #     # Read in incrementally to prevent memory overload
        #     # Python data types are all classes which is exceptionally inefficient - so incrementally read into numpy
        #
        #     while line != "":
        #         line = line.strip()
        #         if not re.match(r"^[\{\}]$", line):
        #             key, value = line.split(":", 1)
        #             line = None
        #
        #             key = re.sub(r"[\"\']", "", key)
        #
        #             if value[-1:] == ",":
        #                 value = value[:-1]
        #
        #             if key in Dataset.data_cols:
        #                 value = re.split(r"\]\],", value.strip()[1:-1])
        #
        #                 dict_val[key] = {}
        #                 for allele_data in value:
        #                     allele_id = allele_data[1:allele_data.find('"', 1)]
        #                     allele_data = allele_data[len(allele_id) + 3:]
        #
        #                     if allele_data[-2:] != "]]":
        #                         allele_data += "]]"
        #
        #                     dict_val[key] = numpy.asarray(json.loads(allele_data))
        #             else:
        #                 value = json.loads(value)
        #                 dict_val[key] = value
        #
        #         line = dataset_file.readline()
        #
        #     # dict_val = json.load(dataset_file)
        #


        # file_data = Dataset(file_data["type"], file_data["working_dir"], file_data["batch_id"], snps, samples,
        #                     file_data["calls"], file_data["read_counts"], file_data["replicates"],
        #                     file_data["replicate_counts"])

        # del dict_val["snps"]
        # del dict_val["samples"]
        # del dict_val["calls"]
        # del dict_val["read_counts"]
        #
        # for k, v in dict_val.items():
        #     setattr(dataset, k, v)

        # log.info("Time to load dataset from file: {:.2f}".format(time.time() - start_time))

        return dataset

    def __repr__(self):
        return self.batch_id + " @ " + self.working_dir
