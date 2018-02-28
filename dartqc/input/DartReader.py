import csv
import json
import os

import logging
from collections import namedtuple

import numpy
import pandas
import sys

import re

import time

from Dataset import SNPDef, Dataset, SampleDef
from PipelineOptions import Input
from SimpleException import SimpleException

log = logging.getLogger(__file__)


class DartInput(Input):
    def get_name(self) -> str:
        return "dart"

    def get_description(self) -> str:
        return "DArT GBS data comes as 2 files: calls & read counts.\n" \
               "JSON mapping files are used to map columns/rows.\n" \
               "Note: Its easiest to run first without the JSON mapping files to generate them - then fix any problems."

    def read(self, working_dir: str, batch_id: str, files: [str], unknown_args: [] = None, **kwargs):
        # Dart reader takes potentially 4 files: calls, read counts & a json mapping for each

        if len(files) < 2:
            log.error(
                "DaRT Reader requires at least 2 files: call data & read counts (you may/should also provide a json mapping file each as well)")
            sys.exit(1)

        files = sorted(files)
        mapping_files = [file for idx, file in enumerate(files) if
                         idx > 0 and file.startswith(files[idx - 1][0:files[idx - 1].rfind(".")])]
        data_files = [file for file in files if file not in mapping_files]

        # Identify which file is the calls/data and wich is the read counts based on file names
        calls_file = data_files[0] if "data" in data_files[0] or "call" in data_files[0] \
                                      or "count" in data_files[1] or "depth" in data_files[1] else data_files[1]
        counts_file = files[1] if calls_file == data_files[0] else data_files[1]

        # Find the provided column/row mapping files (json)
        calls_mapping_file = mapping_files[0] if len(mapping_files) >= 1 and mapping_files[0].startswith(
            calls_file[0:calls_file.rfind(".")]) \
            else mapping_files[1] if len(mapping_files) >= 2 and mapping_files[1].startswith(
            calls_file[0: calls_file.rfind(".")]) else None
        counts_mapping_file = mapping_files[0] if len(mapping_files) >= 1 and mapping_files[0].startswith(
            counts_file[0:counts_file.rfind(".")]) \
            else mapping_files[1] if len(mapping_files) >= 2 and mapping_files[1].startswith(
            counts_file[0: counts_file.rfind(".")]) else None

        # Make relative paths work from the working dir rather than the present working directory (pwd)
        if not os.path.isabs(calls_file):
            calls_file = os.path.join(working_dir, calls_file)
        if not os.path.isabs(counts_file):
            counts_file = os.path.join(working_dir, counts_file)
        if calls_mapping_file is not None and not os.path.isabs(calls_mapping_file):
            calls_mapping_file = os.path.join(working_dir, calls_mapping_file)
        if counts_mapping_file is not None and not os.path.isabs(counts_mapping_file):
            counts_mapping_file = os.path.join(working_dir, counts_mapping_file)

        # Identify if this is excel instead of CSV (+get the sheet #)
        calls_excel_sheet = unknown_args[
            unknown_args.index("calls_sheet") + 1] if "calls_sheet" in unknown_args else None
        counts_excel_sheet = unknown_args[unknown_args.index("counts_sheet") + 1] if "counts_sheet" in unknown_args \
            else None

        # Convert files from excel to CSV
        if calls_excel_sheet is not None:
            output_path = calls_file[0: calls_file.rfind(".")] + ".csv"
            calls_file = DartInput._excel_to_csv(calls_file, output_path, calls_excel_sheet)

        if counts_excel_sheet is not None:
            output_path = counts_file[0: counts_file.rfind(".")] + ".csv"
            counts_file = DartInput._excel_to_csv(counts_file, output_path, counts_excel_sheet)

        calls_mapping = None
        counts_mapping = None
        if calls_mapping_file is None:
            # Attempt to auto create mapping files if they don't exist
            log.info("Guessing call column/row mappings")
            calls_mapping = DartMapping.guess_mappings(calls_file)

            # Save as a template for editing externally
            calls_mapping.write_json(calls_file)
        else:
            calls_mapping = DartMapping.read_json(calls_mapping_file)

        if counts_mapping_file is None:
            # Attempt to auto create mapping files if they don't exist
            log.info("Guessing read count column/row mappings")
            counts_mapping = DartMapping.guess_mappings(counts_file)

            # Save as a template for editing externally
            counts_mapping.write_json(counts_file)
        else:
            counts_mapping = DartMapping.read_json(counts_mapping_file)

        # Read the CSV files in
        call_snp_defs, call_sample_defs, calls = DartInput._read_double_row_file(calls_file, calls_mapping)
        count_snp_defs, count_sample_defs, read_counts = DartInput._read_double_row_file(counts_file, counts_mapping)

        call_snps_ids = [snp.allele_id for snp in call_snp_defs]
        call_sample_names = [sample.id for sample in call_sample_defs]

        count_snps_ids = [snp.allele_id for snp in count_snp_defs]
        count_sample_names = [sample.id for sample in count_sample_defs]

        counts_missing_snps = [snp for snp in call_snps_ids if snp not in count_snps_ids]
        call_missing_snps = [snp for snp in count_snps_ids if snp not in call_snps_ids]

        counts_missing_samples = [sample for sample in call_sample_names if sample not in count_sample_names]
        call_missing_samples = [sample for sample in count_sample_names if sample not in call_sample_names]

        log.info("Call vs Read Counts miss-match summary\n" +
                 "SNPs in calls but not in counts: " + str(counts_missing_snps) + "\n" +
                 "SNPs in counts but not in calls: " + str(call_missing_snps) + "\n" +
                 "Samples in calls but not in counts: " + str(counts_missing_samples) + "\n" +
                 "Samples in calls but not in counts: " + str(call_missing_samples))

        # Delete missing SNPs from call & read count data
        for snp in counts_missing_snps:
            del calls[snp]
        for snp in call_missing_snps:
            del read_counts[snp]

        # Delete missing samples from call & read count data
        for sample in reversed(counts_missing_samples):
            for allele_id in calls:
                del calls[allele_id][call_sample_names.index(sample)]
        for sample in reversed(call_missing_samples):
            for allele_id in read_counts:
                del read_counts[allele_id][count_sample_names.index(sample)]

        num_count_samples = len(count_sample_names) - len(call_missing_samples)
        num_call_samples = len(call_sample_names) - len(counts_missing_samples)
        num_reps = num_count_samples - num_call_samples
        log.info("There are {} ({:.0f}%) replicate samples\n".format(num_reps, num_reps / num_count_samples * 100))

        # Get the union of SNP & sample definitions
        snp_defs = [snp for snp in call_snp_defs if snp.allele_id in count_snps_ids]
        sample_defs = [sample for sample in call_sample_defs if sample.id in count_sample_names]

        # Validate that all call values are one of: (-,-), (1,0), (0,1), (1,1)
        invalid_call_values = []
        for allele_id, snp_calls in calls.items():
            for idx, call in enumerate(snp_calls):
                if call != Dataset.missing and call != Dataset.heterozygous and call != Dataset.homozygous_minor \
                        and call != Dataset.homozygous_major:
                    invalid_call_values.append("{}:{}".format(allele_id, sample_defs[idx].id))

        if len(invalid_call_values) > 0:
            log.warning("Invalid call data for: Is this a read counts/call file mis-match?"
                        "\n\t{}".format("\n\t".join(invalid_call_values)))

        # Validate that all read counts are are numeric
        invalid_count_values = []
        for allele_id, snp_counts in read_counts.items():
            for idx, counts in enumerate(snp_counts):
                if not re.match(r"\d+", counts[0]) or not re.match(r"\d+", counts[0]):
                    invalid_count_values.append("{}:{}".format(allele_id, sample_defs[idx].id))

        if len(invalid_count_values) > 0:
            log.warning("Invalid count data for: \n\t{}".format("\n\t".join(invalid_call_values)))

        # Collapse read count replicates (calls are already collapsed) and make sure they calls and counts have same order
        replicates = {}
        for idx, sample_id in enumerate(count_sample_names):
            if sample_id not in replicates:
                replicates[sample_id] = [idx]

        collapsed_counts = DartInput._collapse_replicates(read_counts, replicates, snp_defs, sample_defs)

        # Create and return the dataset - this is the full representation of this genotype for both calls & read counts.
        return Dataset(self.get_name(), working_dir, batch_id, snp_defs, sample_defs, calls, collapsed_counts)

    @staticmethod
    def _collapse_replicates(counts, replicates, snps, samples):
        # TODO:  I'm sure this could be faster using numpy etc. (currently ~30s for 100MB file)
        # start = time.time()

        results = {}
        for snp in snps:
            results[snp.allele_id] = []

            for idx, sample in enumerate(samples):
                val1 = 0
                val2 = 0

                for counts_idx in replicates[sample.id]:
                    val1 += int(counts[snp.allele_id][counts_idx][0])
                    val2 += int(counts[snp.allele_id][counts_idx][1])

                results[snp.allele_id].append((val1, val2))

        # log.debug("Time to collapse replicates: " + str(round((time.time() - start), 2)) + "s")

        return results

    @staticmethod
    def _read_double_row_file(file_path: str, mapping):
        row_index = 0

        numeric = False

        was_error = False

        with open(file_path, 'r') as file:
            reader = csv.reader(file)

            snp_count = 0

            allele_id = None
            allele_index = 1
            clone_id = None

            pops = []

            missing_call_regex = re.compile(r"^[- ]|(na)(n/a)$", re.IGNORECASE)

            snp_defs = []
            data = {}

            # 1st row data needed for interpreting second row
            snp_def = None
            data_row_1 = None

            sample_names = None
            pops = None
            headers = None

            data_row_idx = 0

            # header = []
            for row in reader:
                if row_index <= mapping.data_row - 1:  # Don't include description header
                    headers = row[:mapping.data_column]

                if row_index == mapping.sample_row:
                    sample_names = row[mapping.data_column:]

                elif row_index == mapping.pop_row:
                    pops = row[mapping.data_column:]

                # Get reduced data by unique allele ID in double Rows (K: Allele ID, V: Data)
                # Implement Error checks for conformity between both alleles: SNP Position, Number of Individuals

                # Data Rows, read from specified row and only if it contains data in at least one field (remove empties)
                elif row_index >= mapping.data_row and any(row):
                    allele_id = row[mapping.allele_column]
                    clone_id = row[mapping.clone_column]
                    clone_id = clone_id.split("|")[0]  # Remove any extra crap (sometimes allele ID is given...)

                    if allele_id is None or clone_id is None or len(allele_id) == 0 or len(clone_id) == 0:
                        log.error("Row is missing allele ID or clone ID\n"
                                  "Row index: " + str(row_index) +
                                  "\nContent: " + row)
                        was_error = True
                        continue

                    if data_row_idx % 2 == 0:
                        # First row for each SNP - read in full SNP details + call

                        all_headers = {header: row[idx] for idx, header in enumerate(headers)}
                        snp_def = SNPDef(clone_id, allele_id, sequence_ref=row[mapping.sequence_column],
                                         rep_average=row[mapping.rep_avg_column], all_headers=all_headers)

                        # Read the value
                        data_row_1 = row[mapping.data_column:]
                        if numeric:
                            data_row_1 = [0 if val is None or val is "" or re.match(missing_call_regex, val)
                                          else int(val) for val in data_row_1]

                    else:
                        # Second row for each SNP - add call & add this SNP to the list.

                        # Get the SNP code (eg. XX:C>A) and compare
                        # This checks that it isn't two clones of the same SNP
                        snp_code_1 = snp_def.allele_id[snp_def.allele_id.rfind("-") + 1:]
                        snp_code_2 = allele_id[allele_id.rfind("-") + 1:]

                        allele_clone_id_1 = snp_def.allele_id[: snp_def.allele_id.index("|")]
                        allele_clone_id_2 = allele_id[: allele_id.index("|")]

                        # Check that the two rows are for the same SNP!
                        if snp_def.clone_id != clone_id or allele_clone_id_1 != allele_clone_id_2 \
                                or snp_code_1 != snp_code_2:
                            raise SimpleException(
                                "Miss matched rows: " + snp_def.clone_id + "(" + snp_def.allele_id + ") -> "
                                + clone_id + "(" + allele_id + ").\n"
                                + "\t\t- Edit " + file_path + " so all allele's have 2 rows & are sequential\n"
                                + "\t\t- Check to make sure clone and allele IDs match correctly\n"
                                + "\t\t- Check the SNP matches (ending of the id such as 20:G>A)")

                        # Set SNP def values that are only available on the second row.
                        snp_def.sequence = row[mapping.sequence_column]
                        snp_def.snp = row[mapping.snp_column]

                        # Read in the second rows data
                        data_row_2 = row[mapping.data_column:]
                        if numeric:
                            data_row_2 = [0 if val is None or val is "" or re.match(missing_call_regex, val)
                                          else int(val) for val in data_row_2]

                        snp_defs.append(snp_def)
                        data[snp_def.allele_id] = list(zip(data_row_1, data_row_2))
                        snp_def = None
                        data_row_1 = None

                    data_row_idx += 1

                row_index += 1

            if data_row_1 is not None or snp_def is not None:
                log.error("Last allele only has 1 row!  "
                          "Edit " + file_path + " to add the missing row or remove the single row for clone: "
                          + clone_id)
                was_error = True

            # Create the sample definitions
            if len(pops) != len(sample_names):
                log.error(
                    "Miss-matching sample names and populations - makes sure each sample has a pop and ID code set")
                was_error = True

            sample_defs = [SampleDef(sample_id, pops[idx]) for idx, sample_id in enumerate(sample_names)]

            if was_error:
                raise SimpleException("Dart reader aborted due to errors")

            return snp_defs, sample_defs, data

    @staticmethod
    def _excel_to_csv(file_path, output_path, excel_sheet):
        log.info("Converting from Excel to CSV: " + file_path + "[" + excel_sheet + "]")

        data_xls = pandas.read_excel(file_path, excel_sheet, index_col=None)
        name, ext = os.path.splitext(os.path.basename(file_path))
        outfile = os.path.join(output_path, name + ".csv")

        data_xls.to_csv(outfile, encoding='utf-8', index=False)

        return outfile


class DartMapping:
    CONFIG = {
        "clone_column": [
            "CloneID"
        ],
        "allele_column": [
            "AlleleID"
        ],
        "sequence_column": [
            "AlleleSequence"
        ],
        "rep_avg_column": [
            "RepAvg"
        ],
        "snp_column": [
            "SNP"
        ]
    }

    def __init__(self):
        """
        guess the specifications of input file and converting to input for DartReader

        limit: number of beginning rows to read that likely contain the meta data, usually no more than 30

        """

        self.sample_row = None
        self.data_row = None
        self.pop_row = None

        self.clone_column = None
        self.allele_column = None
        self.sequence_column = None
        self.rep_avg_column = None
        self.data_column = None
        self.snp_column = None
        # self._reindex()

    @staticmethod
    def guess_mappings(file_path):
        mappings = DartMapping()

        """
        Guess the row indices for samples and data assuming:
            - header row begins after rows starting with "*"
            - data row starts after header row
            - samples are specified in the header row (above calls or raw counts)
        """

        top = pandas.read_csv(file_path, header=None, nrows=30)
        header_row = []

        for i, row in top.iterrows():
            if row[0] != "*":
                header_row = row
                mappings.sample_row = int(i)
                mappings.pop_row = int(i - 1)
                mappings.data_row = int(i + 1)
                break

        """
       Guess column indices for key words defined in the config file (excel_scheme.json), also get
       column indices for start of data assuming:
           - data starts in first column after placeholder ("*")
           - guessed from the row above the header (sample_row-1)
       """
        i = 0
        for col in top:
            data = top[col]
            if data[mappings.sample_row - 1] != "*":
                mappings.data_column = i
                break
            i += 1

        has_error = False
        for column, aliases in DartMapping.CONFIG.items():
            column_indices = [pandas.Index(header_row).get_loc(header) for header in aliases]

            if len(column_indices) < 1:
                has_error = True
                log.error("Could not find column header", column, "in dataset header. "
                                                                  "Please provide a different alias in dartQC/schemes/excel_scheme.json.")
            elif len(column_indices) >= 2:
                has_error = True
                log.error("Found more than one column header out of:", aliases, "make sure column "
                                                                                "headers are not repeated in data. You can change aliases that this command is "
                                                                                "looking for in dartQC/schemes/excel_scheme.json")
            else:
                setattr(mappings, column, column_indices[0])

        if has_error:
            log.error("Dart reader aborted due to mapping errors")
            sys.exit(1)

        return mappings

    # def _reindex(self):
    #     """ Reindex values for non-pythonic input to DartReader (better for users) """
    #
    #     scheme = {k: int(v + 1) for k, v in scheme.items()}
    #
    #     for k, v in sorted(scheme.items(), key=operator.itemgetter(1)):
    #         log.info(k, "=", v)
    #
    #     log.info("Please check these values in your data to ensure correct input for DartQC.")

    @staticmethod
    def read_json(file_path: str):
        with open(file_path, "r") as file:
            mappings_dict = json.load(file)
            mappings = DartMapping()
            for key in mappings_dict:
                setattr(mappings, key, mappings_dict[key])

            return mappings

    def write_json(self, mapped_file_path):
        path, ext = os.path.splitext(mapped_file_path)
        out_path = path + "_scheme.json"

        log.info("Writing scheme to:" + out_path)

        with open(out_path, "w") as outfile:
            json.dump(self.__dict__, sort_keys=True, fp=outfile, indent=4)


DartInput()
