import csv
import json
import os

import logging
from collections import OrderedDict

import numpy
import pandas
import sys

import re

import time

from dartqc.Dataset import SNPDef, Dataset, SampleDef
from dartqc.PipelineOptions import Input
from dartqc.SimpleException import SimpleException

log = logging.getLogger(__file__)


class DartInput(Input):
    def get_name(self) -> str:
        return "dart"

    def get_description(self) -> str:
        return "DArT GBS data comes as 2 files: calls & read counts. " \
               "JSON mapping files are used to map columns/rows. " \
               "Note: Its easiest to run first without the JSON mapping files to generate them - then fix any problems."

    def read(self, working_dir: str, batch_id: str, files: [str], unknown_args: [] = None, **kwargs):
        # Dart reader takes potentially 4 files: calls, read counts & a json mapping for each

        # Work out which file is which based on file names
        calls_file, counts_file, calls_mapping_file, counts_mapping_file, pops_file = DartInput._identify_files(
            working_dir, files, silent="--generate_mappings" in unknown_args)

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

        # Read the mapping/scheme files in or auto-create new ones
        calls_mapping = None
        counts_mapping = None
        if calls_mapping_file is None and calls_file is not None:
            # Attempt to auto create mapping files if they don't exist
            log.info("Guessing call column/row mappings")
            calls_mapping = DartMapping.guess_mappings(calls_file)

            # Save as a template for editing externally
            calls_mapping.write_json(calls_file)
        elif calls_mapping_file is not None:
            calls_mapping = DartMapping.read_json(calls_mapping_file)

        if counts_mapping_file is None and counts_file is not None:
            # Attempt to auto create mapping files if they don't exist
            log.info("Guessing read count column/row mappings")
            counts_mapping = DartMapping.guess_mappings(counts_file)

            # Save as a template for editing externally
            counts_mapping.write_json(counts_file)
        elif counts_mapping_file is not None:
            counts_mapping = DartMapping.read_json(counts_mapping_file)

        # Custom Dart parameter to allow generation of mappings only
        if "--generate_mappings" in unknown_args:
            sys.exit(0)

        if len(files) < 2:
            log.error(
                "DaRT Reader requires at least 2 files: call data & read counts (you may/should also provide a json mapping file each as well)")
            sys.exit(1)

        start_time = time.time()

        # Read the CSV files in
        call_snp_defs, call_sample_defs, calls = DartInput._read_double_row_file(calls_file, calls_mapping)
        count_snp_defs, count_sample_defs, read_counts = DartInput._read_double_row_file(counts_file, counts_mapping,
                                                                                         numeric=True)

        log.info("Time to read files: {}\n".format(time.time() - start_time))

        # call_snps_ids = [snp.allele_id for snp in call_snp_defs]
        # count_snps_ids = [snp.allele_id for snp in count_snp_defs]

        # counts_missing_snps = [snp for snp in call_snps_ids if snp not in count_snps_ids]
        # call_missing_snps = [snp for snp in count_snps_ids if snp not in call_snps_ids]

        # Validate that samples & SNPs in the calls and read count files match
        call_sample_names = [sample.id for sample in call_sample_defs]
        count_sample_names = [sample.id for sample in count_sample_defs]

        counts_missing_snps = list(calls.keys() - read_counts.keys())
        call_missing_snps = list(read_counts.keys() - calls.keys())

        counts_missing_samples = [sample for sample in call_sample_names if sample not in count_sample_names]
        call_missing_samples = [sample for sample in count_sample_names if sample not in call_sample_names]

        log.info("Call vs Read Counts miss-match summary\n"
                 "{} SNPs in calls but not in counts: {}\n"
                 "{} SNPs in counts but not in calls: {}\n"
                 "{} Samples in calls but not in counts: {}\n"
                 "{} Samples in counts but not in calls: {}\n".format(
            len(counts_missing_snps), counts_missing_snps,
            len(call_missing_snps), call_missing_snps,
            len(counts_missing_samples), counts_missing_samples,
            len(call_missing_samples), call_missing_samples))

        # Delete missing SNPs from call & read count data
        for snp in counts_missing_snps:
            del calls[snp]
        for snp in call_missing_snps:
            del read_counts[snp]

        # Delete missing samples from call & read count data
        count_missing_sample_idxs = [idx for idx, sample_id in enumerate(count_sample_names) if
                                     sample_id in counts_missing_samples]
        for allele_id in calls:
            numpy.delete(calls[allele_id], count_missing_sample_idxs, axis=0)

        call_missing_sample_idxs = [idx for idx, sample_id in enumerate(call_sample_names) if
                                    sample_id in call_missing_samples]
        for allele_id in calls:
            numpy.delete(read_counts[allele_id], call_missing_sample_idxs, axis=0)

        # Get the union of SNP & sample definitions
        snp_defs = [snp for snp in call_snp_defs if snp.allele_id in read_counts.keys()]
        sample_defs = [sample for sample in call_sample_defs if sample.id in count_sample_names]

        log.info("Data contains {} samples and {} SNPs and {} calls\n".format(len(sample_defs), len(snp_defs),
                                                                              len(sample_defs) * len(snp_defs)))

        if len(sample_defs) == 0 or len(snp_defs) == 0:
            log.warning("There is no data (0 SNPs or samples) to filter - exiting")
            sys.exit(0)

        # Validate that all call values are one of: (-,-), (1,0), (0,1), (1,1)
        invalid_call_values = []
        # cmpl_count = 0
        for allele_id, snp_calls in calls.items():
            # if cmpl_count % 500 == 0:
            #     log.debug("Completed {} of {}".format(cmpl_count, len(calls)))
            # cmpl_count += 1

            valid_values = ["-", "0", "1"]
            valid_calls = numpy.in1d(snp_calls, valid_values)
            if not numpy.all(valid_calls):
                bad_idxs = set([numpy.math.ceil(idx / 2.0) for idx, good in enumerate(valid_calls) if not good])
                log.warning("Invalid call data for SNP {} samples {}".format(allele_id, [
                    "{} ({})".format(sample_defs[idx].id, snp_calls[idx]) for idx in bad_idxs]))

                for idx in bad_idxs:
                    calls[allele_id][idx] = ["-", "-"]

        # Collapse read count replicates (calls are already collapsed) and make sure they calls and counts have same order
        replicates = {}
        for idx, sample_id in enumerate(count_sample_names):
            if sample_id not in replicates:
                replicates[sample_id] = [idx]
            else:
                replicates[sample_id].append(idx)

        num_replicates = sum([len(idxs) - 1 for sample_id, idxs in replicates.items()])
        log.info("Collapsing read counts for {} ({:.0f}%) replicate sample"
                 .format(num_replicates, (num_replicates / len(count_sample_names)) * 100))

        collapsed_counts = DartInput._collapse_replicates(read_counts, replicates, snp_defs, sample_defs,
                                                          count_sample_names)

        log.info("Extracting replicate counts (stored separately from dataset.read_counts)")

        replicated_samples = [sample_id for sample_id, idxs in replicates.items() if len(idxs) > 1]
        replicate_counts = {
            allele_id: numpy.asarray(
                [counts[replicates[sample_id]] for sample_id in replicated_samples if len(replicates[sample_id]) > 0])
            for allele_id, counts in read_counts.items()}

        log.info("Finished reading data - creating dataset")

        if pops_file is not None and os.path.exists(pops_file):
            sample_defs = DartInput._read_pops_file(sample_defs, pops_file)

        # Create and return the dataset - this is the full representation of this genotype for both calls & read counts.
        return Dataset(self.get_name(), working_dir, batch_id, snp_defs, sample_defs, calls, collapsed_counts,
                       replicated_samples, replicate_counts)

    @staticmethod
    def _identify_files(working_dir: str, files: [str], silent: bool = False) -> [str]:
        files = sorted(files)
        pops_file = None

        for file in files:
            if file == "pops.csv" or file == "populations.csv":
                pops_file = file
                break

        mapping_files = [file for idx, file in enumerate(files) if
                         idx > 0 and os.path.basename(file).startswith(os.path.basename(files[idx - 1])[0:files[idx - 1].rfind(".")])]
        data_files = [file for file in files if file not in mapping_files and (pops_file is None or file != pops_file)]

        if len(data_files) != 2 and not silent:
            log.error("Possible error working out correct files:\n"
                      "\tMake sure scheme/mapping file names start with their related data files name (excluding extension)\n"
                      "\tName the population file pops.csv or populations.csv\n"
                      "\tMake sure the calls file name contains data or call\n"
                      "\tMake sure the read counts file contains count or depth\n\n")

        calls_file = None
        counts_file = None

        # Identify which file is the calls/data and wich is the read counts based on file names
        for idx in reversed(range(len(data_files))):
            file = data_files[idx]
            file_name = file[file.find("\\"):] if "\\" in file else file[file.find("/"):] if "/" in file else file

            if "data" in file_name or "call" in file_name:
                calls_file = file
                data_files.remove(file)
            elif "count" in file_name or "depth" in file_name:
                counts_file = file
                data_files.remove(file)

        if calls_file is None and len(data_files) > 0:
            calls_file = data_files.pop()

        if counts_file is None and len(data_files) > 0:
            counts_file = data_files.pop()

        if (counts_file is None or calls_file is None) and not silent:
            log.error("Missing calls {} and/or counts file {}\n"
                      "Either the files CMD line option is missing or the files are incorrectly named"
                      .format(calls_file, counts_file))

        if len(data_files) > 0:
            log.warning("Unknown/un-used input file(s): {}".format(data_files))

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
        if calls_file is not None and not os.path.isabs(calls_file):
            calls_file = os.path.join(working_dir, calls_file)
        if counts_file is not None and not os.path.isabs(counts_file):
            counts_file = os.path.join(working_dir, counts_file)
        if pops_file is not None and not os.path.isabs(pops_file):
            pops_file = os.path.join(working_dir, pops_file)
        if calls_mapping_file is not None and not os.path.isabs(calls_mapping_file):
            calls_mapping_file = os.path.join(working_dir, calls_mapping_file)
        if counts_mapping_file is not None and not os.path.isabs(counts_mapping_file):
            counts_mapping_file = os.path.join(working_dir, counts_mapping_file)

        return calls_file, counts_file, calls_mapping_file, counts_mapping_file, pops_file

    @staticmethod
    def _read_pops_file(sample_defs: [SampleDef], pops_file: str) -> [SampleDef]:
        with open(pops_file, "r") as file:
            line = file.readline()  # Ignore first line -> headers
            line = file.readline()
            while line is not None and line != "":
                try:
                    sample_id, pop = re.split(r"[,\t;]", line.strip())

                    found = False
                    for sample_def in sample_defs:
                        if sample_def.id == sample_id:
                            sample_def.population = pop
                            found = True

                    if not found:
                        log.warning("Sample in pop file doesn't exist in dataset: " + sample_id)

                except Exception as e:
                    log.warning(
                        "Invalid populations row: {} - each row should contain: <sample_id>,<population>".format(line))

                line = file.readline()

        return sample_defs

    @staticmethod
    def _collapse_replicates(counts, replicates, snps, sample_defs, counts_sample_names):
        start_time = time.time()

        numpy_matrix = numpy.asarray([counts[snp.allele_id] for snp in snps])

        numpy_matrix = numpy.stack(
            [numpy.sum(numpy_matrix[:, replicates[sample.id]], axis=1) for sample in sample_defs], axis=1)

        results = {snp.allele_id: numpy_matrix[idx] for idx, snp in enumerate(snps)}

        # results = {snp_def.allele_id: numpy.asarray([sum(counts[snp_def.allele_id][replicates[sample.id]]) for sample in sample_defs]) for snp_def in snps}
        # for snp_def in snps:
        #     results[snp_def.allele_id] = \
        #         numpy.asarray([sum(counts[snp_def.allele_id][replicates[sample_id]]) for sample_id in replicates])

        # for allele_id, snp_counts in counts.items():
        #     results[allele_id] = numpy.asarray([[sum([item for item in snp_counts[replicates[sample.id]]])] for sample in samples])

        # sample_names = [sample.id for sample in sample_defs]
        #
        # # non_rep = {k: v for k, v in replicates.items() if len(v) == 1}
        # reps = {k: v for k, v in replicates.items() if len(v) > 1}
        #
        # unique_idx_list = [counts_sample_names.index(name) for name in sample_names]
        #
        # results = {snp.allele_id: counts[snp.allele_id][unique_idx_list] for snp in snps}
        #
        # # log.info("Time: {}".format(time.time() - start_time))
        #
        # for sample_id, idxs in reps.items():
        #     for allele_id in results:
        #         results[allele_id][sample_names.index(sample_id)] = sum(counts[allele_id][idxs])
        #         # for idx in idxs[1:]:
        #         #     results[allele_id][sample_names.index(sample_id)] += counts[allele_id][idx]

        log.info("Time to collapse: {}\n".format(time.time() - start_time))

        return results

    @staticmethod
    def _read_double_row_file(file_path: str, mapping, numeric: bool = False) -> (
            [SNPDef], [SampleDef], {str: [numpy.ndarray]}):
        row_index = 0

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
                    try:
                        allele_id = row[mapping.allele_column]
                        clone_id = row[mapping.clone_column]
                        clone_id = clone_id.split("|")[0]  # Remove any extra crap (sometimes allele ID is given...)

                        if clone_id not in allele_id:
                            log.warning("Clone ID {} doesn't match allele ID {} on row {}".format(clone_id, allele_id,
                                                                                                  row_index))

                        if allele_id is None or clone_id is None or len(allele_id) == 0 or len(clone_id) == 0:
                            log.error("Row {} is missing allele ID or clone ID: {}".format(row_index, row))
                            was_error = True
                            continue

                        if data_row_idx % 2 == 0:
                            # First row for each SNP - read in full SNP details + call

                            all_headers = OrderedDict()
                            for idx, header in enumerate(headers):
                                all_headers[header] = row[idx]

                            snp_def = SNPDef(clone_id, allele_id, sequence_ref=row[mapping.sequence_column],
                                             rep_average=row[mapping.replication_column], all_headers=all_headers)

                            # Read the value
                            if numeric:
                                try:
                                    # data_row_1 = [int(val) for val in row[mapping.data_column:]]
                                    data_row_1 = (numpy.asarray(row[mapping.data_column:], dtype=int))
                                except:
                                    data_row_1 = [0 for sample in sample_names]
                                    log.error("Invalid data {}: {}".format(row_index, row[mapping.data_column:]))
                            else:
                                data_row_1 = row[mapping.data_column:]

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

                            snp_def.snp = snp_code_1
                            if snp_def.snp not in snp_def.allele_id:
                                log.warning(
                                    "SNP {} not in allele_ID {} for row {}".format(snp_def.snp, snp_def.allele_id,
                                                                                   row_index))
                            snp_def.all_headers["SNP"] = snp_def.snp

                            # Read in the second rows data

                            if numeric:
                                try:
                                    # data_row_2 = [int(val) for val in row[mapping.data_column:]]
                                    data_row_2 = (numpy.asarray(row[mapping.data_column:], dtype=int))
                                except:
                                    data_row_2 = [0 for sample in sample_names]
                                    log.error("Invalid data {}: {}".format(row_index, row[mapping.data_column:]))
                            else:
                                data_row_2 = row[mapping.data_column:]

                            snp_defs.append(snp_def)

                            # data[snp_def.allele_id] = list(numpy.array([data_row_1, data_row_2]).T)
                            if numeric:
                                data[snp_def.allele_id] = numpy.dstack((data_row_1, data_row_2))[0]
                            else:
                                data[snp_def.allele_id] = numpy.dstack((data_row_1, data_row_2))[0]

                            snp_def = None
                            data_row_1 = None
                    except Exception as e:
                        exc_type, exc_obj, exc_tb = sys.exc_info()
                        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                        log.error("({}:{}) {} - Exception reading {} row {}"
                                  .format(fname, exc_tb.tb_lineno, str(e), file_path, row_index))
                        was_error = True

                    data_row_idx += 1

                row_index += 1

                # Keep cleaning up memory so we don't run out.
                if data_row_idx > 0 and data_row_idx % 5000 == 0:
                    log.debug("Reading {} completed row {}".format(file_path, data_row_idx))

            if data_row_1 is not None or snp_def is not None:
                log.error("Last allele only has 1 row!  "
                          "Edit " + file_path + " to add the missing row or remove the single row for clone: "
                          + clone_id)
                was_error = True

            # Create the sample definitions
            if pops is not None and len(pops) != len(sample_names):
                log.error(
                    "Miss-matching sample names and populations - makes sure each sample has a pop and ID code set")
                was_error = True

            sample_defs = [SampleDef(sample_id, pops[idx] if pops is not None else "pop") for idx, sample_id in
                           enumerate(sample_names)]

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

    def __init__(self, sample_row=None, data_row=None, pop_row=None, clone_column=None, allele_column=None,
                 sequence_column=None, replication_column=None, data_column=None, **kwargs):
        """
        guess the specifications of input file and converting to input for DartReader

        limit: number of beginning rows to read that likely contain the meta data, usually no more than 30

        """

        self.sample_row = sample_row
        self.data_row = data_row
        self.pop_row = pop_row

        self.clone_column = clone_column
        self.allele_column = allele_column
        self.sequence_column = sequence_column
        self.replication_column = replication_column
        self.data_column = data_column
        # self.snp_column = snp_column
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

    def _zero_index(self):
        """ Reindex values for non-pythonic input to DartReader (better for users) """

        for key in self.__dict__:
            if getattr(self, key) is not None:
                setattr(self, key, getattr(self, key) - 1)

    def _excel_index(self):
        """ Reindex values for non-pythonic input to DartReader (better for users) """

        for key in self.__dict__:
            if getattr(self, key) is not None:
                setattr(self, key, getattr(self, key) + 1)

    @staticmethod
    def read_json(file_path: str):
        with open(file_path, "r") as file:
            mappings_dict = json.load(file)

            mappings = DartMapping(**mappings_dict)
            # for key in mappings_dict:
            #     setattr(mappings, key, mappings_dict[key])

            mappings._zero_index()

            return mappings

    def write_json(self, mapped_file_path):
        path, ext = os.path.splitext(mapped_file_path)
        out_path = path + "_scheme.json"

        log.info("Writing scheme to:" + out_path)

        self._excel_index()

        with open(out_path, "w") as outfile:
            json.dump(self.__dict__, sort_keys=True, fp=outfile, indent=4)

        self._zero_index()


DartInput()
