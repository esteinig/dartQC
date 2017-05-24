import os
import csv
import json


class DartReader:

    """Class for reading raw calls."""

    def __init__(self):

        self.project = "Monodon"
        self.verbose = True

        self.out_path = ""

        # Parsing raw data

        self.raw_file = ''              # File name with raw data
        self.data = {}                  # Holds initial unfiltered data
        self.header = []                # Holds the lines before the actual header for statistics and data

        self.sample_names = []
        self.sample_size = 0
        self.snp_number = 0

        self.format = "double"

        # Row numbers (non-pythonic) in Excel Spreadsheet

        self._data_row = 7              # Start of Sequences / Data
        self._sample_row = 5            # Sample Identification
        self._pop_row = 0               # Sample Populations

        # Column numbers (non-pythonic) in Excel Spreadsheet

        self._id = 1
        self._clone = 2
        self._seq = 3
        self._snp = 4
        self._snp_position = 5
        self._call_rate_dart = 6
        self._one_ratio_ref = 7
        self._one_ratio_snp = 8
        self._freq_homozygous_ref = 9
        self._freq_homozygous_snp = 10
        self._freq_heterozygous = 11
        self._pic_ref = 12
        self._pic_snp = 13
        self._pic_average = 14
        self._read_count_ref = 15
        self._read_count_snp = 16
        self._rep_average = 17
        self._call = 18
        self._sample_column = 18

        # Meta Data by Individuals

        self.meta = {}

        self._id_meta = 1
        self._pop_meta = 2

        self.split_clone = False
        self.clone_split = "|"

        # Encoding Scheme

        self.homozygous_major = ("1", "0")
        self.homozygous_minor = ("0", "1")
        self.heterozygous = ("1", "1")
        self.missing = ("-", "-")

        self._dart_qc_encoding = {self.missing: "-", self.heterozygous: "0", self.homozygous_minor: "1",
                                  self.homozygous_major: "2"}

    def set_options(self, project="DartQC", verbose=True, out_path=os.getcwd(), homozygous_major=("1", "0"),
                    homozygous_minor=("0", "1"), heterozygous=("1", "1"), missing=("-", "-"), pop_row=0, sample_row=5,
                    data_start_row=7, id_col=1, clone_col=2, seq_col=3, snp_col=4, snp_position_col=5, call_rate_col=6,
                    one_ratio_ref_col=7, one_ratio_snp_col=8, freq_homozygous_ref_col=9, freq_homozygous_snp_col=10,
                    freq_heterozygous_col=11, pic_ref_col=12, pic_snp_col=13, pic_average=14, read_count_ref_col=15,
                    read_count_snp_col=16, rep_average_col=17, call_start_col=18, sample_start_col=18, clone_split="|",
                    split_clone=False, scheme=None):

        self.project = project
        self.verbose = verbose

        self.out_path = out_path

        self._data_row = data_start_row
        self._sample_row = sample_row
        self._pop_row = pop_row

        self.homozygous_major = homozygous_major
        self.homozygous_minor = homozygous_minor
        self.heterozygous = heterozygous
        self.missing = missing

        self._id = id_col
        self._clone = clone_col
        self._seq = seq_col
        self._snp = snp_col
        self._snp_position = snp_position_col
        self._call_rate_dart = call_rate_col
        self._one_ratio_ref = one_ratio_ref_col
        self._one_ratio_snp = one_ratio_snp_col
        self._freq_homozygous_ref = freq_homozygous_ref_col
        self._freq_homozygous_snp = freq_homozygous_snp_col
        self._freq_heterozygous = freq_heterozygous_col
        self._pic_ref = pic_ref_col
        self._pic_snp = pic_snp_col
        self._pic_average = pic_average
        self._read_count_ref = read_count_ref_col
        self._read_count_snp = read_count_snp_col
        self._rep_average = rep_average_col
        self._call = call_start_col
        self._sample_column = sample_start_col

        self.split_clone = split_clone
        self.clone_split = clone_split

        self._dart_qc_encoding = {self.missing: "-", self.heterozygous: "0", self.homozygous_minor: "1",
                                  self.homozygous_major: "2"}

        if scheme is not None:
            # Only basic supported for now.
            with open(scheme, "r") as infile:
                config = json.load(infile)

                self._id = config["allele_column"]
                self._clone = config["clone_column"]
                self._seq = config["sequence_column"]
                self._rep_average = config["replication_column"]
                self._call = config["data_column"]
                self._sample_column = config["data_column"]
                self._data_row = config["data_row"]
                self._sample_row = config["sample_row"]

    def read_pops(self, file, sep=','):

        """ Read file with header and two columns: 1 - ID, 2 - Population. ID must be the same as in Data. """

        meta_head = []

        with open(file, 'r') as infile:
            reader = csv.reader(infile, delimiter=sep)
            for row in reader:
                if meta_head:
                    self.meta[row[self._id_meta-1]] = row[self._pop_meta-1]
                else:
                    meta_head = row

    def read_double_row(self, file, encode=True, split_char="|", numeric=False, basic=True):

        """"Read data in double row format"""

        self.raw_file = file
        self.format = "double"
        self.clone_split = split_char

        with open(file, 'r') as data_file:
            reader = csv.reader(data_file)

            row_index = 1  # Non-pythonic numbering for Excel Users
            snp_count = 0

            allele_id = None
            allele_index = 1

            pops = []

            for row in reader:

                if row_index <= self._data_row-2:  # Don't include description header
                    self.header.append(row)

                if row_index == self._sample_row:
                    self.sample_names = row[self._sample_column-1:]
                    self.sample_size = len(self.sample_names)

                if row_index == self._pop_row:
                    pops = row[self._sample_column-1:]

                # Data Rows, read from specified row and only if it contains data in at least one field (remove empties)
                if row_index >= self._data_row and any(row):

                    # Get reduced data by unique allele ID in double Rows (K: Allele ID, V: Data)
                    # Implement Error checks for conformity between both alleles: SNP Position, Number of Individuals
                    if allele_index == 1:

                        allele_id = row[self._id-1]
                        clone_id = row[self._clone-1]

                        call_1 = row[self._call-1:]

                        if numeric:
                            call_1 = [int(call) for call in call_1]

                        if self.split_clone:
                            clone_id = clone_id.split(self.clone_split)[0]

                        if basic:
                            entry = {"allele_id": allele_id,
                                     "clone_id": clone_id,
                                     "allele_seq_ref": row[self._seq-1],
                                     "rep_average": float(row[self._rep_average-1]),
                                     "calls": [call_1]}  # Add allele calls 1
                        else:

                            entry = {"allele_id": allele_id,
                                     "clone_id": clone_id,
                                     "allele_seq_ref": row[self._seq-1],
                                     "snp_position": row[self._snp_position-1],
                                     "call_rate_dart": row[self._call_rate_dart-1],
                                     "one_ratio_ref": row[self._one_ratio_ref-1],
                                     "one_ratio_snp": row[self._one_ratio_snp-1],
                                     "freq_homozygous_ref": row[self._freq_homozygous_ref-1],
                                     "freq_homozygous_snp": row[self._freq_homozygous_snp-1],
                                     "freq_heterozygous": row[self._freq_heterozygous-1],
                                     "pic_ref": row[self._pic_ref-1],
                                     "pic_snp": row[self._pic_snp-1],
                                     "pic_average": row[self._pic_average-1],
                                     "read_count_ref": float(row[self._read_count_ref-1]),
                                     "read_count_snp": float(row[self._read_count_snp-1]),
                                     "rep_average": float(row[self._rep_average-1]),
                                     "calls": [call_1]}  # Add allele calls 1

                        self.data[allele_id] = entry

                        snp_count += 1
                        allele_index = 2
                    else:
                        # Add sequence and allele calls 2:

                        call_2 = row[self._call-1:]

                        if numeric:
                            call_2 = [int(call) for call in call_2]

                        self.data[allele_id]["calls"].append(call_2)
                        self.data[allele_id]["allele_seq_snp"] = row[self._seq-1]
                        self.data[allele_id]["snp"] = row[self._snp-1]

                        if len(self.data[allele_id]["calls"]) != 2:
                            raise(ValueError("Error. Genotype of", allele_id, "does not contain two alleles.",
                                  "Check if starting row for data is correctly specified."))
                        if encode:
                            self.data[allele_id]["calls"] = self._encode_dart(self.data[allele_id]["calls"])
                        else:
                            self.data[allele_id]["calls"] = list(zip(self.data[allele_id]["calls"][0],
                                                                     self.data[allele_id]["calls"][1]))

                        self.snp_number += 1
                        allele_index = 1

                row_index += 1

        # Check if reader picked up populations, if not generate generic names:
        if not pops:
            pops = ["Pop" for i in range(self.sample_size)]

        self.meta = dict(zip(self.sample_names, pops))

    def read_single_row(self, file):

        self.raw_file = file
        self.format = "single"

        pass

    def read_json(self, data_file, attribute_file):

        with open(data_file) as data_in:
            data = json.load(data_in)

        with open(attribute_file) as attr_in:
            attributes = json.load(attr_in)

        return data, attributes

    def get_data(self):

        attributes = {

            "project": self.project,
            "sample_size": self.sample_size,
            "sample_names": self.sample_names,
            "pops": self.meta,
            "missing": self._dart_qc_encoding[self.missing],
            "heterozygous": self._dart_qc_encoding[self.heterozygous],
            "homozygous_minor": self._dart_qc_encoding[self.homozygous_minor],
            "homozygous_major": self._dart_qc_encoding[self.homozygous_major],
            "out_path": self.out_path,
            "file": self.raw_file,
            "snps": len(self.data),
            "modules": {}

        }

        return self.data, attributes

    def _encode_dart(self, calls):

        if self.format == "double":
            calls = zip(calls[0], calls[1])

        return [self._dart_qc_encoding[snp_call] for snp_call in calls]