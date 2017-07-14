import os
import json
import csv
import shutil
from subprocess import call

import pandas
import operator
import re

from DartModules import RedundancyModule
from dartqc.DartUtils import stamp


class DartFileValidator:
    """
    Class for validating data in the input files before pre-processing but after the file schemas are created
    (currently cd-hit-2d clusters and validates clone ID's)
    """

    def __init__(self, data, attributes, read_counts, ids_file_path):
        self.data = data
        self.attributes = attributes
        self.read_counts = read_counts
        self.ids_file_path = ids_file_path

        self.seq_vals = []
        self.seq_renames = {}

        self.tmp_path = os.path.join(self.attributes["out_path"], "tmp")
        self.cdhit_path = attributes["args"]["cdhit_path"]

        self.id_fasta_path = ""
        self.data_fasta_path = ""
        self.cluster_path = ""

    def do_validations(self):
        os.makedirs(self.tmp_path, exist_ok=True)

        self.create_ids_fasta()
        self.create_data_fasta()
        self.cluster()

        self.validate_sequences()
        self.write_seq_vals()
        self.rename_sequences()

        shutil.rmtree(self.tmp_path, ignore_errors=True)

    def create_data_fasta(self):
        fasta_writer = RedundancyModule(data=self.data, attributes=self.attributes)

        fasta_writer._write_fasta()

        self.data_fasta_path = os.path.join(self.tmp_path, self.attributes["project"] + "_Seqs.fasta")

    def create_ids_fasta(self, ids_file_path=None, fasta_file_path=None):
        # Use the classes file path if it isn't passed in now.
        if ids_file_path is None:
            ids_file_path = self.ids_file_path

        # Read in the official list of clone ID's and convert it into the fasta format
        fasta_str = ""
        with open(ids_file_path, "r") as ids_file:
            ids_reader = csv.reader(ids_file, delimiter=',')

            for row in ids_reader:
                if ids_reader.line_num > 1:
                    fasta_str = fasta_str + ">" + row[0] + "\n" + row[1][:60] + "\n" + (
                        row[1][60:] + "\n" if len(row[1]) > 60 else "")

        # Write the fasta file
        if fasta_file_path is None:
            fasta_file_path = os.path.abspath(
                os.path.join(self.tmp_path, self.attributes["project"] + "_ids.fasta"))

        self.id_fasta_path = fasta_file_path

        with open(fasta_file_path, "w") as fasta_file:
            fasta_file.write(fasta_str)
            fasta_file.flush()

    def cluster(self):
        """ Run CDHIT-EST for sequences, install with sudo apt install cd-hit on Ubuntu """

        # self.messages.get_cdhit_message(identity)

        cdhit_path = self.attributes["args"]["cdhit_path"]
        if cdhit_path is None:
            cdhit_path = "cd-hit-2d"

        if "cd-hit-est-2d" not in cdhit_path and "cd-hit-est" in cdhit_path:
            cdhit_path = cdhit_path.replace("cd-hit-est", "cd-hit-est-2d")

        file_name = self.attributes["project"] + "_clustered.fasta"

        out_file = os.path.join(self.tmp_path, file_name)

        with open(os.devnull, "w") as devnull:
            # -i    This flags for the set list you want, ie from additional csv file we will input
            # -i2   This flags for the fasta generated from the genotype file
            # -o    Produces the output file, default *.clstr, force to have extension, *.fasta
            # -c    Is the similarity threshold, set to 0.99 for ID matching purposes
            # -n    Word length (for matching algorithm), leave at 10
            # -d    SNP ID name length carried from input to output file, min can't be less than 25 otherwise SNP_ID become truncated
            # -M    Allocated memory
            # -T    Threads allocated
            call([cdhit_path, "-i", self.id_fasta_path, "-i2", self.data_fasta_path, "-o", out_file, "-c", "0.99", "-n",
                  "10", "-d", "20", "-M", "16000", "-T", "1"], stdout=devnull)

        self.cluster_path = out_file + ".clstr"

    def validate_sequences(self):
        clusters = []
        with open(self.cluster_path, "r") as clstr_file:
            lines = clstr_file.readlines()

            cluster = []
            for row in lines:
                if row.startswith(">Cluster"):
                    # If this is a new cluster, write add the previous cluster to the data array
                    if len(cluster) > 0:
                        clusters.append(cluster)
                        cluster = []
                else:
                    # Match and save cluster ID parts
                    regex_match = re.search("(?:[!>]*)(>[^-]*)(?:-+)(.*)(?:\.\.\..*)", row)
                    id = regex_match.group(1)
                    seq_pos = regex_match.group(2)

                    cluster.append((id, seq_pos))

            # Add the last cluster in the file (there isn't a following >Cluster line to promt the addition)
            clusters.append(cluster)

        self.seq_vals = []
        for index, cluster in enumerate(clusters):
            seq_list = []
            self.seq_renames = {"12949790": "62949790"}

            ref_seq = cluster[0]
            ref_seq_str = str(cluster[0][0]) + "--" + str(cluster[0][1])

            for seq in cluster[1:]:
                perfect_match = ref_seq == seq
                diff_id = ref_seq[0] != seq[0]
                diff_pos = ref_seq[1] != seq[1]

                seq_str = str(seq[0]) + "--" + str(seq[1])
                rename_str = None

                match_type = None
                if perfect_match:
                    match_type = "GOOD"
                elif not diff_pos and diff_id:
                    match_type = "BAD ID"
                    rename_str = str(ref_seq[0]) + "--" + str(seq[1])
                    self.seq_renames[seq[0]] = ref_seq[0]
                elif diff_pos and not diff_id:
                    match_type = "BAD LOC"
                elif not diff_pos and not diff_id:
                    match_type = "BAD ID & LOC"

                seq_list.append([match_type, seq_str, rename_str])

            self.seq_vals.append({
                "cluster_num": index,
                "ref_seq": ref_seq,
                "ref_seq_str": ref_seq_str,
                "sequences": seq_list
            })

    def write_seq_vals(self, out_file=None):
        # Add the CSV headers?
        output_csv = [["Cluster #", "Ref Seq", "Cluster Sequences..."]]

        for seq in self.seq_vals:
            row_data = [str(seq["cluster_num"]), seq["ref_seq_str"]] + [str(item) for sublist in seq["sequences"] for
                                                                        item in sublist]
            output_csv.append(row_data)

        if out_file is None:
            out_file = os.path.abspath(
                os.path.join(self.attributes["out_path"], self.attributes["project"] + "_seq_vals.csv"))
        with open(out_file, 'w') as vals_out:
            csv_writer = csv.writer(vals_out, delimiter=",", lineterminator='\n')
            csv_writer.writerows(output_csv)

        stamp("Sequence ID filtering info written to ", out_file)

    def rename_sequences(self):
        with open(self.attributes["args"]["raw_scheme"], "r") as infile:
            config = json.load(infile)
            counts_allele_col = config["allele_column"] - 1
            counts_clone_col = config["clone_column"] - 1
            counts_data_row = config["data_row"] - 1

        with open(self.attributes["args"]["call_scheme"], "r") as infile:
            config = json.load(infile)
            call_allele_col = config["allele_column"] - 1
            call_clone_col = config["clone_column"] - 1
            call_data_row = config["data_row"] - 1

        renamed_file_out = os.path.abspath(
            os.path.join(self.attributes["out_path"], self.attributes["project"] + "_data_validated.csv"))
        with open(self.attributes["args"]["call_file"], "r") as infile:
            csv_reader = csv.reader(infile)

            with open(renamed_file_out, "w") as val_out_file:
                csv_writer = csv.writer(val_out_file, delimiter=",", lineterminator='\n')

                for row in csv_reader:
                    if csv_reader.line_num > call_data_row:
                        current_clone_id = row[call_clone_col]
                        clone_id = current_clone_id if "|" not in current_clone_id else current_clone_id[
                                                                                        0:current_clone_id.index("|")]

                        if clone_id in self.seq_renames:
                            row[call_clone_col] = row[call_clone_col].replace(clone_id, self.seq_renames[clone_id])
                            row[call_allele_col] = row[call_allele_col].replace(clone_id, self.seq_renames[clone_id])

                    csv_writer.writerows([row])

        renamed_file_out = os.path.abspath(
            os.path.join(self.attributes["out_path"], self.attributes["project"] + "_read_counts_validated.csv"))
        with open(self.attributes["args"]["raw_file"], "r") as infile:
            csv_reader = csv.reader(infile)

            with open(renamed_file_out, "w") as val_out_file:
                csv_writer = csv.writer(val_out_file, delimiter=",", lineterminator='\n')

                for row in csv_reader:
                    if csv_reader.line_num > counts_data_row:
                        current_clone_id = row[counts_clone_col]
                        clone_id = current_clone_id if "|" not in current_clone_id else current_clone_id[
                                                                                        0:current_clone_id.index("|")]

                        if clone_id in self.seq_renames:
                            row[counts_clone_col] = row[counts_clone_col].replace(clone_id, self.seq_renames[clone_id])
                            row[counts_allele_col] = row[counts_allele_col].replace(clone_id,
                                                                                    self.seq_renames[clone_id])

                    csv_writer.writerows([row])

    def get_data(self):
        return self.data, self.attributes, self.read_counts
