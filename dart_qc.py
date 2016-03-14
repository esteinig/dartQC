__author__ = 'esteinig'

import os
import csv
import time
import operator

from subprocess import call

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

class DartReader:

    """

    Class for reading data from DArT.

    """

    def __init__(self):

        self.project = "Monodon"

        # Parsing raw data from DArT

        self.raw_file = ''              # File name with raw data from DArT
        self.data = {}                  # Holds initial unfiltered data from DArT

        self.header = []                # Holds the lines before the actual header for statistics and data

        self.sample_names = []
        self.sample_number = 0

        # Row numbers (non-pythonic) in Excel Spreadsheet

        self._data_row = 7              # Start of Sequences / Data
        self._sample_row = 5            # Sample Identification

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
        self._average_pic = 14
        self._read_count_ref_column = 15
        self._read_count_snp_column = 16
        self._replication_average = 17
        self._call_column = 18
        self._sample_column = 18

        self.get_clone_id = False
        self._clone_split = '|'

        # Parsing CD HIT

        self.identity_clusters = {}
        self.identity_snps = 0

        # Initiating helper Tricoder

        self.tricoder = Tricoder()

    def parse_cdhit(self, file):

        """
        Parses the CDHIT cluster output and retains only cluster with more than one sequence for selection and
        removal from total SNPs.

        """

        with open(file, "r") as clstr_file:
            clstr_reader = csv.reader(clstr_file, delimiter=' ')

            cluster_id = 0
            ids = []

            for row in clstr_reader:

                if row[0] == ">Cluster":
                    if len(ids) > 1:
                        self.identity_clusters[cluster_id] = ids
                        self.identity_snps += len(ids)
                    ids = []
                    cluster_id += 1
                else:
                    allele_id = self.tricoder.find_between(row[1], ">", "...")
                    ids.append(allele_id)

    def read_double_row(self, file):

        """"

        AB: (1,1)
        AA: (1,0)
        BB: (0,1)
        NA: (-,-) ... a bit sleepy, little Owl?

        """

        self.raw_file = file

        with open(file, 'r') as data_file:
            reader = csv.reader(data_file)

            row_index = 1  # Non-pythonic for Excel Users
            snp_count = 0

            allele_id = None
            allele_index = 1

            for row in reader:

                if row_index <= self._data_row-2:  # Don't include description header
                    self.header.append(row)

                if row_index == self._sample_row:
                    self.sample_names = row[self._sample_column-1:]
                    self.sample_number = len(self.sample_names)

                # Data Rows
                if row_index >= self._data_row and any(row):

                    # Get reduced data by uniquesample_name allele ID in double Rows (K: Allele ID, V: Data)
                    # Implement Error checks for conformity between both alleles: SNP Position, Number of Individuals
                    if allele_index == 1:

                        allele_id = row[self._id-1]
                        clone_id = row[self._clone-1]

                        if self.get_clone_id:
                            clone_id = clone_id.split(self._clone_split)[0]

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
                                 "average_pic": row[self._average_pic-1],
                                 "average_read_count_ref": float(row[self._read_count_ref_column-1]),
                                 "average_read_count_snp": float(row[self._read_count_snp_column-1]),
                                 "average_replication": float(row[self._replication_average-1]),
                                 "calls": [row[self._call_column-1:]]}  # Add allele 1

                        self.data[allele_id] = entry

                        snp_count += 1
                        allele_index = 2
                    else:
                        # Add sequence and calls of second allele and calculate MAF of SNP
                        self.data[allele_id]["calls"].append(row[self._call_column-1:])
                        self.data[allele_id]["allele_seq_snp"] = row[self._seq-1]
                        self.data[allele_id]["snp"] = row[self._snp-1]

                        # Manual calculation of MAF and Call Rate, much easier than parsing and fast enough.
                        self.data[allele_id]["maf"] = self.tricoder.calculate_maf(self.data[allele_id]["calls"],
                                                                                  self.sample_number)
                        self.data[allele_id]["call_rate"] = self.tricoder.calculate_call_rate(self.data[allele_id]["calls"],
                                                                                              self.sample_number)

                        if len(self.data[allele_id]["calls"]) != 2:
                            raise(ValueError("Error. Genotype of", allele_id, "does not contain two alleles.",
                                  "Check if starting row for data is correctly specified."))

                        allele_index = 1

                row_index += 1

class DartWriter:

    """

    Class for writing output from DART QC.

    """

    def __init__(self, dart_control):

        self.qc = dart_control
        self.raw = dart_control.raw

        self.fasta_path = None

    def write_fasta(self, path=os.getcwd(), filtered=False):

        file_name = self.qc.project + "_DArT_Seqs"

        if filtered:
            seqs = [SeqRecord(Seq(data["allele_seq_ref"], IUPAC.unambiguous_dna), id=snp_id, name="", description="")
                    for snp_id, data in self.qc.snps_filtered.items()]
            file_name += "_Filtered.fasta"
        else:
            seqs = [SeqRecord(Seq(data["allele_seq_ref"], IUPAC.unambiguous_dna), id=snp_id, name="", description="")
                    for snp_id, data in self.qc.snps_total.items()]
            file_name += "_Total.fasta"

        print("\n---------------------------------------------------------------------------------")
        print("Writing", len(seqs), "sequences to file (.fasta): ", file_name)
        print("---------------------------------------------------------------------------------")

        self.fasta_path = os.path.join(path, file_name)

        with open(self.fasta_path, "w") as fasta_file:
            SeqIO.write(seqs, fasta_file, "fasta")

    def write_filtered(self, path=os.getcwd(), mode="dart"):

        """
        Write the filtered SNPs in various formats:

        1. Original: Reads the original and only retains alleles with allele ID in filtered SNPs.
        2. PLINK
        3. Structure

        """

        out_file = os.path.join(path, self.qc.project + "_dartQC_" + time.strftime("[%H:%M:%S]") + '.csv')

        if mode == "dart":

            head = [["AlleleID", "CloneID", "AlleleSequence", "SNP", "SnpPosition", "CallRate", "OneRatioRef", "OneRatioSnp",
                     "FreqHomRef", "FreqHomSnp", "FreqHets",	"PICRef", "PICSnp", "AvgPIC", "AvgCountRef", "AvgCountSnp",
                     "RepAvg"] + self.raw.sample_names]

            out_head = self.raw.header + head

            with open(out_file, 'w') as outfile:
                writer = csv.writer(outfile)

                writer.writerows(out_head)

                order = sorted(self.qc.snps_filtered.keys())

                for allele_id in order:
                    entry = self.qc.snps_filtered[allele_id]

                    allele_one = [allele_id, entry["clone_id"], entry["allele_seq_ref"], '-', entry["snp_position"],
                               entry["call_rate_dart"], entry["one_ratio_ref"], entry["one_ratio_snp"],
                               entry["freq_homozygous_ref"], entry["freq_homozygous_snp"], entry["freq_heterozygous"],
                               entry["pic_ref"], entry["pic_snp"], entry["average_pic"], entry["average_read_count_ref"],
                               entry["average_read_count_snp"], entry["average_replication"]]

                    for c in entry["calls"][0]:
                        allele_one.append(c)

                    allele_two = allele_one.copy()[:17]
                    allele_two[2] = entry["allele_seq_snp"]
                    allele_two[3] = entry["snp"]

                    for c in entry["calls"][1]:
                        allele_two.append(c)

                    writer.writerows([allele_one, allele_two])


class DartControl:

    """

    DArT Quality Control Pipeline

    Prototype

    """

    def __init__(self, dart):

        ### Data Storage ###

        self.raw = dart                 # Make the DataReader object accessible for DartWriter

        self.project = dart.project
        self.data = dart.data           # Holds intital unfiltered data from DArT

        self.clones_duplicate = {}      # Duplicate clones by clone ID (Key) and Dict: Number / List: SNP IDs (Value)
        self.cluster_duplicate = {}     # Identity clusters (Key) and List: SNP IDs (Value) SNPs

        self.cluster_snps = 0           # Number of SNPs in Identity Clusters

        self.snps_unique = {}           # SNPs unique after removing duplicate clone IDs
        self.snps_duplicate = {}        # SNPs with duplicate clone IDs after selecting best SNPs
        self.snps_clustered = {}        # SNPs in identity clusters and removed from total after selecting best SNPs

        self.snps_total = {}            # Final SNPs at each stage after removing duplicates and clustered SNPs
        self.snps_filtered = {}         # SNPs at each stage after processing through Filters

        self.filter_count = 1           # Number of applied filters in sequence, resets when filtering on total SNPs

        ### File IO ###

        self.raw_file = dart.raw_file

        self._data_row = dart._data_row
        self._id = dart._id
        self._clone = dart._clone

        ### Operations ###

        self.tricoder = Tricoder()

        self._tmp_path = os.path.join(os.getcwd(), "dart_qc_tmp")
        os.makedirs(self._tmp_path, exist_ok=True)

    def find_duplicate_clones(self):

        """
        Search for duplicate clone IDs in SNPs and

        """

        clone_counts = {}

        for k, v in self.data.items():

            clone_id = v["clone_id"]
            if clone_id not in clone_counts.keys():
                clone_counts[clone_id] = {"count": 1, "allele_ids": [k]}
            else:
                clone_counts[clone_id]["count"] += 1
                clone_counts[clone_id]["allele_ids"].append(k)

        self.clones_duplicate = {k: v for (k, v) in clone_counts.items() if v["count"] > 1}

        self.snps_unique = self.data.copy()

        for k, v in self.clones_duplicate.items():
            # Remove clones from unique SNPs, which are initially a copy of raw data (see above)
            for allele_id in v["allele_ids"]:
                self.snps_duplicate[allele_id] = self.snps_unique.pop(allele_id)

        print("\nSEARCHING FOR DUPLICATES")
        print("---------------------------------------------------------------------------------")
        print("Number of duplicated clone IDs:", len(self.clones_duplicate))
        print("Number of duplicated SNPs:", len(self.snps_duplicate))
        print("Number of unique SNPs:", len(self.snps_unique))
        print("---------------------------------------------------------------------------------")

    def find_identity_clusters(self, identity=0.95, word_size=10, description_length=0, cdhit_path=None):

        """
        Clusters the reference allele sequences with CDHIT-EST and parses the clusters for selecting and
        retaining best sequences.

        CD-HIT returns slightly different cluster configurations for each run due to greedy incremental algorithm,
        but little variation observed between runs in the data for P. monodon. Know thyself!

        """

        dart_writer = DartWriter(self)
        dart_writer.write_fasta()

        if cdhit_path is None:
            cdhit_path = "cdhit-est"

        file_name = self.project + "_IdentityClusters_" + str(identity)

        out_file = os.path.join(self._tmp_path, file_name)
        cluster_file = os.path.join(self._tmp_path, file_name + '.clstr')

        print("\nRUNNING CDHIT-EST")
        print("---------------------------------------------------------------------------------")
        call([cdhit_path, "-i", dart_writer.fasta_path, "-o", out_file, "-c", str(identity), "-n", str(word_size),
              "-d", str(description_length)])
        print("---------------------------------------------------------------------------------")

        clstr_reader = DartReader()
        clstr_reader.parse_cdhit(cluster_file)

        self.cluster_duplicate = clstr_reader.identity_clusters

        self.cluster_snps = clstr_reader.identity_snps

        print("\nCLUSTERING SEQUENCES BY IDENTITY")
        print("---------------------------------------------------------------------------------")
        print("Generated CD-HIT identity clusters of sequences at a threshold of", str(identity*100) + "%.")
        print("Number of identity clusters:", len(self.cluster_duplicate))
        print("Number of duplicate identity SNPs:", self.cluster_snps)
        print("---------------------------------------------------------------------------------")

    def select_best_identity_seqs(self, selector="maf"):

        """
        Selects the best sequences (SNPs) from the identity clusters (containing more than one sequence) according
        to 'selector', removes the best sequence from the cluster and removes bad sequences from the total
        SNPs, while retaining the best-pick SNP.

        """

        before = len(self.snps_total)

        error_msg = "You need to find identity clusters first, then remove    sequences."

        if not self.cluster_duplicate:
            raise(ValueError(error_msg))

        if not self.snps_total:
            raise(ValueError(error_msg))

        for cluster, ids in self.cluster_duplicate.items():
            # Find the best SNP and remove it form ID list, then remove from total SNPs
            best_snp = self._compare_entries(ids, selector=selector)
            ids.remove(best_snp)
            for i in ids:
                del self.snps_total[i]
                # Add the identification to the final removed clustered SNPs
                self.snps_clustered[i] = cluster

        print("\nRETAINING BEST CLUSTERED SEQUENCES")
        print("---------------------------------------------------------------------------------")
        print("Total number of SNPs before removing clustered SNPs:", before)
        print("Number of SNPs found in identity clusters:", self.cluster_snps)
        print("Number of clustered SNPs retained after selection with", selector.upper() + ":",
              len(self.cluster_duplicate))
        print("Number of clustered SNPs removed:", (self.cluster_snps - len(self.cluster_duplicate)))
        print("Total number of SNPs after removing clustered SNPs:", len(self.snps_total))
        print("---------------------------------------------------------------------------------")

    def select_best_clones(self, selector="maf"):

        """
        From each duplicate clone ID select the best according to entry in data dictionary (i.e. MAF,
        Average Replication or Call Rate) and add them to a copy of the unique SNPs for constructing
        the total unique SNPs.

        """

        self.snps_total = self.snps_unique.copy()

        for clone, clone_data in self.clones_duplicate.items():
            best_snp = self._compare_entries(clone_data["allele_ids"], selector=selector)
            # Popping them from the duplicate SNPs, so that leftovers are the actual duplicated ones,
            # excluding the best picks.
            self.snps_total[best_snp] = self.snps_duplicate.pop(best_snp)

        print("\nRETAINING BEST DUPLICATES")
        print("---------------------------------------------------------------------------------")
        print("Number of total SNPs after inclusion of best duplicate SNPs:", len(self.snps_total))
        print("---------------------------------------------------------------------------------")

    def _compare_entries(self, ids, selector="maf"):

        """
        Gets data from dictionary for each duplicate SNP (MAF, Average Replication, Call Rate) according to 'selector'
        and returns the allele identification of the best entry.

        Later rank the data by multiple categories (e.g. if MAF > ... and Call Rate > ...).

        """

        entries_stats = [[i, self.data[i][selector]] for i in ids]
        entries_ranked = sorted(entries_stats, key=operator.itemgetter(1), reverse=True)

        return entries_ranked[0][0]

    def filter_snps(self, selector="maf", threshold=0.02, comparison="<=", data="total"):

        """"
        Filter SNPs either on total SNPs or already filtered SNPs. Comparison sign is reversed because the list
        comprehensions retain SNPs, while the function asks dor removal of SNPs.

        """

        if comparison not in ["<=", ">=", "=="]:
            raise(ValueError("Comparison must be one of: <=, >=, =="))

        if data == "total":
            before = len(self.snps_total)

            # Switch around comparison symbols, since we are getting the retained SNPs, but ask for Filter
            if comparison == "<=":
                snps_retained = {k: v for (k, v) in self.snps_total.items() if v[selector] >= threshold}
            elif comparison == ">=":
                snps_retained = {k: v for (k, v) in self.snps_total.items() if v[selector] <= threshold}
            else:
                snps_retained = {k: v for (k, v) in self.snps_total.items() if v[selector] == threshold}

            self.snps_filtered = snps_retained.copy()

            self.filter_count = 1

        elif data == "filtered" and self.snps_filtered:
            before = len(self.snps_filtered)

            if comparison == "<=":
                snps_retained = {k: v for (k, v) in self.snps_filtered.items() if v[selector] >= threshold}
            elif comparison == ">=":
                snps_retained = {k: v for (k, v) in self.snps_filtered.items() if v[selector] <= threshold}
            else:
                snps_retained = {k: v for (k, v) in self.snps_filtered.items() if v[selector] == threshold}

            self.snps_filtered = snps_retained.copy()

            # Keep track of applied filters to already filtered SNPs:
            self.filter_count += 1

        else:
            raise(ValueError("Data category '" + data + "' or previously filtered SNPs not available."))

        after = len(self.snps_filtered)

        print("\nFILTERING SNPs")
        print("---------------------------------------------------------------------------------")
        print("Filter", str(self.filter_count) + ":", selector.upper(), comparison, str(threshold), "on", data.upper())
        print("Removed", before - after, "SNPs from a total of", before, "SNPs, retaining", after, "SNPs.")
        print("---------------------------------------------------------------------------------")

class Tricoder:

    """
    Class for calculating statistics for genotypes and SNPs. Also contains a variety of helper functions for use
    in Dart. Do not actually need to be in a class, but keeps the code nice and tidy.

    """

    def __init__(self):

        self.missing = '-'
        self.present = '1'
        self.absent = '0'

    def calculate_call_rate(self, calls, sample_number):

        """
        Calculates call rate across samples for a single SNP. Pass a list with two lists containing
        calls for A1 and A2.

        """

        geno = list(zip(calls[0], calls[1]))

        return 1 - (geno.count((self.missing, self.missing))/sample_number)

    def calculate_maf(self, calls, sample_number):

        """
        Calculates minor allele frequency for a single SNP. Pass a list with two lists containing calls for A1 and A2.
        Returns the minimum allele frequency for processing.

        """

        if sample_number != len(calls[0]) or sample_number != len(calls[1]):
            print("Warning: Number of samples does not correspond number of allele calls.")

        geno = list(zip(calls[0], calls[1]))

        adjusted_samples = sample_number - geno.count((self.missing, self.missing))

        het_count = geno.count((self.present, self.present))

        allele_one = geno.count((self.present, self.absent)) + (het_count/2)
        allele_two = geno.count((self.absent, self.present)) + (het_count/2)

        freq_allele_one = allele_one / adjusted_samples
        freq_allele_two = allele_two / adjusted_samples

        return min(freq_allele_one, freq_allele_two)

    ### Helper Functions ###

    def find_between(self, s, first, last):

        """Finds the substring between two strings."""

        try:
            start = s.index(first) + len(first)
            end = s.index(last, start)
            return s[start:end]
        except ValueError:
            return ""

