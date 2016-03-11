__author__ = 'esteinig'

import os
import csv
import operator

from subprocess import call

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

class Tricoder:

    """
    Class for calculating statistics for genotypes and SNPs. Also contains a variety of helper functions for use
    in Dart.

    """

    def __init__(self):
        pass

    def calculate_call_rate(self, calls, sample_number, missing='-'):

        """
        Calculates call rate across samples for a single SNP. Pass a list with two lists containing calls for A1 and A2.

        """

        geno = list(zip(calls[0], calls[1]))
        return 1 - (geno.count((missing, missing))/sample_number)

    def calculate_maf(self, calls, sample_number, present='1', missing="-"):

        """
        Calculates minor allele frequency for a single SNP. Pass a list with two lists containing calls for A1 and A2.
        Remove missing calls from sample number, a bit slower with zipping the alleles first, but more robust than
        counting missing only in one allele and assuming the missing call occurs also on the other allele.

        """

        if sample_number != len(calls[0]) or sample_number != len(calls[1]):
            print("Warning: Number of samples does not correspond number of allele calls.")

        geno = list(zip(calls[0], calls[1]))
        adjusted_samples = sample_number - geno.count((missing, missing))

        freq_allele_one = calls[0].count(present) / adjusted_samples
        freq_allele_two = calls[1].count(present) / adjusted_samples

        return min(freq_allele_one, freq_allele_two)

    ### Helper Functions ###

    def find_between(self, s, first, last):

        """
        Finds the substring between two strings.

        """

        try:
            start = s.index(first) + len(first)
            end = s.index(last, start)
            return s[start:end]
        except ValueError:
            return ""

class DartReader:

    def __init__(self):

        self.project = "Monodon"
        self.data = {}                  # Holds intital unfiltered data from DArT

        self.identity_clusters = {}
        self.identity_snps = 0

        self.sample_name = {}
        self.sample_number = 0

        self._data_row = 6              # Row number (non-pythonic) in Excel Spreadsheet: Start of Sequences / Data
        self._sample_row = 5            # Row number (non-pythonic) in Excel Spreadsheet: Sample Identification
        self._call_column = 17          # Column number (pythonic) in Excel Spreadsheet: Start of SNP Calls

        self.statistics = {}

        self.tricoder = Tricoder()

    def parse_cdhit(self, file):

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
                    id = self.tricoder.find_between(row[1], ">", "...")
                    ids.append(id)

    def read_double_row(self, file):

        """"
        Read data file - requires format specific to Monodon. Will be incorporated into a DArT Reader class to parse
        DArTs and SNPs in various formats. Currently supports only SNPs with double row and encoding:

        AB: (1,1)
        AA: (1,0)
        BB: (0,1)
        NA: (-,-) ... a bit sleepy, little Owl?

        """

        with open(file, 'r') as data_file:
            reader = csv.reader(data_file)

            row_index = 1
            snp_count = 0
            allele_id = None
            allele_index = 1

            for row in reader:

                if row_index == self._sample_row:
                    sample_names = [r for r in row if r != '*']
                    self.sample_names = {k: v for (k, v) in enumerate(sample_names)}
                    self.statistics["samples.total"] = len(sample_names)
                    self.sample_number = len(sample_names)

                # Data Rows
                if row_index > self._data_row:

                    # Get reduced data by unique allele ID in double Rows (K: Allele ID, V: Data)
                    # Implement Error checks for conformity between both alleles: SNP Position, Number of Individuals
                    if allele_index == 1:

                        allele_id = row[0]

                        entry = {"allele_id": row[0], "clone_id": row[1], "allele_seq": row[2], "snp": row[3],
                                 "snp_position": int(row[4]), "average_read_count_reference": float(row[14]),
                                 "average_read_count_snp": float(row[15]), "average_replication": float(row[16]),
                                 "calls": [row[self._call_column:]]}  # Add allele 1

                        self.data[allele_id] = entry

                        snp_count += 1
                        allele_index = 2
                    else:
                        # Add allele 2 and calculate MAF of SNP
                        self.data[allele_id]["calls"].append(row[self._call_column:])

                        # Manual calculation    of MAF and Call Rate, much easier than parsing and fast enough,
                        # values correspond to pre-calculated ones from DArT
                        self.data[allele_id]["maf"] = self.tricoder.calculate_maf(self.data[allele_id]["calls"],
                                                                                  self.sample_number)
                        self.data[allele_id]["call_rate"] = self.tricoder.calculate_call_rate(self.data[allele_id]["calls"],
                                                                                              self.sample_number)

                        allele_index = 1

                row_index += 1

            self.statistics["snps.total"] = snp_count

class DartWriter:

    def __init__(self, dart_control):

        self.project = dart_control.project

        self.sample_names = dart_control.sample_names
        self.sample_number = dart_control.sample_number

        self.snps_total = dart_control.snps_total
        self.snps_filtered = dart_control.snps_filtered

        self.fasta_path = None

    def write_fasta(self, path=os.getcwd(), filtered=False):

        file_name = self.project + "_DArT_Seqs"

        if filtered:
            seqs = [SeqRecord(Seq(data["allele_seq"], IUPAC.unambiguous_dna), id=snp_id, name="", description="")
                    for snp_id, data in self.snps_filtered.items()]
            file_name += "_Filtered.fasta"
        else:
            seqs = [SeqRecord(Seq(data["allele_seq"], IUPAC.unambiguous_dna), id=snp_id, name="", description="")
                    for snp_id, data in self.snps_total.items()]
            file_name += "_Total.fasta"

        print("---------------------------------------------------------------------------------")
        print("Writing", len(seqs), "sequences to file (.fasta): ", file_name)
        print("---------------------------------------------------------------------------------")

        self.fasta_path = os.path.join(path, file_name)

        with open(self.fasta_path, "w") as fasta_file:
            SeqIO.write(seqs, fasta_file, "fasta")


class DartControl:

    """

    DArT Quality Control Pipeline

    Prototype


    """

    def __init__(self, dart):

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

        self.sample_names = dart.sample_names            # Number (Key) and sample ID (Value)
        self.sample_number = dart.sample_number          # Total number of samples

        self.filter_count = 1           # Number of applied filters in sequence, resets when filtering on total SNPs

        self.statistics = dart.statistics   # Dictionary of basic statistics for writing summary QC

        self.tricoder = Tricoder()      # Inititate Tricoder

        # Make temporary output directory
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
            for allele_id in v["allele_ids"]:
                self.snps_duplicate[allele_id] = self.snps_unique.pop(allele_id)

        self.statistics["clones.duplicate"] = len(self.clones_duplicate)
        self.statistics["snps.duplicate"] = len(self.snps_duplicate)
        self.statistics["snps.unique"] = len(self.snps_unique)

        print("---------------------------------------------------------------------------------")
        print("Number of duplicated clone IDs:", len(self.clones_duplicate))
        print("Number of duplicated SNPs:", len(self.snps_duplicate))
        print("Number of unique SNPs:", len(self.snps_unique))
        print("---------------------------------------------------------------------------------")

    def find_identity_clusters(self, identity=0.95, word_size=10, description_length=0, cdhit_path=None):

        """
        Clusters the reference allele sequences with CDHIT-EST and parses the clusters for selecting best sequences
        with selectors in second function.

        Check out if it is ok only to use the reference allele, not the one containing the SNP. Also, CD-HIT returns
        slightly different cluster configurations for each run, be wary and check it out.

        """

        dart_writer = DartWriter(self)
        dart_writer.write_fasta()

        if cdhit_path is None:
            cdhit_path = "cdhit-est"

        file_name = self.project + "_IdentityClusters_" + str(identity)

        out_file = os.path.join(self._tmp_path, file_name)
        cluster_file = os.path.join(self._tmp_path, file_name + '.clstr')

        call([cdhit_path, "-i", dart_writer.fasta_path, "-o", out_file, "-c", str(identity), "-n", str(word_size),
              "-d", str(description_length)])

        clstr_reader = DartReader()
        clstr_reader.parse_cdhit(cluster_file)

        self.cluster_duplicate = clstr_reader.identity_clusters

        self.cluster_snps = clstr_reader.identity_snps

        print("\n\nCLUSTERING SEQUENCES")
        print("---------------------------------------------------------------------------------")
        print("Generated CD-HIT identity clusters of sequences at a threshold of", str(identity*100) + "%.")
        print("Number of identity clusters:", len(self.cluster_duplicate))
        print("Number of duplicate identity SNPs:", self.cluster_snps)
        print("---------------------------------------------------------------------------------")

    def select_best_identity_seqs(self, selector="maf"):

        before = len(self.snps_total)

        error_msg = "You need to find identity clusters first, then remove sequences."

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

        print("REMOVING IDENTITY CLUSTERED SEQUENCES")
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

        print("Number of total SNPs after inclusion of best duplicate SNPs:", len(self.snps_total))

    def _compare_entries(self, ids, selector="maf"):

        """
        Gets data from dictionary for each duplicate SNP (MAF, Average Replication, Call Rate) according to 'selector'
        and returns the allele identification of the best entry.

        Later rank the data by multiple categories (e.g. if MAF > ... and Call Rate > ...).

        """

        entries_stats = [[i, self.data[i][selector]] for i in ids]
        entries_ranked = sorted(entries_stats, key=operator.itemgetter(1), reverse=True)

        return entries_ranked[0][0]

    def filter_snps(self, selector="maf", threshold=0.02, comparison="<", data="total"):

        if data == "total":
            before = len(self.snps_total)

            # Switch around comparison symbols, since we are getting the retained SNPs, but ask for Filter
            if comparison == "<":
                snps_retained = {k: v for (k, v) in self.snps_total.items() if v[selector] >= threshold}
            elif comparison == ">":
                snps_retained = {k: v for (k, v) in self.snps_total.items() if v[selector] <= threshold}
            else:
                snps_retained = {k: v for (k, v) in self.snps_total.items() if v[selector] == threshold}

            self.snps_filtered = snps_retained.copy()

            self.filter_count = 1

        elif data == "filtered" and self.snps_filtered:
            before = len(self.snps_filtered)

            if comparison == "<":
                snps_retained = {k: v for (k, v) in self.snps_filtered.items() if v[selector] >= threshold}
            elif comparison == ">":
                snps_retained = {k: v for (k, v) in self.snps_filtered.items() if v[selector] <= threshold}
            else:
                snps_retained = {k: v for (k, v) in self.snps_filtered.items() if v[selector] == threshold}

            self.snps_filtered = snps_retained.copy()

            # Keep track of applied filters to already filtered SNPs:
            self.filter_count += 1

        else:
            raise(ValueError("Data category '" + data + "' or previously filtered SNPs not available."))

        after = len(self.snps_filtered)

        print("---------------------------------------------------------------------------------")
        print("Filter", str(self.filter_count) + ":", selector.upper(), comparison, str(threshold), "on", data.upper())
        print("Removed", before - after, "SNPs from a total of", before, "SNPs, retaining", after, "SNPs.")
        print("---------------------------------------------------------------------------------")


