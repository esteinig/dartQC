__author__ = 'esteinig'

import os
import csv
import time
import numpy
import argparse
import operator
import textwrap

from subprocess import call

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

class CommandLine:

    def __init__(self):

        self.parser = argparse.ArgumentParser(description='DArT Quality Control', add_help=True)
        self.setParser()
        self.args = self.parser.parse_args()

        self.arg_dict = vars(self.args)

    def setParser(self):

        data_type = self.parser.add_mutually_exclusive_group(required=True)

        ### Required Options Input ###

        data_type.add_argument('-f', "--file", dest='data_file', required=True, type=str,
                                help="Name of input file from DArT")
        data_type.add_argument('-c', "--config", dest='config_file', required=False, type=str,
                                help="Name of configuration file for DartQC")

        ### Output ###

        self.parser.add_argument('-o', dest='output_format', default='plink', required=False, type=str,
                                 help="Output format, one of: dart, plink, structure")

        ### Population Designations for Output ###

        self.parser.add_argument('-p', '--pop', dest='pop_file', required=False, type=str,
                                 help="Name of file containing IDs and Populations")

        ### Encoding Schemes ###

        self.parser.add_argument('--major', dest='homozygous_major', default=("1", "0"), required=False, type=tuple,
                                 help="Homozygous major call encoding in DArT")
        self.parser.add_argument('--minor', dest='homozygous_minor', default=("0", "1"), required=False, type=tuple,
                                 help="Homozygous minor call encoding in DArT")
        self.parser.add_argument('--hetero', dest='heterozygous', default=("1", "1"), required=False, type=tuple,
                                 help="Heterozygous call encoding in DArT")
        self.parser.add_argument('--missing', dest='homozygous_minor', default=("-", "-"), required=False, type=tuple,
                                 help="Homozygous minor encoding in DArT")



class Tricoder:

    """
    Class for calculating statistics for genotypes and SNPs. Also contains a variety of helper functions for use
    in DartReaders and Writers.

    MAF
    Call Rate
    SNP Decoder / Encoder
    Sub-string Finder

    """

    def __init__(self):

        self.homozygous_major = ("1", "0")
        self.homozygous_minor = ("0", "1")
        self.heterozygous = ("1", "1")
        self.missing = ("-", "-")

        self.encoding_scheme = {}

    ### Calculations ###

    def calculate_call_rate(self, calls, sample_number):

        """
        Calculates call rate across samples for a single SNP. Pass a list with two lists containing
        calls for A1 and A2.

        """

        geno = list(zip(calls[0], calls[1]))

        return 1 - (geno.count(self.missing)/sample_number)

    def calculate_maf(self, calls, sample_number):

        """
        Calculates minor allele frequency for a single SNP. Pass a list with two lists containing calls for A1 and A2.
        Returns the minimum allele frequency for processing.

        """

        if sample_number != len(calls[0]) or sample_number != len(calls[1]):
            print("Warning: Number of samples does not correspond number of allele calls.")

        geno = list(zip(calls[0], calls[1]))

        adjusted_samples = sample_number - geno.count(self.missing)

        het_count = geno.count(self.heterozygous)

        allele_one = geno.count(self.homozygous_major) + (het_count/2)
        allele_two = geno.count(self.homozygous_minor) + (het_count/2)

        freq_allele_one = allele_one / adjusted_samples
        freq_allele_two = allele_two / adjusted_samples

        return min(freq_allele_one, freq_allele_two)

    ### SNP Operations ###

    def decode(self, snp_data):

        """ Decodes DArT encoding to specified format, requires list of zipped, one-row calls for SNPs"""

        decoded_data = []
        for snp in snp_data:
            new_snp = []
            for c in snp:
                try:
                    new_snp.append(self.encoding_scheme[c])
                except KeyError:
                    print("Could not find appropriate encoding scheme for:", c)

            decoded_data.append(new_snp)

        return decoded_data

    def get_encoding(self, dart_reader):

        """
        Gets the encoding scheme for homozygotes and heterozygotes when initiated in function for
        calculating MAF and Call Rate.

        """

        self.homozygous_major = dart_reader.homozygous_major
        self.homozygous_minor = dart_reader.homozygous_minor
        self.heterozygous = dart_reader.heterozygous
        self.missing = dart_reader.missing

    def get_encoding_scheme(self, dart_writer):

        self.encoding_scheme = {dart_writer.from_homozygous_major: dart_writer.to_homozygous_major,
                                dart_writer.from_homozygous_minor: dart_writer.to_homozygous_minor,
                                dart_writer.from_heterozygous: dart_writer.to_heterozygous,
                                dart_writer.from_missing: dart_writer.to_missing}

    ### Helper Functions ###

    def find_between(self, s, first, last):

        """Finds the substring between two strings."""

        try:
            start = s.index(first) + len(first)
            end = s.index(last, start)
            return s[start:end]
        except ValueError:
            return ""

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
        self.snp_number = 0

        # Row numbers (non-pythonic) in Excel Spreadsheet

        self._data_row = 7              # Start of Sequences / Data
        self._sample_row = 5            # Sample Identification
        self._pop_row = None

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

        # Meta Data by Individuals

        self.meta = {}
        self.meta_head = []

        self.pops = []                                      # IMPLEMENT POP READER

        self._id_meta = 1
        self._pop_meta = 2

        # Encoding Scheme

        self.homozygous_major = ("1", "0")
        self.homozygous_minor = ("0", "1")
        self.heterozygous = ("1", "1")
        self.missing = ("-", "-")

        # Parsing CD HIT

        self.identity_clusters = {}
        self.identity_snps = 0

        # Initiating helper Tricoder and setting

    def parse_cdhit(self, file):

        """
        Parses the CDHIT cluster output and retains only cluster with more than one sequence for selection and
        removal from total SNPs.

        """

        tricoder = Tricoder()

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
                    allele_id = tricoder.find_between(row[1], ">", "...")
                    ids.append(allele_id)

    def read_pops(self, file, sep=','):

        """ Read file with header and two columns: 1 - ID, 2 - Population. ID must be the same as in DArT Data. """

        with open(file, 'r') as infile:
            reader = csv.reader(infile, delimiter=sep)
            for row in reader:
                if self.meta_head:
                    self.meta[row[self._id_meta-1]] = row[self._pop_meta-1]
                else:
                    self.meta_head = row

    def read_double_row(self, file):

        """" Read DArT data in double row format. """

        self.raw_file = file

        tricoder = Tricoder()
        tricoder.get_encoding(self)

        with open(file, 'r') as data_file:
            reader = csv.reader(data_file)

            row_index = 1  # Non-pythonic for Excel Users
            snp_count = 0

            allele_id = None
            allele_index = 1

            for row in reader:

                if row_index <= self._data_row-2:  # Don't include description header
                    self.header.append(row)

                if self._pop_row is not None and row_index == self._pop_row:  # Add populations if in DArT (for Shannon)
                    self.pops = row[self._sample_column-1:]

                if row_index == self._sample_row:
                    self.sample_names = row[self._sample_column-1:]
                    self.sample_number = len(self.sample_names)

                # Data Rows, read from specified row and only if it contains data in fields
                if row_index >= self._data_row and any(row):

                    # Get reduced data by unique allele ID in double Rows (K: Allele ID, V: Data)
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
                        self.data[allele_id]["maf"] = tricoder.calculate_maf(self.data[allele_id]["calls"],
                                                                             self.sample_number)
                        self.data[allele_id]["call_rate"] = tricoder.calculate_call_rate(self.data[allele_id]["calls"],
                                                                                         self.sample_number)

                        if len(self.data[allele_id]["calls"]) != 2:
                            raise(ValueError("Error. Genotype of", allele_id, "does not contain two alleles.",
                                  "Check if starting row for data is correctly specified."))

                        self.snp_number += 1

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

        self.from_homozygous_major = dart_control.raw.homozygous_major
        self.from_homozygous_minor = dart_control.raw.homozygous_minor
        self.from_heterozygous = dart_control.raw.heterozygous
        self.from_missing = dart_control.raw.missing

        self.to_homozygous_major = ("A", "A")
        self.to_homozygous_minor = ("B", "B")
        self.to_heterozygous = ("A", "B")
        self.to_missing = ("-9", "-9")

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

        self.fasta_path = os.path.join(path, file_name)

        with open(self.fasta_path, "w") as fasta_file:
            SeqIO.write(seqs, fasta_file, "fasta")

    def write_snps(self, path=os.getcwd(), mode="dart", filtered=True, sep=','):

        """
        Write the filtered SNPs in various formats:

        1. DArT
        2. PLINK
        3. Structure

        """

        out_file = os.path.join(path, self.qc.project + "_dartQC_" + time.strftime("[%H:%M:%S]"))

        if filtered:
            n = len(self.qc.snps_filtered)
            out_data = self.qc.snps_filtered
        else:
            n = len(self.qc.snps_total)
            out_data = self.qc.snps_total

        order = sorted(out_data.keys())

        ### DArT Format ###

        if mode == "dart":

            out_file += '.csv'

            head = [["AlleleID", "CloneID", "AlleleSequence", "SNP", "SnpPosition", "CallRate", "OneRatioRef", "OneRatioSnp",
                     "FreqHomRef", "FreqHomSnp", "FreqHets",	"PICRef", "PICSnp", "AvgPIC", "AvgCountRef", "AvgCountSnp",
                     "RepAvg"] + self.raw.sample_names]

            out_head = self.raw.header + head

            with open(out_file, 'w') as outfile:
                dart_writer = csv.writer(outfile, delimiter=sep)

                dart_writer.writerows(out_head)

                for allele_id in order:
                    entry = out_data[allele_id]

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

                    dart_writer.writerows([allele_one, allele_two])

        ### PLINK and STRUCTURE Format ###

        elif mode == 'plink' or 'structure':

            # Inititate Tricoder and get the encoding scheme for translating calls from DArT

            tricoder = Tricoder()
            tricoder.get_encoding_scheme(self)

            # Zip SNPs by Allele ID

            snp_rows = [list(zip(out_data[allele_id]["calls"][0], out_data[allele_id]["calls"][1]))
                        for allele_id in order]

            # Decode SNP Rows
            snp_rows_decoded = tricoder.decode(snp_rows)

            # Turn into Numpy Array
            snp_rows_numpy = numpy.asarray(snp_rows_decoded)

            # Transpose Array
            snps_by_sample = snp_rows_numpy.transpose(1, 0, 2)

            genotypes = [sample.flatten().tolist() for sample in snps_by_sample]

            names = self.raw.sample_names

            if self.raw.meta:
                try:
                    pops = [self.raw.meta[name] for name in names]
                except KeyError:
                    raise(KeyError("Sample names from DArT file do not correspond to sample names from meta data file."))
            else:
                # Placeholder
                pops = ["PopEye" for name in names]

            # Remove all whitespace from names, since it messes with PLINK
            names = ["_".join(name.split()) for name in names]
            pops = ["_".join(pop.split()) for pop in pops]

            if mode == 'plink':

                ped_file = out_file + '.ped'
                map_file = out_file + '.map'

                paternal = ["0"] * len(names)
                maternal = ["0"] * len(names)
                sex = ["0"] * len(names)
                phenotype = ["-9"] * len(names)

                # Zip all necessary Data for PLINK

                plink = zip(pops, names, paternal, maternal, sex, phenotype, genotypes)

                # PED Formatting

                ped_data = []
                for row in plink:
                    new_row = list(row[:6])
                    for geno in row[6]:
                        new_row.append(geno)
                    ped_data.append(new_row)

                with open(ped_file, 'w') as ped_out:
                    ped_writer = csv.writer(ped_out, delimiter=' ')
                    ped_writer.writerows(ped_data)

                # MAP Formatting

                map_data = [["0", snp_id, "0", "0"] for snp_id in out_data.keys()]

                with open(map_file, 'w') as map_out:
                    ped_writer = csv.writer(map_out, delimiter=' ')
                    ped_writer.writerows(map_data)

            elif mode == "structure":

                structure_file = out_file + '.str'

                # One row format STRUCTURE, no Header

                structure = zip(names, pops, genotypes)

                structure_data = []
                for row in structure:
                    new_row = list(row[:2])
                    for geno in row[2]:
                        new_row.append(geno)
                    structure_data.append(new_row)

                with open(structure_file, 'w') as str_file:
                    str_writer = csv.writer(str_file, delimiter=' ')
                    str_writer.writerows(structure_data)

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

        self._tmp_path = os.path.join(os.getcwd(), "dart_qc_tmp")

        ### Logger ###

        self.starting_parameters = textwrap.dedent("""
        DART QC
        ---------------------------------------------
        DArT Quality Control Pipeline
        James Cook University
        ARC Hub for Advanced Prawn Breeding
        Eike J. Steinig, Jarrod Guppy, Kyall Zenger
        ---------------------------------------------


        Initializing DArT Quality Control
        ----------------------------------

        Project:  {0}
        SNPs:     {1}
        Samples:  {2}
        """ .format(self.project, self.raw.snp_number, self.raw.sample_number))

        self.log = [self.starting_parameters]

        self.verbose = True

        if self.verbose:
            print(self.starting_parameters)

    def find_duplicate_clones(self):

        """Search for duplicate clone IDs in SNPs and remove from a copy of the raw data. """

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

        clone_msg = textwrap.dedent("""
        SEARCHING FOR DUPLICATE CLONE IDs
        ------------------------------------------

        Number of duplicated clone IDs:    {0}
        Number of duplicated SNPs:         {1}
        Number of unique SNPs:             {2}
        """ .format(len(self.clones_duplicate), len(self.snps_duplicate),len(self.snps_unique)))

        if self.verbose:
            print(clone_msg)

    def find_identity_clusters(self, identity=0.95, word_size=10, description_length=0, cdhit_path=None):

        """
        Clusters the reference allele sequences with CDHIT-EST and parses the clusters for selecting and
        retaining best sequences.

        CD-HIT returns slightly different cluster configurations for each run due to greedy incremental algorithm,
        but little variation observed between runs in the data for P. monodon. Know thyself!

        """

        os.makedirs(self._tmp_path, exist_ok=True)

        dart_writer = DartWriter(self)
        dart_writer.write_fasta()

        if cdhit_path is None:
            cdhit_path = "cdhit-est"

        file_name = self.project + "_IdentityClusters_" + str(identity)

        out_file = os.path.join(self._tmp_path, file_name)
        cluster_file = os.path.join(self._tmp_path, file_name + '.clstr')

        if self.verbose:
            print("...")
        with open(os.devnull, "w") as devnull:
            call([cdhit_path, "-i", dart_writer.fasta_path, "-o", out_file, "-c", str(identity), "-n", str(word_size),
                  "-d", str(description_length)], stdout=devnull)

        clstr_reader = DartReader()
        clstr_reader.parse_cdhit(cluster_file)

        self.cluster_duplicate = clstr_reader.identity_clusters

        self.cluster_snps = clstr_reader.identity_snps

        cluster_msg = textwrap.dedent("""
        CLUSTERING REFERENCE ALLELE SEQUENCES
        -----------------------------------------------

        CDHIT-EST
        Threshold: {0}%

        Warning: clusters are non-deterministic and
        may vary slightly between runs. Please see
        the manual for CD-HIT.

        Number of identity clusters:            {1}
        Number of SNPs in identity clusters:    {2}
        """ .format(identity*100, len(self.cluster_duplicate), self.cluster_snps))

        self.log.append(cluster_msg)

        if self.verbose:
            print(cluster_msg)

    def select_best_identity_seqs(self, selector="maf"):

        """
        Selects the best sequences (SNPs) from the identity clusters (containing more than one sequence) according
        to 'selector', removes the best sequence from the cluster and removes bad sequences from the total
        SNPs, while retaining the best-pick SNP.

        """

        before = len(self.snps_total)

        error_msg = "Error. You need to find identity clusters first, then remove sequences."

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

        retain_cluster_msg = textwrap.dedent("""
        RETAINING BEST CLUSTERED SEQUENCES
        --------------------------------------------------------------

        Retention criterion:                                    {0}
        Number of SNPs found in identity clusters:              {1}
        Total number of SNPs before removing clustered SNPs:    {2}
        Number of clustered SNPs removed:                       {3}
        Total number of SNPs after removing clustered SNPs:     {4}
        """ .format(selector.upper(), before, self.cluster_snps, self.cluster_snps - len(self.cluster_duplicate),
                    len(self.snps_total)))

        self.log.append(retain_cluster_msg)

        if self.verbose:
            print(retain_cluster_msg)

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

        retain_clones_msg = textwrap.dedent("""
        RETAINING BEST DUPLICATE CLONES
        ----------------------------------

        Retention criterion:     {0}
        Retained SNPs:           {1}
        Total unique SNPs:       {2}
        """ .format(selector.upper(), len(self.clones_duplicate), len(self.snps_total)))

        self.log.append(retain_clones_msg)

        if self.verbose:
            print(retain_clones_msg)

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

        filter_msg = textwrap.dedent("""
        FILTERING SNPs #{0}
        ----------------------

        {1} {2} {3} on {4}

        Initial:    {6}
        Removed:    {5}
        Retained:   {7}

        ----------------------
        """ .format(self.filter_count, selector.upper(), comparison, threshold, data.upper(),
                    before - after, before, after))

        self.log.append(filter_msg)

        if self.verbose:
            print(filter_msg)
