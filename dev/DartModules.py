import os
import csv
import numpy
import operator

from subprocess import call
from scipy import stats

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

from DartMessages import DartMessages

import time

from subprocess import run


class SummaryModule:

    def __init__(self, data=None, attributes=None, out_path=None):

        self.data = data
        self.attributes = attributes

        if out_path is None:
            self.out_path = attributes["out_path"]
        else:
            self.out_path = out_path

        os.makedirs(self.out_path, exist_ok=True)

    def write_snp_summary(self, file="snp_summary.csv", summary_parameters=None, sort=False):

        if summary_parameters is None:
            summary_parameters = ["maf", "read_count_ref", "hwe"]

        out_file = os.path.join(self.out_path, file)

        out_data = [["id"] + summary_parameters]

        snps = [[snp] + [data[parameter] for parameter in summary_parameters] for snp, data in self.data.items()]

        if sort:
            snps = sorted(snps, key=operator.itemgetter(*[i for i in range(1, len(summary_parameters)+1)]),
                                reverse=True)

        out_data += snps

        with open(out_file, "w") as snp_summary:
            writer = csv.writer(snp_summary)
            writer.writerows(out_data)

    def write_matrix(self, combination_matrix, r_matrix=None, file="combination_table.csv", r_file="r_matrix.csv"):

        out_file = os.path.join(self.out_path, file)

        with open(out_file, "w") as table_file:
            writer = csv.writer(table_file)
            writer.writerows(combination_matrix)

        if r_matrix is not None:
            out_r = os.path.join(self.out_path, r_file)
            with open(out_r, "w") as out_r_file:
                writer = csv.writer(out_r_file)
                writer.writerows(r_matrix)

########################################################################################################################


class QualityControl:

    def __init__(self, data, attributes):

        self.data = data                                # Dictionary holds data from DartReader
        self.attributes = attributes

        self.messages = DartMessages()

        self.sample_size = attributes["sample_size"]
        self.sample_names = attributes["sample_names"]

        self.missing = attributes["missing"]
        self.homozygous_major = attributes["homozygous_major"]
        self.homozygous_minor = attributes["homozygous_minor"]
        self.heterozygous = attributes["heterozygous"]

        self.project = attributes["project"]
        self.out_path = attributes["out_path"]

########################################################################################################################


class PopulationModule(QualityControl):

    def __init__(self, data, attributes):

        QualityControl.__init__(self, data, attributes)

        self.samples = attributes["samples"]                  # Dictionary of individual parameters by name, unique.
        self.sample_names = attributes["sample_names"]        # List of ordered unique names, same order as calls for SNPs

        self.populations = {}                                 # Dictionary of populations and list of sample indices
        self.monomorphics = {}
        self.genotypes = numpy.asarray([[]])

        self.filters = {"mind": {}, "mono": {}}                 # Filters for removing: Samples(Mind), SNPs (Mono)

        self._get_sample_indices()

    def get_data(self, threshold=None, parameter="mind", recalculate=True):

        """
        Get filtered data either for missingness per individual or monomorphic SNPs per Population, e.g.
        Samples < 0.7 MIND or SNP monomorphic in 3 populations
        """

        if threshold is None:
            return self.data, self.attributes

        if parameter == "mind":

            filtered_data = {}
            for snp, data in self.data.items():
                calls = [snp_call for i, snp_call in enumerate(data["calls"])
                         if i not in self.filters["mind"][threshold]]
                filtered_data[snp]["calls"] = calls

            # Adjust attribute file:
            attributes = self.attributes.copy()
            attributes["sample_names"] = [name for i, name in enumerate(self.sample_names)
                                          if i not in self.filters["mind"][threshold]]
            attributes["samples"] = {k: v for (k, v) in self.samples.items() if k in attributes["sample_names"]}
            attributes["sample_size"] = self.attributes["sample_size"] - len(self.filters["mind"][threshold])

            # Re-calculate statistics for SNPs, except Monomorphic (need Fix):
            if recalculate:
                marker = MarkerModule(filtered_data, attributes)
                filtered_data, attributes = marker.get_data(threshold=None)

            return filtered_data, attributes

        elif parameter == "mono":

            return {k: v for (k, v) in self.data.items() if v["mono"] != threshold}, self.attributes

    def filter_data(self, thresholds, parameter="mind", comparison=">"):

        """
        Create the filter dictionary with missingness per individual > threshold, or monomorphic in number of
        populations == threshold...

        """

        self._calculate_mind()
        self._calculate_monomorphics()

        for threshold in thresholds:

            if parameter == "mind":
                if comparison == ">":
                    filtered = [self.sample_names.index(name) for name, data in self.samples.items()
                                if data["mind"] > threshold]
                elif comparison == "<=":
                    filtered = [self.sample_names.index(name) for name, data in self.samples.items()
                                if data["mind"] <= threshold]
                elif comparison == "==":
                    filtered = [self.sample_names.index(name) for name, data in self.samples.items()
                                if data["mind"] == threshold]
                try:
                    self.filters["mind"][threshold] = filtered
                except KeyError:
                    self.filters["mind"] = {threshold: filtered}

            elif parameter == "mono":
                filtered = [snp for snp, data in self.data.items() if data["mono"] == threshold]

                try:
                    self.filters["mono"][threshold] = filtered
                except KeyError:
                    self.filters["mono"] = {threshold: filtered}

    def _get_sample_indices(self):

        """
        Make a dictionary from the individual parameters and ordered name list to get indices of samples for each
        population.

        """

        for name, parameters in self.samples.items():
            name_index = self.sample_names.index(name)
            if parameters["pop"] not in self.populations:
                self.populations[parameters["pop"]] = [name_index]
            else:
                self.populations[parameters["pop"]].append(name_index)

    def _transform(self):

        """ Sort the markers alphabetically and transform the data to a genotype array (row per individual) """

        snp_order = sorted(self.data.keys())                # SNPs are ordered alphabetically.

        snp_rows = [self.data[allele_id]["calls"] for allele_id in snp_order]

        return numpy.asarray(snp_rows).transpose()       # Same order as input names and populations.

    def _calculate_mind(self):

        """
        Calculate rate of missing markers per individual and add to individual dictionary.
        """

        genotypes = self._transform()

        for i in range(len(genotypes)):
            missing = genotypes[i].tolist().count(self.missing)/len(genotypes[i])
            self.samples[self.sample_names[i]]["mind"] = missing

    def _calculate_monomorphics(self):

        for pop in self.populations.keys():
            self.monomorphics[pop] = []

        for snp, data in self.data.items():
            number = 0
            for pop, indices in self.populations.items():

                calls = set([snp_call for i, snp_call in enumerate(data["calls"])
                             if i in indices and snp_call != self.missing])
                if len(calls) == 1:
                    self.monomorphics[pop].append(snp)
                    number += 1
            self.data[snp]["mono"] = number

    def _calculate_maf(self, calls):

        """
        Calculates minor allele frequency for a single SNP. Pass a one-row format list (zip on double-row)
        of allele calls. Returns the minimum allele frequency for processing.
        """

        adjusted_samples = self.sample_size - calls.count(self.missing)
        het_count = calls.count(self.heterozygous)
        try:
            freq_allele_one = (calls.count(self.homozygous_major) + (het_count/2)) / adjusted_samples
            freq_allele_two = (calls.count(self.homozygous_minor) + (het_count/2)) / adjusted_samples
        except ZeroDivisionError:
            return 0

        return min(freq_allele_one, freq_allele_two)

########################################################################################################################


class RedundancyModule(QualityControl):

    def __init__(self, data, attributes):

        QualityControl.__init__(self, data, attributes)

        self.duplicates = {}
        self.retained_duplicates = []
        self.removed_duplicates = []

        self.clusters = {}
        self.retained_sequences = []
        self.removed_sequences = []

        self.identity = 0.95            # Identity used in Cluster Removal

    def get_data(self, duplicates=True, clusters=True, redundant=False):

        data = self.data

        if redundant:
            duplicate_list = self.retained_duplicates
            cluster_list = self.retained_sequences
        else:
            duplicate_list = self.removed_duplicates
            cluster_list = self.removed_sequences

        # Check if logic is correct:

        if duplicates and clusters:
            filters = [("duplicates", duplicate_list), ("clusters", cluster_list)]
        elif clusters and not duplicates:
            filters = [("clusters", cluster_list)]
        elif duplicates and not clusters:
            filters = [duplicate_list]
        else:
            filters = list()

        if filters:
            for type, filter_list in filters:
                before = len(data)
                data = {k: v for (k, v) in data.items() if k not in filter_list}
                after = len(data)
                print(self.messages.get_redundancy_message(type, before, before-after, after))

        return data, self.attributes

    def remove_duplicates(self, selector="maf", target="clone_id", selector_list=None):

        """Search for duplicate clone IDs in SNPs and remove from data. Wrapper. """

        self._find_duplicates(target=target)
        self._select_duplicates(selector=selector, selector_list=selector_list)

    def remove_clusters(self, identity=0.95, target="allele_seq_ref", selector="maf",
                        selector_list=None, cdhit_path=None):

        """ Search for sequence identity clusters and remove them from data. Wrapper. """

        self.identity = identity

        self._find_clusters(target=target, identity=identity, cdhit_path=cdhit_path)
        self._select_clusters(selector=selector, selector_list=selector_list)

    # Private functions for Redundancy Module #

    def _find_clusters(self, target="allele_seq_ref", identity=0.95, word_size=10, description_length=0,
                       cdhit_path=None):

        """
        Clusters the reference allele sequences with CDHIT-EST and parses the clusters for selecting and
        retaining best sequences.

        CD-HIT returns slightly different cluster configurations for each run due to greedy incremental algorithm,
        but little variation observed between runs in the data for P. monodon. Know thyself!

        """

        os.makedirs(self.out_path, exist_ok=True)

        fasta_path = self._write_fasta(target)

        cluster_path = self._run_cdhit(fasta_path=fasta_path, identity=identity, word_size=word_size,
                                       description_length=description_length, cdhit_path=cdhit_path)

        self.clusters = self._parse_cdhit(cluster_path)

    def _select_clusters(self, selector="maf", selector_list=None):

        """ Select best markers from clusters by selector. """

        for cluster, cluster_members in self.clusters.items():
            best_sequence = self._compare_entries(cluster_members, selector=selector, selector_list=selector_list)
            self.retained_sequences.append(best_sequence)
            self.removed_sequences += [member for member in cluster_members if member != best_sequence]

    def _write_fasta(self, target="allele_seq_ref"):

        """ Write fasta file of sequences with SNP IDs for CD-HIT. """

        file_name = os.path.join(self.out_path, self.project + "_Seqs")

        seqs = [SeqRecord(Seq(data[target], IUPAC.unambiguous_dna), id=snp_id, name="", description="")
                for snp_id, data in self.data.items()]
        file_name += ".fasta"

        with open(file_name, "w") as fasta_file:
            SeqIO.write(seqs, fasta_file, "fasta")

        return file_name

    def _run_cdhit(self, fasta_path, identity=0.95, word_size=10, description_length=0, cdhit_path=None):

        """ Run CDHIT-EST for sequences, install with sudo apt install cd-hit on Ubuntu """

        print(self.messages.get_cdhit_message(identity))

        if cdhit_path is None:
            cdhit_path = "cdhit-est"

        file_name = self.project + "_IdentityClusters_" + str(identity)

        out_file = os.path.join(self.out_path, file_name)
        cluster_path = os.path.join(self.out_path, file_name + '.clstr')

        with open(os.devnull, "w") as devnull:
            call([cdhit_path, "-i", fasta_path, "-o", out_file, "-c", str(identity), "-n", str(word_size),
                  "-d", str(description_length)], stdout=devnull)
            # print("CD Hit result: " + str(retVal))
            # run(["C:\\Users\\jc229483\\Desktop\\Repositories\\DartQC\\cdhit.bat"], shell=True)
            # cmdCall = cdhit_path + " -i " + fasta_path + " -o " + out_file + " -c " + str(identity) + " -n " + str(word_size) + " -d " + str(description_length)
            # print("CD Hit call: " + cmdCall)
            # time.sleep(10)
            # os.system(cmdCall)

        return cluster_path

    def _parse_cdhit(self, file):

        """
        Parses the CDHIT cluster output and retains only cluster with more than one sequence for selection and
        removal from total SNPs in dictionary: {cluster_id: snp_ids}

        """

        identity_clusters = {}

        with open(file, "r") as clstr_file:
              nrows = len(list(clstr_file))

        with open(file, "r") as clstr_file:
            clstr_reader = csv.reader(clstr_file, delimiter=' ')

            cluster_id = 0
            ids = []

            row_index = 1
            for row in clstr_reader:
                if row[0] == ">Cluster":
                    if len(ids) > 1:
                        identity_clusters[cluster_id] = ids
                    ids = []
                    cluster_id += 1
                    row_index += 1
                else:
                    allele_id = self._find_between(row[1], ">", "...")
                    ids.append(allele_id)
                    # EOF
                    if row_index == nrows:
                        if len(ids) > 1:
                            identity_clusters[cluster_id] = ids
                    row_index += 1

        return identity_clusters

    @staticmethod
    def _find_between(s, first, last):

        """Finds the substring between two strings."""

        try:
            start = s.index(first) + len(first)
            end = s.index(last, start)
            return s[start:end]
        except ValueError:
            return ""

    def _find_duplicates(self, target="clone_id"):

        """ Count duplicates in target category (usually clone ID) """

        clone_counts = {}

        for k, v in self.data.items():

            clone_id = v[target]
            if clone_id not in clone_counts.keys():
                clone_counts[clone_id] = {"count": 1, "allele_ids": [k]}
            else:
                clone_counts[clone_id]["count"] += 1
                clone_counts[clone_id]["allele_ids"].append(k)

        self.duplicates = {k: v for (k, v) in clone_counts.items() if v["count"] > 1}

    def _select_duplicates(self, selector="maf", selector_list=None):

        """ Select best SNP from duplicate clusters. """

        for clone, clone_data in self.duplicates.items():
            best_clone = self._compare_entries(clone_data["allele_ids"], selector=selector, selector_list=selector_list)
            self.retained_duplicates.append(best_clone)
            self.removed_duplicates += [marker for marker in clone_data["allele_ids"] if marker != best_clone]

    def _compare_entries(self, ids, selector="maf", selector_list=None):

        """
        Gets data from dictionary for each duplicate SNP according to 'selector'
        and returns the allele identification of the best entry.

        Selector list currently sorts descending, that is all selector values must be ranked highest value ("best") -
        this is the case for MAF, Call Rate, Rep, Read Counts ...

        Later rank the data by QC Score.

        """

        if selector_list is None:
            entries_stats = [[i, self.data[i][selector]] for i in ids]
            entries_ranked = sorted(entries_stats, key=operator.itemgetter(1), reverse=True)
        else:
            entries_stats = [[i] + [self.data[i][selector] for selector in selector_list] for i in ids]
            entries_ranked = sorted(entries_stats, key=operator.itemgetter(*[i for i in range(1, len(selector_list)+1)]),
                                    reverse=True)

        return entries_ranked[0][0]

########################################################################################################################


class CombinationModule:

    """ Combinator module for MarkerModule. Assess and visualize combinations of parameters for QC. """

    def __init__(self, marker_module):

        if not marker_module.filters:
            ValueError("Data must have been assessed and filtered with MarkerModule.")

        self.filters = marker_module.filters
        self.data = marker_module.data

    def get_matrix(self, parameter_one, parameter_two, values_one, values_two,
                   parameter_three=None, pairwise=False):

        """"
        Pairwise analysis of two values given two parameters and value vectors, outputs data for visualization of
        matrix in Excel and R.

        """

        result_matrix = [[parameter_one + "/" + parameter_two] + values_one]
        r_matrix = []

        try:
            for value_y in values_two:
                result_row = [value_y]
                filtered_y = self.filters[parameter_two][value_y]
                for value_x in values_one:
                    filtered_x = self.filters[parameter_one][value_x]
                    filter_combined = set(filtered_x + filtered_y)
                    r_matrix += [[str(value_y), str(value_x), len(self.data) - len(filter_combined), 100, 100]]
                    result_row.append(len(self.data) - len(filter_combined))        # Number of retained SNPs
                result_matrix.append(result_row)

            return result_matrix, r_matrix

        except KeyError:
            "Can't find parameter or specified value in marker filters, please use MarkerModule"


class MarkerModule(QualityControl):

    """ Analysis module for markers, calculate parameters and filter SNPs. """

    def __init__(self, data, attributes):

        QualityControl.__init__(self, data, attributes)

        self.filters = {}                   # {"maf" : {0.5 : [ "SNP1" , "SNP2" ...] ...} ...}

        self._calculate_parameters()

    def filter_data(self, thresholds, parameter="maf", comparison="<="):

        """
        Filter data by a list of thresholds and parameter / comparsion.
        Filtered marker IDs are stored in filter attribute and can be returned with .get_data().
        """

        if comparison not in ["<=", ">=", "=="]:
            raise(ValueError("Comparison must be one of: <=, >=, =="))

        for threshold in thresholds:

            if comparison == "<=":
                filtered = [k for (k, v) in self.data.items() if v[parameter] <= threshold]
            elif comparison == ">=":
                filtered = [k for (k, v) in self.data.items() if v[parameter] >= threshold]
            else:
                filtered = [k for (k, v) in self.data.items() if v[parameter] == threshold]

            try:
                self.filters[parameter][threshold] = filtered
            except KeyError:
                self.filters[parameter] = {threshold: filtered}

    def get_data(self, threshold=None, parameter="maf", multiple=None):

        """
        Multiple is a dictionary of parameter keys and parameter values (e.g. {"maf": 0.5, "hwe": 0.0001}) returning
        data filtered by multiple parameters.

        """

        if multiple is not None:

            data = self.data
            for parameter, threshold in multiple:
                before = len(data)
                data = {k: v for (k, v) in data.items() if k not in self.filters[parameter][threshold]}
                after = len(data)

                print(self.messages.get_filter_message(parameter, threshold, before, before-after, after))

            return data, self.attributes

        if threshold is None:
            return self.data, self.attributes
        else:

            data = {k: v for (k, v) in self.data.items() if k not in self.filters[parameter][threshold]}
            print(self.messages.get_filter_message(parameter, threshold, len(self.data), len(self.data)-len(data),
                                                   len(data)))

            return data, self.attributes

    # Private functions for parameter calculations #

    def _calculate_maf(self, snp, genotypes):

        """
        Calculates minor allele frequency for a single SNP. Pass a one-row format list (zip on double-row)
        of allele calls. Returns the minimum allele frequency for processing.
        """

        if self.sample_size != len(genotypes):
            print("Warning: Number of samples does not correspond number of allele calls.")

        adjusted_samples = self.sample_size - genotypes.count(self.missing)

        het_count = genotypes.count(self.heterozygous)

        try:
            freq_allele_one = (genotypes.count(self.homozygous_major) + (het_count/2)) / adjusted_samples
            freq_allele_two = (genotypes.count(self.homozygous_minor) + (het_count/2)) / adjusted_samples
        except ZeroDivisionError:
            #print("Detected complete missing data in SNP:", snp)
            return 0

        return min(freq_allele_one, freq_allele_two)

    def _calculate_hwe(self, snp, genotypes):

        """
        Calculates p-value for HWE using ChiSquare statistic: remove missing, get observed counts, get observed
        frequencies, get expected counts, calculate test values using (O-E)**2 / E and return ChiSquare probability
        with 1 degree of Freedom (bi-allelic SNP).

        """

        adjusted_samples = self.sample_size - genotypes.count(self.missing)
        hetero_obs, major_obs, minor_obs = self._get_observed(genotypes)

        try:
            p = (major_obs + (hetero_obs/2)) / adjusted_samples
            q = (minor_obs + (hetero_obs/2)) / adjusted_samples
        except ZeroDivisionError:
            #print("Detected complete missing data in SNP:", snp)
            return 0

        if (p + q) != 1:
            ValueError("Sum of observed allele frequencies (p + q) does not equal one.")

        hetero_exp, major_exp, minor_exp = self._get_expected(p, q, adjusted_samples)

        try:
            hetero_test = ((hetero_obs-hetero_exp)**2)/hetero_exp
            major_test = ((major_obs-major_exp)**2)/major_exp
            minor_test = ((minor_obs-minor_exp)**2)/minor_exp
        except ZeroDivisionError:
            return 0

        return stats.chisqprob(sum([hetero_test, major_test, minor_test]), 1)

    def _get_observed(self, genotypes):

        """" Get observed counts in the genotype for each SNP """

        return genotypes.count(self.heterozygous), genotypes.count(self.homozygous_major),\
               genotypes.count(self.homozygous_minor)

    @staticmethod
    def _get_expected(p, q, adjusted_samples):

        """ Get expected counts under HWE """

        return adjusted_samples*(2*p*q), adjusted_samples*(p**2), adjusted_samples*(q**2)

    def _calculate_call(self, genotypes):

        """ Calculates call rate across samples for a single SNP. """

        return 1 - (genotypes.count(self.missing)/self.sample_size)

    def _calculate_parameters(self):

        """
        Calculate minor allele frequency, call rate and chi-square probability for HWE for each SNP and decorate
        the data dictionary for statistics across SNPs.
        """

        for snp, data in self.data.items():

            self.data[snp]["maf"] = self._calculate_maf(snp, data["calls"])
            self.data[snp]["call_rate"] = self._calculate_call(data["calls"])
            self.data[snp]["hwe"] = self._calculate_hwe(snp, data["calls"])

########################################################################################################################
