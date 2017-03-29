import numpy
from DartReader import DartReader


class Preprocessor(DartReader):

    def __init__(self, call_data, call_attributes):

        """
        Give the Preprocessor the call data and inherit from DartReader. Use the attributes of DartReader
        (self.data, self.attributes) to process the read count matrix and collapse replicate read counts for samples
        present in the call data.
        """

        DartReader.__init__(self)

        self.call_data = call_data
        self.call_attributes = call_attributes

        self.call_names = call_attributes["sample_names"]

        self.replicates = {}

    def get_replicates(self):

        """
        Get replicate indices across individuals for summation of read calls.

        """

        for i, sample in enumerate(self.sample_names):
            if sample not in self.replicates.keys():
                self.replicates[sample] = [i]
            else:
                self.replicates[sample] += [i]

    def read_count_data(self, file):

        """
        Alternative call to reading double row format for inputting and transforming the read count matrix,
        read the dat without encoding, i.e. in list of tuple format for read counts per call.
        """

        self.read_double_row(file=file, encode=False, numeric=True)

    def check_concordance(self):

        if set(self.sample_names) != set(self.call_names):
            print("Sample names from the read count file are not the same as sample names from the data file.")
            print("Difference:", set(self.sample_names).difference(set(self.call_names)))
        else:
            print("Concordance between sample names in call and count data, all is good.")

        if len(self.call_data) != len(self.data):

            diff = set(self.call_data.keys()).difference(set(self.data.keys()))
            inter = set(self.call_data.keys()).intersection(set(self.data.keys()))

            print("Number of SNPs are different, there are:", len(self.call_data), "SNPs in the call data and",
                  len(self.data), "SNPs in the count data.", len(diff),
                  "SNPs have a different ID. Keeping the intersection of", len(inter), "SNPs...")

            self.data = {k: v for (k, v) in self.data.items() if k in inter}
            self.call_data = {k: v for (k, v) in self.call_data.items() if k in inter}

        if set(self.call_data.keys()) != set(self.data.keys()):
            print("SNP IDs are not the same, removal not effective, please re-format your data.")

    def get_missing(self):

        """

        Get count of missing alleles in call data, to subtract from total / replaced in filtering read counts

        """

        missing = 0

        for k, v in self.call_data.items():
            missing += v["calls"].count("-")

        return missing

    def filter_read_counts(self, threshold=7):

        """
        1. Transform read count matrix to numpy array, ordered by allele IDs.
        2. Sum-collapse replicate columns in the order of sample names from the original data (sample_names)
        3. Construct the reduced array and assign each call in a dictionary the allele ID
        4. For each allele in the dictionary, construct a boolean vector with True if the sum of the two allele counts
           is smaller than the threshold value, otherwise False
        5. Use this vector in the same iteration to assign missing to all calls in the original data

        """

        self.check_concordance()
        call_missing = self.get_missing()
        print("Number of missing in call data:", call_missing)

        snp_order = sorted(self.data.keys())
        reduced_counts = {}

        print("Finding replicate columns...")

        self.get_replicates()

        print("Ordering count data by SNPs...")

        counts = [self.data[snp]["calls"] for snp in snp_order]

        count_array = numpy.asarray(counts)

        print("Sum-collapsing replicates...")

        columns = [numpy.sum(count_array[:, self.replicates[sample]], axis=1).tolist() for sample in self.call_names]

        reduced_array = list(zip(*columns))

        for i, snp in enumerate(snp_order):
            reduced_counts[snp] = reduced_array[i]

        replaced = 0
        total = 0

        print("Replacing low counts with missing...")

        for snp, counts in reduced_counts.items():
            filter_vector = [False if sum(allele_counts) <= threshold else True for allele_counts in counts]
            self.call_data[snp]["calls"] = [call if filter_vector[i] else "-" for i, call in
                                            enumerate(self.call_data[snp]["calls"])]

            replaced += filter_vector.count(False)
            total += len(filter_vector)

        replaced -= call_missing

        print("Pre-processing found an additional", replaced, "out of a total of", total, "calls = ", format((replaced/total)*100, ".2f"),
                                                              "% (excluding calls already designated as missing) "
                                                              " with an "
                                                              "insufficient total read depth across "
                                                              "two calls ( <=", threshold, ").\n")

    def get_data(self):

        return self.call_data, self.call_attributes
