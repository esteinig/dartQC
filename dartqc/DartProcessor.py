import numpy

from dartqc.DartReader import DartReader
from dartqc.DartUtils import stamp


class Preprocessor(DartReader):

    def __init__(self, call_data, call_attributes):

        """
        Give the Preprocessor the call data and inherit from DartReader. Use the attributes of DartReader
        (self.data, self.attributes) to process the read count matrix and collapse replicate read counts for samples
        present in the call data.
        """

        DartReader.__init__(self)

        self.name = "preprocessor"

        self.call_data = call_data
        self.call_attributes = call_attributes

        self.call_names = call_attributes["sample_names"]

        self.replicates = {}

        self._set_log()

    def _set_log(self):

        self.call_attributes["modules"][self.name] = {
            "results": {},
            "settings": {},
            "states": {}  # States are other parameters of interest not necessary results or settings.
        }

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
        read the data without encoding, i.e. in list of tuple format for read counts per call.
        """

        self.read_double_row(file=file, encode=False, numeric=True)

    def check_concordance(self):

        if set(self.sample_names) != set(self.call_names):
            stamp("Sample names from the read count file are not the same as sample names from the data file.")
            stamp("Sample difference, present in one but not the other data:")
            for sample in set(self.sample_names).difference(set(self.call_names)):
                stamp(sample)
            exit(0)
        else:
            stamp("Concordance between sample names in call and count data, all is good.")

        if len(self.call_data) != len(self.data):

            diff = set(self.call_data.keys()).difference(set(self.data.keys()))
            inter = set(self.call_data.keys()).intersection(set(self.data.keys()))

            stamp("Number of SNPs are different, there are:", len(self.call_data), "SNPs in the called set and",
                  len(self.data), "SNPs in the raw set.")

            stamp(len(diff), "SNPs have a different ID. Keeping the intersection of", len(inter), "SNPs...")

            self.data = {k: v for (k, v) in self.data.items() if k in inter}
            self.call_data = {k: v for (k, v) in self.call_data.items() if k in inter}

        if set(self.call_data.keys()) != set(self.data.keys()):
            stamp("SNP IDs are not the same, removal not effective, please re-format your data.")

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
        stamp("Number of missing in call data:", call_missing)

        snp_order = sorted(self.data.keys())
        reduced_counts = {}

        stamp("Finding replicate columns...")

        self.get_replicates()

        stamp("Ordering count data by SNPs...")

        counts = [self.data[snp]["calls"] for snp in snp_order]

        count_array = numpy.asarray(counts)

        stamp("Sum-collapsing replicates...")

        columns = [numpy.sum(count_array[:, self.replicates[sample]], axis=1).tolist() for sample in self.call_names]

        reduced_array = list(zip(*columns))

        for i, snp in enumerate(snp_order):
            reduced_counts[snp] = reduced_array[i]

        replaced = 0
        total = 0

        stamp("Replacing low counts with missing...")

        for snp, counts in reduced_counts.items():
            filter_vector = [False if sum(allele_counts) <= threshold else True for allele_counts in counts]
            self.call_data[snp]["calls"] = [call if filter_vector[i] else "-" for i, call in
                                            enumerate(self.call_data[snp]["calls"])]

            replaced += filter_vector.count(False)
            total += len(filter_vector)

        replaced -= call_missing

        stamp("Pre-processing silenced {r}/{t} calls {p}".format(r=replaced, t=total,
                                                                 p=format((replaced/total)*100, ".2f")))

        self.call_attributes["modules"][self.name]["results"] = {
            "total_calls": total,
            "replaced_calls": replaced,
            "before_missing": call_missing,
            "after_missing": call_missing + replaced
        }

        self.call_attributes["modules"][self.name]["settings"] = {
            "read_count_sum_threshold": threshold
        }

    def get_data(self):

        return self.call_data, self.call_attributes
