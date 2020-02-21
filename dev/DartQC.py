import time
import numpy

from matplotlib import pyplot

import Configuration
from DartModules import MarkerModule, RedundancyModule, SummaryModule
from DartProcessor import Preprocessor
from DartReader import DartReader
from DartWriter import DartWriter
from Graphs import Graphs


class DartQC:
    def __init__(self, project, working_dir):
        self.project = project
        self.working_dir = working_dir

        self.data = None
        self.attrs = None
        self.filtered_data = None
        self.filtered_attrs = None

    def validate_files(self, data_file, read_counts_file):
        pass

    # Read the genotype data from files into JSON structures then threshold/filter genotype calls by read counts.
    def pre_process(self, data_file, data_mapping, read_counts_file, read_counts_mapping, read_counts_threshold=5):
        print("Pre-Processing...")
        # Reading initial double-row read calls with settings adjusted t input file format:
        dart_reader = DartReader()
        dart_reader.set_options(project=self.project, clone_col=data_mapping["clone_col"], id_col=data_mapping["id_col"],
                                sample_row=data_mapping["sample_row"], data_start_row=data_mapping["data_start_row"],
                                out_path=self.working_dir)
        dart_reader.read_double_row(file=self.working_dir + data_file, split_clone=True)
        data, attributes = dart_reader.get_data()

        # Read the read call data into the Preprocessor:
        pp = Preprocessor(call_data=data, call_attributes=attributes)

        # Setting options for raw read count file:
        pp.set_options(project=self.project, clone_col=read_counts_mapping["clone_col"], id_col=read_counts_mapping["id_col"],
                       sample_row=read_counts_mapping["sample_row"], data_start_row=read_counts_mapping["data_start_row"],
                       call_start_col=read_counts_mapping["call_start_col"],
                       sample_start_col=read_counts_mapping["sample_start_col"])

        # Reading the raw read counts:
        pp.read_count_data(self.working_dir + read_counts_file)

        self.read_data = pp.get_read_data()

        # Create bar graphs and scatter plots before thresholding
        self.create_static_plots(data, self.read_data)
        self.create_plots(data, self.read_data, attributes, "Original")

        # Set all calls to missing < threshold (sum of minor and major read counts for SNP)
        pp.filter_read_counts(threshold=read_counts_threshold)
        self.counts = pp.get_counts()

        # Export data and attributes for further use in the filtering modules...
        data, attributes = pp.get_data()

        # Create new set of plots after thresholding based on read counts.
        self.create_plots(data, self.read_data, attributes, "Thresholded")

        self.data = data
        self.attrs = attributes

        # Writing these data to JSON, as pre-processing may take a while...
        dart_writer = DartWriter(data, attributes)
        dart_writer.write_json(self.project)

    # Do MAF, Call rate, Rep and read ref filtering
    def filter(self, maf=.02, call_rate=.5, read_ref=0, rep=.9):
        print("Filtering...")

        # Run provided MAF, Call rate, rep and read ref filters
        # Reading the data from JSON after Preprocessing:
        data = self.data
        attributes = self.attrs

        if data is None:
            dart_reader = DartReader()
            data, attributes = dart_reader.read_json(data_file=self.working_dir + "//" + self.project + "_data.json",
                                                     attribute_file=self.working_dir + "//" + self.project + "_attr.json")

        # Initializing MarkerModule, which handles filtering:

        mm = MarkerModule(data=data, attributes=attributes)

        # Indexing all filter values defined above
        # True/False for retaining the SNP across all filter values and SNPs

        mm.filter_data([maf], parameter="maf", comparison="<=")
        mm.filter_data([call_rate], parameter="call_rate", comparison="<=")
        mm.filter_data([read_ref], parameter="read_count_ref", comparison="<=")
        mm.filter_data([rep], parameter="rep_average", comparison="<=")

        # Deploying filter across pre-indexed values (values given here must be present in lists above)
        # Filtered data is exported from module for further use:
        data, attributes = mm.get_data(multiple=[("maf", maf),
                                                 ("call_rate", call_rate),
                                                 ("rep_average", rep),
                                                 ("read_count_ref", read_ref)
                                                 ])

        self.filtered_data = data
        self.filtered_attrs = attributes

        # Writing these data to JSON
        dart_writer = DartWriter(data, attributes)
        dart_writer.write_json(self.project + "_filtered")

    # Remove duplicates and clusters
    def clusters_and_dups(self):
        print("Removing duplications and clusters...")
        # Reading the data from JSON after Preprocessing:

        data = self.filtered_data
        attributes = self.filtered_attrs

        if data is None:
            dart_reader = DartReader()
            data, attributes = dart_reader.read_json(
                data_file=self.working_dir + "//" + self.project + "_filtered_data.json",
                attribute_file=self.working_dir + "//" + self.project + "_filtered_attr.json")

        # Need CD-HIT (cdhit-est) in PATH
        # Ubuntu: sudo apt install cd-hit
        # Initializing redundancy module
        rm = RedundancyModule(data=data, attributes=attributes)

        # Indexing duplicate and identity clusters:
        rm.remove_duplicates(selector_list=["maf", "read_count_ref"])
        rm.remove_clusters(selector_list=["maf", "read_count_ref"], cdhit_path=Configuration.CDHIT_PATH)

        # Export data with duplicates and clustered SNPs removed:
        data, attributes = rm.get_data(duplicates=True, clusters=True)

        self.create_plots(data, self.read_data, attributes, "Final")

        # Summary module for writing a summary of SNP parameters:
        sm = SummaryModule(data=data, attributes=attributes)
        sm.write_snp_summary(summary_parameters=["maf", "call_rate", "rep_average", "read_count_ref"])

        # TODO: Create output plots

        dart_writer = DartWriter(data, attributes)
        dart_writer.write_json(self.project + "_final")

        dart_writer.write_plink(self.project)

    # Some graphs won't be changed based on filtering (eg. read counts, repAvg & freq Hetz)
    def create_static_plots(self, data, read_data):
        print("Creating static plots")

        Graphs.read_counts_per_individ(read_data, self.working_dir + self.project + "_IndividReadCounts" + ".jpg")
        Graphs.avg_reads_per_snp(read_data, self.working_dir + self.project + "_AvgReadsPerSNP" + ".jpg")
        Graphs.avg_rep_across_snp(data, self.working_dir + self.project + "_AvgRepAcrossSNP" + ".jpg")
        Graphs.het_across_snp(data, self.working_dir + self.project + "_HetAcrossSNP" + ".jpg")

    def create_plots(self, data, read_data, attrs, name):
        print("Creating " + name + " plots")

        mm = MarkerModule(data=data, attributes=attrs)
        snpMAF = {k: v["maf"] for (k, v) in data.items()}

        Graphs.call_rates_across_snp(data, read_data, self.working_dir + self.project + "_" + name + "_CallRatesAcrossSNP" + ".jpg")
        Graphs.call_rate_across_individ(data, read_data, self.working_dir + self.project + "_" + name + "_CallRatesAcrossIndivid" + ".jpg")
        Graphs.maf_across_snp(snpMAF, self.working_dir + self.project + "_" + name + "_MAFAcrossSNP" + ".jpg")
        Graphs.maf_to_read_count(snpMAF, read_data, self.working_dir + self.project + "_" + name + "_MAFToReadCount" + ".jpg")
        Graphs.call_rate_to_maf(data, snpMAF, self.working_dir + self.project + "_" + name + "_CallRateToMAF" + ".jpg")

    @staticmethod
    def create_file_mapping(clone_col, id_col, sample_row, data_start_row, call_start_col=None, sample_start_col=None):
        mapping = {
            "clone_col": clone_col,
            "id_col": id_col,
            "sample_row": sample_row,
            "data_start_row": data_start_row,
            "call_start_col": call_start_col,
            "sample_start_col": sample_start_col,
        }
        return mapping


program_start = time.time()

dartQc = DartQC("prawns", "./testData/")
dartQc.pre_process("prawn_data_double.csv", DartQC.create_file_mapping(1,2,7,8), "prawn_read_counts.csv", DartQC.create_file_mapping(2,1,4,7,33,33))
dartQc.filter()
dartQc.clusters_and_dups()

# DartQC.create_bar_graph([4, 6, 8, 5], "Test Title", ["a", "b", "c", "d"], "X Label", "Y Label", 'foo.png')
# DartQC.create_scatter_plot([1, 2, 3, 4], [4, 6, 8, 5], "Test Title", ["a", "b", "c", "d"], "X Label", "Y Label",
#                            'bar.png')

print("Program run time: " + str(round((time.time() - program_start), 2)) + "s")