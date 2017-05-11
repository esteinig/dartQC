from DartProcessor import *
from DartReader import *
from DartWriter import *
from DartModules import *

import Configuration



def test_preprocessing():

    """

    Prawn processing of the new data, my apologies if this is still a script!

    So, first step: Preprocessing by masking calls whose read count sum is below a certain threshold.

    1. Reading the data in with DartReader - split_clones is set to True, which just means that the actual (non-unique)
       clone ID is extracted from the column Clone ID (e.g. 1234567 from 1234567|F|36--G-C ...)

    2. Running the preprocessor which takes the data and reads in the read count data ('prawn_read_counts.csv') and
       first collapses replicate columns, sums the read counts from the replicates, then sums the read counts for each
       allele (e.g. A:5 + B:4 = 9) and replaces all counts below the 'threshold' value (<=) as 'missing'.

    3. At this stage, the preprocessed data is returned as a temporary (.json) file, which can be read in again for
       filtering - you may notice that the preprocessing takes a while! I think it's in summing the replicates,
       which already uses a numpy array structure and should be pretty fast; will be fixed later.

    """

    # Reading initial double-row read calls with settings adjusted t input file format:
    dart_reader = DartReader()
    dart_reader.set_options(project="PrawnFinal", clone_col=1, id_col=2, sample_row=7, data_start_row=8, out_path="testData")
    dart_reader.read_double_row(file="./testData/prawn_data_double.csv", split_clone=True)
    data, attributes = dart_reader.get_data()

    # Read the read call data into the Preprocessor:
    pp = Preprocessor(call_data=data, call_attributes=attributes)

    # Setting options for raw read count file:
    pp.set_options(project="PrawnPreprocessor", clone_col=2, id_col=1, sample_row=4, data_start_row=7,
                   call_start_col=33, sample_start_col=33)

    # Reading the raw read counts:
    pp.read_count_data("./testData/prawn_read_counts.csv")

    # Set all calls to missing < threshold (sum of minor and major read counts for SNP)
    pp.filter_read_counts(threshold=0)

    # Export data and attributes for further use in the filtering modules...
    data, attributes = pp.get_data()

    # Writing these data to JSON, as pre-processing may take a while...
    dart_writer = DartWriter(data, attributes)
    dart_writer.write_json("prawns_pp_0")


def test_filtering():

    """

    Ok, second step, filtering of preprocessed data:

    1. On the top, define the lists of values that you would like to filter.

    2. Reading the temporary (.json) files so they are available for processing.

    3. Initializing the marker module with the data and attributes: when we call the '.filter_data' the data is indexed
       not actually filtered. So, we construct an index of SNPs that would be filtered at the given values, but only
       when we return the data with '.get_data' we specify which values you want the data to be filtered by. You can
       at this point get multiple 'data' and 'attribute' objects back with different filter values applied if you want
       to write them to file (e.g. data2, data3, data4...)

    4. The filtered data is put into the redundancy module, the 'selector_list' just says when you find a cluster or
       duplicate by clone ID, rank the SNPs by these parameters (i.e. first MAF) and break ties with the next parameter
       (i.e. second Read Count Ref).

    5. The summary module is initialized, which just writes a summary file ("snp_summary.csv") of the given parameters
       to disk.

    """

    # Set up lists of parameters for SNP indexing in the MarkerModule

    maf = [0.02, 0.05, 0.1, 0.15, 0.2]
    call_rate = [0.3, 0.5, 0.6, 0.65, 0.7, 0.8, 0.9]
    rep = [0.80, 0.90, 0.95]
    read_ref = [0, 5, 6, 7, 8, 9, 10, 15, 20]

    # Reading the data from JSON after Preprocessing:
    dart_reader = DartReader()
    data, attributes = dart_reader.read_json(data_file=".//testData//prawns_pp_0_data.json",
                                             attribute_file=".//testData//prawns_pp_0_attr.json")

    # Initializing MarkerModule, which handles filtering:

    mm = MarkerModule(data=data, attributes=attributes)

    # Indexing all filter values defined above
    # True/False for retaining the SNP across all filter values and SNPs

    mm.filter_data(maf, parameter="maf", comparison="<=")
    mm.filter_data(call_rate, parameter="call_rate", comparison="<=")
    mm.filter_data(read_ref, parameter="read_count_ref", comparison="<=")
    mm.filter_data(rep, parameter="rep_average", comparison="<=")

    # Deploying filter across pre-indexed values (values given here must be present in lists above)
    # Filtered data is exported from module for further use:
    data, attributes = mm.get_data(multiple=[("maf", 0.02),
                                             ("call_rate", 0.5),
                                             ("rep_average", 0.90),
                                             ("read_count_ref", 0)
                                             ])

    # Need CD-HIT (cdhit-est) in PATH
    # Ubuntu: sudo apt install cd-hit
    # Initializing redundancy module
    rm = RedundancyModule(data=data, attributes=attributes)

    # Indexing duplicate and identity clusters:
    rm.remove_duplicates(selector_list=["maf", "read_count_ref"])
    rm.remove_clusters(selector_list=["maf", "read_count_ref"], cdhit_path=Configuration.CDHIT_PATH)

    # Export data with duplicates and clustered SNPs removed:
    data, attributes = rm.get_data(duplicates=True, clusters=True)

    # Summary module for writing a summary of SNP parameters:
    sm = SummaryModule(data=data, attributes=attributes)
    sm.write_snp_summary(summary_parameters=["maf", "call_rate", "rep_average", "read_count_ref"])

# test_preprocessing()
test_filtering()
