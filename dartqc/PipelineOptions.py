from collections import OrderedDict

import Dataset
import FilterResult


# from filters import HardyWeinberg
# from filters import HetComp
# from filters import MinSNPData
# from filters import MinSampleData
# from filters import MinorAlleleFreq
# # from filters import SNPMetric
# from filters import ReadCounts
# from filters import Cluster
# from input import DartReader
# from output import CSV
# from output import Plink
#
# # All filter types available
# # Order matters!  This is the default filter order
# # Each filter should have the signature:
# #       def filter_name(dataset: Dataset, threshold: ??? typ. float, unkown_args: [???], **kwargs)
# #       unknown_args & kwargs are provided so 3rd party implementations can access cmd line args.
# #       threshold is implementation specific based on the cmd line parsing, but it is typically a float.
# filter_types = OrderedDict()
# filter_types[ReadCounts.NAME] = ReadCounts.filter
# # filter_types[SNPMetric.NAME] = SNPMetric.filter
# filter_types[HetComp.NAME] = HetComp.filter
# filter_types[MinSampleData.NAME] = MinSampleData.filter
# filter_types[MinSNPData.NAME] = MinSNPData.filter
# filter_types[Cluster.NAME] = Cluster.filter
# filter_types[MinorAlleleFreq.NAME] = MinorAlleleFreq.filter
# filter_types[HardyWeinberg.NAME] = HardyWeinberg.filter
#
# # All available input types (eg. supported genotype providers)
# input_types = {
#     DartReader.NAME: DartReader.read
# }
#
# # All available output types
# output_types = {
#     CSV.NAME: CSV.output,
#     Plink.NAME: Plink.output
# }


class Filter:
    def __init__(self):
        add_filter(self)

    def get_name(self) -> str:
        raise NotImplemented()

    def get_cmd_type(self):
        raise NotImplemented()

    def get_cmd_help(self) -> str:
        raise NotImplemented()

    def get_order(self) -> int:
        raise NotImplemented()

    def filter(self, dataset: Dataset, threshold, unknown_args: [], **kwargs) -> FilterResult:
        raise NotImplemented()


class Input:
    def __init__(self):
        add_input(self)

    def get_name(self) -> str:
        raise NotImplemented()

    def get_description(self) -> str:
        raise NotImplemented()

    def read(self, working_dir: str, batch_id: str, files: [str], unknown_args: [] = None, **kwargs) -> Dataset:
        raise NotImplemented()


class Output:
    def __init__(self):
        add_output(self)

    def get_name(self) -> str:
        raise NotImplemented()

    def get_description(self) -> str:
        raise NotImplemented()

    def write(self, filter_name: str, folder: str, dataset: Dataset, unknown_args: [], **kwargs) -> None:
        raise NotImplemented()


input_types = {}
output_types = {}
filter_types = {}


def add_filter(filter: Filter):
    filter_types[filter.get_name()] = filter

    sorted_names = sorted(filter_types, key=lambda o: filter_types[o].get_order())
    temp_filters = dict(filter_types)

    filter_types.clear()
    for filter_name in sorted_names:
        filter_types[filter_name] = temp_filters[filter_name]


def add_input(input: Input):
    input_types[input.get_name()] = input


def add_output(output: Output):
    output_types[output.get_name()] = output
