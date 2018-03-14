from collections import OrderedDict

import re

from dartqc.Dataset import Dataset
from dartqc.FilterResult import FilterResult


class Filter:
    LIST_OF_FLOAT = lambda s: [float(item) if len(item) > 0 else None for item in re.sub(r'[()\[\] ]', "", s).split(",")]
    LIST_OF_LISTS = lambda s: [re.sub(r'[()\[\] ]', "", item).split(",") for item in re.split(r"[)\]],[(\[]", re.sub(r'(^[(\[]{2})|([)\]]{2}$)', r"\1", s)) if len(s.strip()) > 0]

    def __init__(self):
        add_filter(self)

    def get_name(self) -> str:
        raise NotImplemented()

    def get_description(self) -> str:
        raise NotImplemented()

    def get_cmd_names(self) -> [str]:
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

    def write(self, filter_name: str, folder: str, encoding: str, dataset: Dataset, unknown_args: [], **kwargs) -> None:
        raise NotImplemented()

if "input_types" not in globals():
    global input_types
    input_types = {}

    global output_types
    output_types = {}

    global filter_types
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
