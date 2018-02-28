import os
import time
import argparse
import textwrap

from subprocess import call, check_output, CalledProcessError
import logging

import re

import PipelineOptions


class CmdLine:
    def __init__(self):
        parser = argparse.ArgumentParser()

        parser.add_argument("--version", "-v", dest="version", action="store_true", help="print version and exit")

        parser.add_argument("--list_inputs", dest="list_inputs", action="store_true",
                            help="print list of implemented data input options and exit (supported genotype providers)")

        parser.add_argument("--list_filters", dest="list_filters", action="store_true",
                            help="print list of implemented filter types")

        parser.add_argument("--list_outputs", dest="list_outputs", action="store_true",
                            help="print list of implemented output data types")

        parser.add_argument("--working_dir", "-w", type=lambda p: os.path.abspath(p),
                            default=os.getcwd(), required=False, dest="working_dir",
                            help="Root folder for all filtering files (defaults to current dir)")

        parser.add_argument("--batch_id", "-b", type=str, default="filter", required=False,
                            dest="batch_id", help="Genotype batch ID - used as prefix for files")

        subparsers = parser.add_subparsers(help='Command-line interface for DartQC')

        prepare_parser = subparsers.add_parser("read")

        prepare_parser.add_argument("--type", "-t", type=str, default="dart", required=False,
                                    dest="type", help="Type of input data (eg. genotype provider).  "
                                                      "Use --print_providers for list of options")

        prepare_parser.add_argument("--files", "-f", type=lambda s: [p for p in s.split(",")],
                                    required=True, dest="files", help="path to input file(s)")

        prepare_parser.set_defaults(subparser='read')

        validate_parser = subparsers.add_parser("validate")

        validate_parser.add_argument("--id_list", "-i", default="id_list.csv", type=lambda p: os.path.abspath(p),
                                     required=True, dest="id_list",
                                     help="path to CSV file with list of official clone IDs that should be used (eg. to fix ID's that Dart outputs wrong)")

        validate_parser.add_argument("--identity", default=0.95, type=float,
                                     dest="identity", help="remove snps in identical sequence clusters")

        validate_parser.set_defaults(subparser='validate')

        filter_parser = subparsers.add_parser("filter")

        # Add all configured filter types
        for filter_name, filter in PipelineOptions.filter_types.items():
            filter_parser.add_argument("--" + filter_name, dest=filter_name,
                                       type=filter.get_cmd_type(), help=filter.get_cmd_help())

        # # Input values must follow [<val>, <val>, ...] - missing values should be empty as in [,,,]
        filter_parser.add_argument("--ignore_pops", "-i", default=False, dest="ignore_pops",
                                   type=bool, help="Don't use populations (originally intended for use with MAF & HWE)")

        filter_parser.add_argument("--pop_blacklist", default=[], dest="pop_blacklist",
                                   type=lambda s: [item for item in re.sub(r"[\[\] ]", "", s).split(",")],
                                   help="Populations to ignore (originally intended for use with MAF & HWE)")
        # filter_parser.add_argument("--maf", default=[], dest="minor_allele_freq",
        #                            type=lambda s: [float(item.strip()) if len(item.strip()) > 0 else None
        #                                            for item in s[1:-1].split(',')],
        #                            help="filter snps <= minor allele frequency")
        # filter_parser.add_argument("--hwe", default=[],
        #                            type=lambda s: [float(item.strip()) if len(item.strip()) > 0 else None for item in
        #                                            s[1:-1].split(',')],
        #                            dest="hwe", help="filter snps <= p-value of hardy-weinberg test")
        # filter_parser.add_argument("--call_rate", default=[],
        #                            type=lambda s: [float(item.strip()) if len(item.strip()) > 0 else None for item in
        #                                            s[1:-1].split(',')], dest="call_rate",
        #                            help="filter snps <= call rate of snp")
        # filter_parser.add_argument("--rep", default=[],
        #                            type=lambda s: [float(item.strip()) if len(item.strip()) > 0 else None for item in
        #                                            s[1:-1].split(',')],
        #                            dest="rep", help="filter snps <= replication average of snp")
        #
        # filter_parser.add_argument("--mind", default=[],
        #                            type=lambda s: [float(item.strip()) if len(item.strip()) > 0 else None for item in
        #                                            s[1:-1].split(',')],
        #                            dest="mind", help="filter samples > missingness per sample")
        # filter_parser.add_argument("--mono", default=None,
        #                            dest="mono", help="filter samples monomorphic in <mono> populations ('all', int)")
        # filter_parser.add_argument("--mono_comparison", default="==",
        #                            dest="mono_comp", help="filter samples monomorphic in >=, <=, == populations ('==')")
        #
        # filter_parser.add_argument("--clusters", default=False, action="store_true",
        #                            dest="remove_clusters", help="remove snps in identical sequence clusters")
        # filter_parser.add_argument("--identity", default=0.95, type=float,
        #                            dest="identity", help="remove snps in identical sequence clusters")

        filter_parser.add_argument("--outputs", default={},
                                   type=lambda s: {k: v.split(",") for item in re.split(r"],", re.sub(r"[{} ]", "", s))
                                                   for k, v in dict([re.sub(r"[\[\] ]", "", item).split(":")]).items()},
                                   required=False,
                                   dest="outputs", help="Which filters should generate genotypes")

        filter_parser.add_argument("--order", "-o", default=[],
                                   type=lambda s: [item for item in s.split(',')],
                                   required=False, dest="data", help="Which filters should generate ped & map files")

        filter_parser.set_defaults(subparser='filter')

        output_parser = subparsers.add_parser("output",
                                              help="Primarily for user interface to allow generation of output types on click")

        validate_parser.add_argument("--types", type=lambda s: [type for type in s.split(",")],
                                     required=True, dest="types",
                                     help="Output type(s) to generate")

        validate_parser.add_argument("--folder", type=str, default=None,
                                     dest="folder",
                                     help="Dataset folder (eg. Filter_1) - leave blank for unfiltered dataset")

        validate_parser.add_argument("--filter", type=str, default=None,
                                     dest="filter",
                                     help="FilterResult to generate output data from - leave blank for final results")

        output_parser.set_defaults(subparser='output')

        self.args, self.unknown = parser.parse_known_args()


def stamp(*args):
    message = str(time.strftime("[%H:%M:%S]")) + " " + " ".join([str(arg) for arg in args])
    print(message)
    logging.info(message)
