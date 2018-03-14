#!/usr/bin/env python
import csv
import datetime
import glob
import inspect

import numpy
import sys
import os

import time

from os.path import dirname, basename, isfile

import re

# Need these unused imports to load the __init__.py scripts which auto load the extension points
from dartqc import filters
from dartqc import output
from dartqc import input

from dartqc import PipelineOptions
from dartqc import Pipeline
from dartqc.CmdLine import CmdLine
import logging
import logging.handlers

from dartqc.Dataset import Dataset
from install import cdhit_config


# file = "C:\\Users\\jc229483\\Desktop\\ListIndexErr\\subsetted__recount_file_2.csv"
#
# with open(file, "r") as test_in:
#     reader = csv.reader(test_in)
#
#     last_clone = None
#     row_index = 0
#     for row in reader:
#         if row_index % 5000 == 0:
#             print("row {}".format(row_index))
#
#         if any(row) and row_index >= 6:
#             if last_clone is None:
#                 last_clone = row[1]
#             else:
#                 if last_clone != row[1]:
#                     test = 2
#
#                 last_clone = None
#
#             for val in row[32:]:
#                 test = int(val)
#
#         row_index += 1
#
# exit(0)


def main():
    start = time.time()

    cmd_line = CmdLine()

    args = vars(cmd_line.args)
    unknown_args = cmd_line.unknown

    if args["version"]:
        print("0.1.5")
        exit(0)

    if args["list_inputs"]:
        for type in PipelineOptions.input_types.values():
            print("{}: {}".format(type.get_name(), type.get_description()))
        exit(0)

    if args["list_filters"]:
        for type in PipelineOptions.filter_types.values():
            print("{} ({}): {}".format(type.get_name(), ",".join(type.get_cmd_names()), type.get_description()))
        exit(0)

    if args["list_outputs"]:
        for type in PipelineOptions.output_types.values():
            print("{}: {}".format(type.get_name(), type.get_description()))
        exit(0)

    # Clean up args so they aren't output later...
    del args["version"]
    del args["list_inputs"]
    del args["list_filters"]
    del args["list_outputs"]

    working_dir = args["working_dir"]
    batch_id = args["batch_id"]

    # Make sure the working dir exists before we start
    os.makedirs(working_dir, exist_ok=True)

    setup_logging(working_dir, batch_id)
    log = logging.getLogger(__file__)

    log.info("Filtering started at " + datetime.date.strftime(datetime.date.today(), "%d/%m/%y")
                 + "\nArgs: " + str(sys.argv) + "\n"
             + "Unknown args: " + str(unknown_args) + "\n")

    # Workflows
    dataset = None
    if args["subparser"] == "read":
        dataset = Pipeline.read_data(working_dir, batch_id, args["type"], args["files"], unknown_args)

    # Load the dataset from json file (must have been read in already)
    dataset_path = os.path.join(working_dir, batch_id + ".npy")
    if dataset is None:
        # Datasets are always located in the same place with the same name.
        # This follows the pattern <working_dir>/<batch_id>.json
        if not os.path.exists(dataset_path):
            logging.error("\n\nMissing dataset!  Generate using read option first\n"
                          "(Read parses the genotype calls & read counts into a single data structure and saves as json\n\n"
                          "This file must exist: " + dataset_path + "\n\n")
            exit(1)

        dataset = Dataset.read_json(dataset_path)

    if args["subparser"] == "rename":
        Pipeline.rename_clone_ids(dataset, args["id_list"])

    if args["subparser"] == "read" or args["subparser"] == "rename":
        log.info("Writing dataset to file: {}".format(dataset_path))
        dataset.write_json(dataset_path)

    if args["subparser"] == "filter":
        Pipeline.filter(dataset, unknown_args, **args)

    if args["subparser"] == "output":
        Pipeline.output(dataset, unknown_args=unknown_args, **args)

    log.debug("\nRun time: " + str(round((time.time() - start), 2)) + "s")

def setup_logging(working_dir, batch_id):
    logger = logging.getLogger("")

    logger.setLevel(logging.DEBUG)

    log_file = os.path.join(working_dir, batch_id + ".log")

    # Log format
    log_formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%H:%M:%S")

    # File logs
    file_handler = logging.handlers.RotatingFileHandler(log_file, maxBytes=(1048576 * 5), backupCount=7)
    file_handler.setFormatter(log_formatter)
    logger.addHandler(file_handler)

    # Console logs
    # Stderr
    consoleErrHandler = logging.StreamHandler(sys.stderr)
    consoleErrHandler.setFormatter(log_formatter)
    consoleErrHandler.setLevel(logging.WARNING)
    logger.addHandler(consoleErrHandler)

    # Stdout
    class InfoFilter(logging.Filter):
        def filter(self, rec):
            return rec.levelno in (logging.DEBUG, logging.INFO)

    consoleOutHandler = logging.StreamHandler(sys.stdout)
    consoleOutHandler.setFormatter(log_formatter)
    consoleOutHandler.addFilter(InfoFilter())
    consoleOutHandler.setLevel(logging.DEBUG)
    logger.addHandler(consoleOutHandler)



main()
