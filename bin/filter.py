#!/usr/bin/env python
import datetime
import glob
import inspect

import numpy
import sys
import os

import time

from os.path import dirname, basename, isfile

import PipelineOptions
import Pipeline
from CmdLine import CmdLine
import logging
import logging.handlers

from Dataset import Dataset
from install import cdhit_config

import filters
import output
import input


def main():
    start = time.time()

    cmd_line = CmdLine()

    args = vars(cmd_line.args)
    unknown_args = cmd_line.unknown

    if args["version"]:
        print("0.1.5")
        exit(0)

    if args["list_inputs"]:
        print(list(PipelineOptions.input_types.keys()))
        exit(0)

    if args["list_filters"]:
        print(list(PipelineOptions.filter_types.keys()))
        exit(0)

    if args["list_outputs"]:
        print(list(PipelineOptions.filter_types.keys()))
        exit(0)

    # Clean up args so they aren't output later...
    del args["version"]
    del args["list_inputs"]
    del args["list_filters"]
    del args["list_outputs"]


    working_dir = args["working_dir"]
    batch_id = args["batch_id"]

    # Make sure the working dir exists before we start
    os.makedirs(args["working_dir"], exist_ok=True)

    setup_logging(working_dir, batch_id)
    log = logging.getLogger(__file__)

    log.info("Filtering started at " + datetime.date.strftime(datetime.date.today(), "%d/%m/%y")
                 + "\nArgs: " + str(sys.argv) + "\n"
             + "Unknown args: " + str(unknown_args) + "\n")

    # Datasets are always located in the same place with the same name.
    # This follows the pattern <working_dir>/<batch_id>.json
    dataset_path = os.path.join(working_dir, batch_id + ".json")

    # Workflows
    if args["subparser"] == "read":
        dataset = Pipeline.read_data(working_dir, batch_id, args["type"], args["files"], unknown_args)
        # log.debug("\nRead time: " + str(round((time.time() - start), 2)) + "s")

        # json_start = time.time()
        dataset.write_json(dataset_path)
        # log.debug("\nJSON write time: " + str(round((time.time() - json_start), 2)) + "s")

    # Load the dataset from json file (must have been read in already)
    else:
        if not os.path.exists(dataset_path):
            logging.error("\n\nMissing dataset!  Generate using read option first\n"
                          "(Read parses the genotype calls & read counts into a single data structure and saves as json\n\n"
                          "This file must exist: " + dataset_path + "\n\n")
            exit(1)

        dataset = Dataset.read_json(dataset_path)
        # log.debug("\nJSON read time: " + str(round((time.time() - start), 2)) + "s")

        if args["subparser"] == "validate":
            Pipeline.fix_clone_ids(dataset, args["id_list"], args["identity"])

        if args["subparser"] == "filter":
            Pipeline.filter(dataset, unknown_args, **args)

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
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(log_formatter)
    logger.addHandler(consoleHandler)



main()
