import csv
import datetime
import json
import logging
import os
import shutil

import sys

import time

import re

from dartqc.FilterResult import FilterResult
from dartqc import PipelineOptions
from dartqc.Dataset import Dataset

log = logging.getLogger("Pipeline")
PARAMS_FILE = "filter_params.txt"


def read_data(working_dir: str, batch_id: str, type: str = "dart", files: [str] = None,
              unknown_args: [] = None) -> Dataset:
    """
    Read all files into the dataset structure & save to JSON

     Note: For each genotype, this should be done once only as the first step.

    :param type: Type of input files (eg. genotype provider specific)
    :param files:
    :return:
    """
    if type not in PipelineOptions.input_types:
        log.error("Invalid input/reader type: " + type + "  Available options include: " +
                  str(list(PipelineOptions.input_types.keys())))
        sys.exit(1)

    log.info("Reading input files into dataset: {}".format(files))
    dataset = PipelineOptions.input_types[type].read(working_dir, batch_id, files, unknown_args)

    return dataset


def rename_clone_ids(dataset, official_ids: str):
    """
    Rename clone IDs using the reference sequence & SNP loc.

    :return:
    """
    if not os.path.isabs(official_ids):
        official_ids = os.path.join(dataset.working_dir, official_ids)

    reference = {}
    with open(official_ids, "r") as ids_file:
        ids_file.readline()  # Skip the first line - headers

        for row in ids_file:
            parts = re.split(r"[,;\t]", row.strip())

            reference[parts[1]] = parts[0][: parts[0].find("|")].strip() if "|" in parts[0] else parts[0].strip()

    renamed = []
    for snp_def in dataset.snps:
        if snp_def.sequence_ref in reference:
            if reference[snp_def.sequence_ref] in snp_def.allele_id:
                continue  # Already has the correct name!

            old_allele_id = snp_def.allele_id
            snp_def.allele_id = reference[snp_def.sequence_ref] + snp_def.allele_id[snp_def.allele_id.find("|"):]
            snp_def.clone_id = reference[snp_def.sequence_ref]

            dataset.calls[snp_def.allele_id] = dataset.calls[old_allele_id]
            del dataset.calls[old_allele_id]

            dataset.read_counts[snp_def.allele_id] = dataset.read_counts[old_allele_id]
            del dataset.read_counts[old_allele_id]

            dataset.replicate_counts[snp_def.allele_id] = dataset.replicate_counts[old_allele_id]
            del dataset.replicate_counts[old_allele_id]

            renamed.append(old_allele_id + "->" + snp_def.allele_id)

    log.info("Renamed {} SNPs: {}".format(len(renamed), renamed[:300] + (["..."] if len(renamed) > 300 else [])))


def filter(dataset, unkown_args: [], **kwargs):
    """
    Complete all filtering

    :return:
    """
    # Default filter order is based off the order added to filter_types dictionary in Filters.py
    num_sets = 0  # Max # params in any filter (# of sets) -> pad out any filters with less params using last value
    enabled_filters = []
    for filter_type in PipelineOptions.filter_types.keys():
        if filter_type in kwargs and kwargs[filter_type] is not None and len(kwargs[filter_type]) > 0:
            enabled_filters.append(filter_type)
            if len(kwargs[filter_type]) > num_sets:
                num_sets = len(kwargs[filter_type])

    # Remove unknown filters (eg. typo) or filters which aren't actually enabled.
    filter_order = [] if "order" not in kwargs else kwargs["order"]
    for aFilter in filter_order:
        if aFilter not in enabled_filters:
            filter_order.remove(aFilter)

    # Add all missing filters to the filter order.
    for aFilter in enabled_filters:
        if aFilter not in filter_order:
            filter_order.append(aFilter)

    count = 0
    set_folders = []
    for idx in range(num_sets):
        dataset.clear_filters()

        filter_folder = None
        while filter_folder is None or os.path.exists(filter_folder):
            filter_folder = os.path.join(dataset.working_dir, "Filter_" + str(count))
            count += 1

        os.mkdir(filter_folder)
        set_folders.append(filter_folder)

        # Open file to save filter params to (document what this filter set is)
        params_file = os.path.join(filter_folder, PARAMS_FILE)
        with open(params_file, "w") as params_out:
            log.info("\n==================== Set {} ====================".format(idx + 1))

            # Dump all args to make sure anything isn't missed.
            params_out.write("Set {} run at {}\n".format(idx + 1, datetime.datetime.now()))
            params_out.write("All known args: {}\n".format(json.dumps(kwargs)))
            params_out.write("Unknown args: {}\n".format(json.dumps(unkown_args)))
            params_out.write("Filter order: {}\n\n\n".format(filter_order))

            # Do the filtering
            reuse_results = False
            for filter_idx, aFilter in enumerate(filter_order):
                # Look at previous sets completed and if exaactly the same, just copy over the filter results.
                for prev_set_num in range(idx):
                    reuse_results = True
                    for filter_name in [name for name in filter_order[:filter_idx + 1]]:
                        this_filter_threshold = kwargs[filter_name][len(kwargs[filter_name]) - 1] if len(kwargs[filter_name]) <= idx else \
                            kwargs[filter_name][idx]

                        # Get the previous sets params & pad to max num sets.
                        prev_filter_threshold = kwargs[filter_name][len(kwargs[filter_name]) - 1] if len(kwargs[filter_name]) <= prev_set_num else \
                            kwargs[filter_name][prev_set_num]

                        if this_filter_threshold != prev_filter_threshold:
                            reuse_results = False
                            break

                    # Everything matches :) - skip the filtering and just copy the results
                    if reuse_results:
                        prev_results_path = os.path.join(set_folders[prev_set_num], aFilter + ".json")
                        this_results_path = os.path.join(filter_folder, aFilter + ".json")
                        dataset.filter(FilterResult.read_json(prev_results_path))
                        shutil.copy(prev_results_path, this_results_path)

                        log.info("Skipped {} (Re-use results found from set {}): {}\n".format(aFilter, prev_set_num, this_results_path))
                        break

                if not reuse_results:
                    start = time.time()
                    log.info("Running {} filtering".format(aFilter))

                    # Get the params & pad to max num sets.
                    filter_threshold = kwargs[aFilter][len(kwargs[aFilter]) - 1] if len(kwargs[aFilter]) <= idx else \
                        kwargs[aFilter][idx]

                    # Record this filter was run and with what params.
                    # params_out.write(aFilter + ": " + str(filter_threshold) + "\n")

                    # Do the filtering
                    results = PipelineOptions.filter_types[aFilter].filter(dataset, filter_threshold, unkown_args, **kwargs)

                    # Print filtering results to log & console
                    results.log(aFilter, filter_threshold, dataset, params_out)

                    # Dump the FilterResult to the datasets folder to allow for future dynamic output generation
                    results_path = os.path.join(filter_folder, aFilter + ".json")
                    log.info("Writing {} filter results to {}".format(aFilter, results_path))
                    results.write_json(results_path)

                    # Record what was filtered in the dataset so silenced calls can be pulled out
                    # Note:  This doesn't actually modify the call data, just allows a filtered copy to be grabbed
                    dataset.filter(results)

                    # Output data at all requested points.
                    for output_type in PipelineOptions.output_types:
                        if "outputs" in kwargs and aFilter in kwargs["outputs"] and output_type \
                                in kwargs["outputs"][aFilter]:
                            PipelineOptions.output_types[output_type].write(aFilter, filter_folder, dataset, unkown_args,
                                                                            **kwargs)

                    log.info("Completed {} filter in: {}s\n".format(aFilter, time.time() - start))

            dataset.filtered.log("Final Results", None, dataset, params_out, True)

            # Write the final total filtering results to JSON (ie. everything silenced for this filtering set)
            results_path = os.path.join(filter_folder, "final.json")
            log.info("Writing final filter results for this set to {}".format(results_path))
            dataset.filtered.write_json(results_path)

            # Generate all requested output types.
            for output_type in PipelineOptions.output_types:
                if "outputs" in kwargs and "final" in kwargs["outputs"] and output_type \
                        in kwargs["outputs"]["final"]:
                    PipelineOptions.output_types[output_type].write("final", filter_folder, dataset, unkown_args,
                                                                    **kwargs)

            params_out.flush()


def output(dataset, types: [str], set: str, filter: str, unkown_args: [], **kwargs):
    folder = dataset.working_dir

    # Add in all filtered snps/samples/calls/changes up to the specified filter
    if set is not None:
        if filter is None:
            filter = "final"

        if set is not None and re.match(r"", set):
            set = "filter_" + set

        folder = os.path.join(dataset.working_dir, set)
        params_path = os.path.join(folder, PARAMS_FILE)

        # Once the folder is identified - read in the filter order to correctly add in all previous filters
        with open(params_path, "r") as params_file:
            params_file.readline()  # Skip line - set & time
            params_file.readline()  # Skip line - known args
            params_file.readline()  # Skip line - unknown args

            filters = json.loads(params_file.readline().split(":")[1])

            for aFilter in filters:
                dataset.filter(FilterResult.read_json(aFilter + ".json"))

                if aFilter == filter:
                    break
    else:
        # If this is in the parent folder, there are no results so filter means nothing
        if filter is not None:
            log.warning("Filter cannot be specified when set is not given (there are no filter results to add)")
            filter = ""

    for output_type in types:
        PipelineOptions.output_types[output_type].write(filter, folder, dataset, unkown_args, **kwargs)
