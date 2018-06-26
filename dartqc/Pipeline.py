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
              id_list: str=None, unknown_args: [] = None) -> Dataset:
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
    dataset = PipelineOptions.input_types[type].read(working_dir, batch_id, files, id_list, unknown_args)

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

    log.info("Renamed {} SNPs: {}".format(len(renamed), renamed))
    # log.info("Renamed {} SNPs: {}".format(len(renamed), renamed[:300] + (["..."] if len(renamed) > 300 else [])))


def filter(dataset, filters: [], unkown_args: [], **kwargs):
    """
    Complete all filtering

    :return:
    """

    num_sets = 0  # Max # params in any filter (# of sets) -> pad out any filters with less params using last value
    for filter_data in filters:
        if len(filter_data["params"]) > num_sets:
            num_sets = len(filter_data["params"])

    #
    # if "outputs" in kwargs and "encoding" not in kwargs["outputs"]:
    #     kwargs["outputs"]["encoding"] = "ACTG"
    # Default filter order is based off the order added to filter_types dictionary in Filters.py
    # num_sets = 0  # Max # params in any filter (# of sets) -> pad out any filters with less params using last value
    # enabled_filters = []
    # for filter_type in PipelineOptions.filter_types.keys():
    #     if filter_type in kwargs and kwargs[filter_type] is not None and len(kwargs[filter_type]) > 0:
    #         enabled_filters.append(filter_type)
    #         if len(kwargs[filter_type]) > num_sets:
    #             num_sets = len(kwargs[filter_type])
    # # Remove unknown filters (eg. typo) or filters which aren't actually enabled.
    # filter_order = [] if "order" not in kwargs else kwargs["order"]
    # for aFilter in filter_order:
    #     if aFilter not in enabled_filters:
    #         filter_order.remove(aFilter)
    #
    # # Add all missing filters to the filter order.
    # for aFilter in enabled_filters:
    #     if aFilter not in filter_order:
    #         filter_order.append(aFilter)

    # Whenever a set changes record which sets/filters cannot be reused going forward.
    reusable = []
    for set_idx in range(num_sets):
        reusable.append([])
        for filter_data in filters:
            reusable[set_idx].append([i for i in range(set_idx) if i != set_idx])

    count = 0
    set_folders = []
    for set_idx in range(num_sets):
        dataset.clear_filters()

        # Find an unused filter folder & create it
        filter_folder = None
        while filter_folder is None or os.path.exists(filter_folder):
            filter_folder = os.path.join(dataset.working_dir, "Filter_" + str(count))
            count += 1

        os.mkdir(filter_folder)
        set_folders.append(filter_folder)

        set_summary = "FILTER,PARAMS,SILENCED SAMPLES,SILENCED SNPS,SILENCED CALLS,CHANGED CALLS\n"

        # Open file to save filter params to (document what this filter set is)
        params_file = os.path.join(filter_folder, PARAMS_FILE)
        with open(params_file, "w") as params_out:
            log.info("\n==================== Set {}: {} ====================".format(set_idx + 1, filter_folder))

            # Dump all args to make sure anything isn't missed.
            params_out.write("Set {} run at {}\n".format(set_idx + 1, datetime.datetime.now()))
            params_out.write("Args: {}\n".format(",".join(sys.argv)))
            params_out.write("Unknown args: {}\n".format(json.dumps(unkown_args)))
            params_out.write("Filters: {}\n\n\n".format([filter_data["name"] for filter_data in filters]))

            # Do the filtering
            for filter_idx, filter_data in enumerate(filters):
                if len(dataset.filtered.samples) == len(dataset.samples) or len(dataset.filtered.snps) == len(dataset.snps):
                    log.info("All data has been silenced - skipping all further filters")
                    break


                # Find the threshold/filtering prarms
                filter_threshold = filter_data["params"][len(filter_data["params"]) - 1] if len(filter_data["params"]) <= set_idx else filter_data["params"][set_idx]

                set_summary += '{},"{}"'.format(filter_data["name"], filter_threshold)

                # Update the re-usability (only check things that are still reusable)
                for aSet in range(num_sets):
                    if aSet != set_idx and aSet in reusable[set_idx][filter_idx]:
                        other_thresh = filter_data["params"][len(filter_data["params"]) - 1] if len(filter_data["params"]) <= aSet else filter_data["params"][aSet]

                        # If the params don't match they are not re-usable going forward
                        if filter_threshold != other_thresh:
                            for i in range(filter_idx, len(filters)):
                                if aSet in reusable[set_idx][i]:
                                    reusable[set_idx][i].remove(aSet)   # Remove re-usability of aSet from this set

                                if set_idx in reusable[aSet][i]:
                                    reusable[aSet][i].remove(set_idx)   # Remove re-usability of this set by aSet

                filter_name = filter_data["name"]
                count = 0
                for idx in range(filter_idx):
                    if filters[idx]["name"] == filter_data["name"]:
                        count += 1

                if count > 0:
                    filter_name = filter_name + str(count)

                filter_data["final_name"] = filter_name

                if len(reusable[set_idx][filter_idx]) > 0:
                    # A previously run filter is reusable!  Just copy the results
                    prev_set_num = reusable[set_idx][filter_idx][0]
                    prev_results_path = os.path.join(set_folders[prev_set_num], filter_data["final_name"] + ".json")
                    this_results_path = os.path.join(filter_folder, filter_name + ".json")
                    results = FilterResult.read_json(prev_results_path)
                    dataset.filter(results)
                    shutil.copy(prev_results_path, this_results_path)

                    log.info("Skipped {} (Re-use results found from set {}): {}\n".format(filter_name, prev_set_num, this_results_path))
                else:
                    start = time.time()
                    log.info("Running {} filtering with {}".format(filter_name, filter_threshold))

                    # Do the filtering
                    results = PipelineOptions.filter_types[filter_data["name"]].filter(dataset, filter_threshold, unkown_args, **kwargs)

                    # Print filtering results to log & console
                    results.log(filter_name, filter_threshold, dataset, params_out)

                    # Dump the FilterResult to the datasets folder to allow for future dynamic output generation
                    results_path = os.path.join(filter_folder, filter_name + ".json")
                    log.info("Writing {} filter results to {}".format(filter_name, results_path))
                    results.write_json(results_path)

                    # Record what was filtered in the dataset so silenced calls can be pulled out
                    # Note:  This doesn't actually modify the call data, just allows a filtered copy to be grabbed
                    dataset.filter(results)

                if len(dataset.filtered.samples) == len(dataset.samples) or len(dataset.filtered.snps) == len(dataset.snps):
                    log.info("All data has been silenced - skipping outputs")
                else:
                    for output_type, encoding in filter_data["outputs"].items():
                        PipelineOptions.output_types[output_type].write(filter_name, filter_folder, encoding, dataset, unkown_args, **kwargs)

                log.info("Completed {} filter in: {:.2f}s\n".format(filter_name, time.time() - start))

                set_summary += ",{},{},{},{}\n".format(len(results.samples), len(results.snps), results.get_tot_silenced_calls(), results.get_tot_changed_calls())

            dataset.filtered.log("Final Results", None, dataset, params_out, True)

            # Write the final total filtering results to JSON (ie. everything silenced for this filtering set)
            results_path = os.path.join(filter_folder, "final.json")
            log.info("Writing final filter results for this set to {}".format(results_path))
            dataset.filtered.write_json(results_path)

            # # Generate all requested output types.
            # for output_type in PipelineOptions.output_types:
            #     if "outputs" in kwargs and "final" in kwargs["outputs"] and output_type \
            #             in kwargs["outputs"]["final"]:
            #         PipelineOptions.output_types[output_type].write("final", filter_folder,
            #                                                         kwargs["outputs"]["encoding"], dataset, unkown_args,
            #                                                         **kwargs)

            params_out.flush()

            set_summary += "TOTAL SILENCED,,{},{},{},{}\n".format(len(dataset.filtered.samples), len(dataset.filtered.snps), dataset.filtered.get_tot_silenced_calls(), dataset.filtered.get_tot_changed_calls())
            set_summary += "TOTAL IN DATASET,,{},{},{},{}".format(len(dataset.samples), len(dataset.snps), len(dataset.samples) * len(dataset.snps), len(dataset.samples) * len(dataset.snps))

            csv_file = os.path.join(filter_folder, "filter_summary.csv")
            with open(csv_file, "w") as csv_out:
                csv_out.write(set_summary)
                csv_out.flush()


def output(dataset, types: [str], set: str, filter: str, unknown_args: [], **kwargs):
    folder = dataset.working_dir

    # Add in all filtered snps/samples/calls/changes up to the specified filter
    if set is not None:
        if filter is None:
            filter = "final"

        if set is not None and not set.lower().startswith("filter_"):
            set = "Filter_" + set

        folder = os.path.join(dataset.working_dir, set)
        params_path = os.path.join(folder, PARAMS_FILE)

        # Once the folder is identified - read in the filter order to correctly add in all previous filters
        try:
            with open(params_path, "r") as params_file:
                params_file.readline()  # Skip line - set & time
                params_file.readline()  # Skip line - known args
                params_file.readline()  # Skip line - unknown args

                filters = params_file.readline().split(":")[1].strip()
                filters = re.sub(r"[\[\]\']", "", filters).split(",")

                filter_cnts = {}

                for aFilter in filters:
                    filter_name = aFilter.strip()
                    if aFilter in filter_cnts:
                        filter_cnts[aFilter] += 1
                        filter_name += filter_cnts[aFilter]
                    else:
                        filter_cnts[aFilter] = 0

                    aFilter_path = os.path.join(folder, filter_name + ".json")
                    dataset.filter(FilterResult.read_json(aFilter_path))

                    if filter_name == filter:
                        break
        except FileNotFoundError as ex:
            log.error("Incorrect working dir or set: {}".format(ex))
            return
    else:
        # If this is in the parent folder, there are no results so filter means nothing
        if filter is not None:
            log.warning("Filter cannot be specified when set is not given (there are no filter results to add)")
            filter = ""

    for output_type in types:
        encoding = "ACTG"
        if ":" in output_type:
            output_type, encoding = output_type.split(":")

        PipelineOptions.output_types[output_type].write(filter, folder, encoding, dataset, unknown_args, **kwargs)
