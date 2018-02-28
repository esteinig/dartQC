import datetime
import json
import logging
import os

import sys

import time

import PipelineOptions
from AlleleIDValidation import IDValidator
from Dataset import Dataset

log = logging.getLogger("Pipeline")


def read_data(working_dir:str, batch_id:str, type: str = "dart", files: [str] = None, unknown_args: []=None) -> Dataset:
    """
    Read all files into the dataset structure - this should be done once only as the first step.

    :param type: Type of input files (eg. genotype provider specific)
    :param files:
    :return:
    """

    if type not in PipelineOptions.input_types:
        log.error("Invalid input/reader type: " + type + "  Available options include: " + str(list(PipelineOptions.input_types.keys())))
        sys.exit(1)

    return PipelineOptions.input_types[type].read(working_dir, batch_id, files, unknown_args)


def fix_clone_ids(dataset, official_ids: str, dist: float):
    """
    Use cd-hit-est-2d to cluster SNPs within dist (0-1) using a provided list of official seqences and rename
    any clustered SNPs to the official ID.

    Note:  Does not remove SNPs and has no effect at all on the filtering - this is a tool to help downstream processing

    :return:
    """
    # TODO:  Validate IDs
    # validator = IDValidator(dataset, official_ids, dist)
    # validator.do_validations()
    pass


def filter(dataset, unkown_args:[], **kwargs):
    """
    Complete all filtering

    :return:
    """
    # Default filter order is based off the order added to filter_types dictionary in Filters.py
    num_sets = 0    # Max # params in any filter (# of sets) -> pad out any filters with less params using last value
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
    for idx in range(num_sets):
        dataset.clear_filters()

        filter_folder = None
        while filter_folder is None or os.path.exists(filter_folder):
            filter_folder = os.path.join(dataset.working_dir, "Filter_" + str(count))
            count += 1

        os.mkdir(filter_folder)

        # Open file to save filter params to (document what this filter set is)
        params_file = os.path.join(filter_folder, "filter_params.txt")
        with open(params_file, "w") as params_out:
            log.info("\n==================== Set {} ====================".format(idx + 1))

            # Dump all args to make sure anything isn't missed.
            params_out.write("Set {} run at {}\n".format(idx + 1, datetime.datetime.now()))
            params_out.write("All known args: {}\n".format(json.dumps(kwargs)))
            params_out.write("Unknown args: {}\n\n\n".format(json.dumps(unkown_args)))

            # Do the filtering
            for aFilter in filter_order:
                # Get the params & pad to max num sets.
                filter_threshold = kwargs[aFilter][len(kwargs[aFilter]) - 1] if len(kwargs[aFilter]) <= idx else kwargs[aFilter][idx]

                # Record this filter was run and with what params.
                # params_out.write(aFilter + ": " + str(filter_threshold) + "\n")

                # Do the filtering
                results = PipelineOptions.filter_types[aFilter].filter(dataset, filter_threshold, unkown_args, **kwargs)

                # Print filtering results to log & console
                results.log(aFilter, filter_threshold, dataset, params_out)

                # Record what was filtered in the dataset so silenced calls can be pulled out
                # Note:  This doesn't actually modify the call data, just allows a filtered copy to be grabbed
                dataset.filter(results)

                # Output data at all requested points.
                for output_type in PipelineOptions.output_types:
                    if "outputs" in kwargs and aFilter in kwargs["outputs"] and output_type in kwargs["outputs"][aFilter]:
                        PipelineOptions.output_types[output_type].write(aFilter, filter_folder, dataset, unkown_args, **kwargs)

                # Dump the FilterResult to the datasets folder to allow for future dynamic output generation
                results_path = os.path.join(filter_folder, aFilter + ".json")
                results.write_json(results_path)

            # Write the final total filtering results to JSON (ie. everything silenced for this filtering set)
            results_path = os.path.join(filter_folder, "final.json")
            dataset.filtered.write_json(results_path)

            dataset.filtered.log("Final Results", None, dataset, params_out, True)

            filtered_calls = dataset.get_filtered_calls()
            params_out.write("\n\nFinal data:\n{}".format("\n\t".join(["{}: {}".format(allele_id, snp_calls) for allele_id, snp_calls in filtered_calls.items()])))

            params_out.flush()

