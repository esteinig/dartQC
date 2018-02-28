import csv
import os

import logging

from Dataset import Dataset
from PipelineOptions import Output

log = logging.getLogger(__file__)


class GraphOutput(Output):
    def get_name(self) -> str:
        return "graph"

    def get_description(self) -> str:
        return "Generate a set of graphs as a summary of the dataset filtering changes"

    def write(self, filter_name: str, folder: str, dataset: Dataset, unknown_args: [], **kwargs) -> None:
        file_path = os.path.join(folder, filter_name + "_" + self.get_name() + ".csv")

        log.info("Outputting Graphs")

        with open(file_path, "w") as csv_out:
            headers = [""] + [sample_def.id for sample_def in dataset.samples]
            csv_out.write(",".join(headers) + "\n")

            for allele_id, snp_calls in dataset.get_filtered_calls().items():
                snp_vals = dataset.get_snp_def(allele_id).snp.split(":")[1].split(">")
                calls = [snp_vals[0] + snp_vals[1] if call == Dataset.heterozygous
                         else snp_vals[0] + snp_vals[0] if call == Dataset.homozygous_major
                         else snp_vals[1] + snp_vals[1] if call == Dataset.homozygous_minor
                         else "--" for call in snp_calls]

                csv_out.write("{},{}\n".format(allele_id, ",".join(calls)))

            csv_out.flush()


GraphOutput()
