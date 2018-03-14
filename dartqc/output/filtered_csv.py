import csv
import os

import logging

from dartqc.Dataset import Dataset
from dartqc.PipelineOptions import Output

log = logging.getLogger(__file__)


class FilteredCSVOutput(Output):
    def get_name(self) -> str:
        return "filtered"

    def get_description(self) -> str:
        return "Output the genotype in CSV format (1 SNP per row, 1 sample per column)"

    def write(self, filter_name: str, folder: str, dataset: Dataset, unknown_args: [], **kwargs) -> None:
        file_path = os.path.join(folder, filter_name + "_" + self.get_name() + ".csv")

        log.info("Outputting CSV")

        with open(file_path, "w") as csv_out:
            headers = [""] + [sample_def.id for sample_def in dataset.samples]
            csv_out.write(",".join(headers) + "\n")

            for allele_id, snp_calls in dataset.get_filtered_calls().items():
                snp_vals = dataset.get_snp_def(allele_id).snp.split(":")[1].split(">")
                calls = [snp_vals[0] + snp_vals[1] if tuple(call) == Dataset.heterozygous
                         else snp_vals[0] + snp_vals[0] if tuple(call) == Dataset.homozygous_major
                else snp_vals[1] + snp_vals[1] if tuple(call) == Dataset.homozygous_minor
                else "--" for call in snp_calls]

                csv_out.write("{},{}\n".format(allele_id, ",".join(calls)))

            csv_out.flush()


FilteredCSVOutput()