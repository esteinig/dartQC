import csv
import os

import logging
import numpy

from Dataset import Dataset
from PipelineOptions import Output

log = logging.getLogger(__file__)


class PlinkOutput(Output):
    def get_name(self) -> str:
        return "plink"

    def get_description(self) -> str:
        return "Output as PED and MAP files"

    def write(self, filter_name: str, folder: str, dataset: Dataset, unknown_args: [], **kwargs) -> None:
        file_path = os.path.join(folder, filter_name + "_" + self.get_name())

        ped_file = file_path + '.ped'
        map_file = file_path + '.map'

        log.info("Writing MAP file")

        # Write .map file - seems like its just a listing of allele ID's?
        map_data = [["0", snp_def, "0", "0"] for snp_def in dataset.snps]
        with open(map_file, 'w') as map_out:
            ped_writer = csv.writer(map_out, delimiter="\t")
            ped_writer.writerows(map_data)

        filtered_calls = dataset.get_filtered_calls()

        # Convert from (0,1) type tuples to ACTG strings (eg. AA, TA, etc)
        snp_call_matrix = []
        for allele_id, snp_calls in filtered_calls.items():
            snp_vals = dataset.get_snp_def(allele_id).snp.split(":")[1].split(">")

            # Convert to ACTG values such as AC, TG, etc. instead of 0,1 or 1,0 etc.
            calls = numpy.asarray([snp_vals[0] + snp_vals[1] if tuple(call) == Dataset.heterozygous
                     else snp_vals[0] + snp_vals[0] if tuple(call) == Dataset.homozygous_major
            else snp_vals[1] + snp_vals[1] if tuple(call) == Dataset.homozygous_minor
            else "00" for call in snp_calls])
            snp_call_matrix.append(calls)

        sample_call_matrix = numpy.dstack(snp_call_matrix)[0]

        plink_rows = []
        for idx, sample_def in enumerate(dataset.samples):
            ped_row = [sample_def.population, sample_def.id, "0", "0", "0", "0"]
            row_calls = [call_letter for call_str in sample_call_matrix[idx] for call_letter in call_str]

            ped_row.extend(row_calls)
            plink_rows.append(ped_row)

        log.info("Writing PED file")
        with open(ped_file, 'w') as ped_out:
            for row in plink_rows:
                ped_out.write("\t".join(row) + "\n")
                ped_out.flush()

                # MAP Formatting


PlinkOutput()
