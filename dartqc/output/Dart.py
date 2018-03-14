import csv
import os

import logging

import numpy

from dartqc.Dataset import Dataset
from dartqc.PipelineOptions import Output

log = logging.getLogger(__file__)


class CSVOutput(Output):
    def get_name(self) -> str:
        return "dart"

    def get_description(self) -> str:
        return "Output the calls in DArT format (2 rows per SNP SNP per row, 1 sample per column)"

    def write(self, filter_name: str, folder: str, encoding: str, dataset: Dataset, unknown_args: [], **kwargs) -> None:
        file_path = os.path.join(folder, filter_name + "_" + self.get_name() + ".csv")

        log.info("Outputting DArT Calls")

        if encoding is not None and encoding != "11":
            log.warning("Only the 11 encoding is supported for the Dart output type - {} ignored".format(encoding))

        filtered_calls = dataset.get_filtered_calls()

        # Get calls as matrix - swap SNP/sample axes & trim back to only first allele's call (much quicker processing)
        numpy_matrix = numpy.asarray([filtered_calls[snp.allele_id] for snp in dataset.snps])
        numpy_matrix = numpy.dstack(numpy_matrix)  # [SNPs][samples][calls] -> [samples][calls][SNPs]
        numpy_matrix = numpy.dstack(numpy_matrix)  # [SNPs][calls][samples] -> [calls][SNPs][samples]
        numpy_matrix = numpy_matrix  # Only get first allele calls (if "-" -> missing)

        del filtered_calls

        with open(file_path, "w") as file_out:
            writer = csv.writer(file_out, lineterminator='\n')

            pops_row = [""] * len(dataset.snps[0].all_headers) + [sample_def.population for sample_def in dataset.samples]
            samples_row = list(dataset.snps[0].all_headers.keys()) + [sample_def.id for sample_def in dataset.samples]

            writer.writerow(pops_row)
            writer.writerow(samples_row)

            for snp_idx, snp_def in enumerate(dataset.snps):
                snp_headers = list(snp_def.all_headers.values())

                writer.writerow(snp_headers + numpy_matrix[0][snp_idx].tolist())

                # Replace the first rows sequence (ref seq.) with the second rows sequence
                if "AlleleSequence" in snp_def.all_headers:
                    snp_headers[list(snp_def.all_headers.keys()).index("AlleleSequence")] = snp_def.sequence

                writer.writerow(snp_headers + numpy_matrix[1][snp_idx].tolist())

            file_out.flush()


CSVOutput()
