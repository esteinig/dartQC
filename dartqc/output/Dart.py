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
        return "Output the calls in DArT format (2 rows per SNP, 1 sample per column) - ready to use as input files"

    def write(self, filter_name: str, folder: str, encoding: str, dataset: Dataset, unknown_args: [], **kwargs) -> None:
        file_path = os.path.join(folder, filter_name + "_" + self.get_name() + ".csv")

        log.info("Outputting DArT Calls")

        if encoding is not None and encoding != "11":
            log.warning("Only the 11 encoding is supported for the Dart output type - {} ignored".format(encoding))

        filtered_calls, filtered_snps, filtered_samples = dataset.get_filtered_calls()

        if len(filtered_calls) == 0:
            log.info("All data has been silenced - nothing to output")
            return

        # Get calls as matrix - swap SNP/sample axes & trim back to only first allele's call (much quicker processing)
        numpy_matrix = numpy.asarray([filtered_calls[snp.allele_id] for snp in filtered_snps])
        numpy_matrix = numpy.dstack(numpy_matrix)  # [SNPs][samples][calls] -> [samples][calls][SNPs]
        numpy_matrix = numpy.dstack(numpy_matrix)  # [SNPs][calls][samples] -> [calls][SNPs][samples]
        # numpy_matrix = numpy_matrix  # Only get first allele calls (if "-" -> missing)

        del filtered_calls

        with open(file_path, "w") as file_out:
            writer = csv.writer(file_out, lineterminator='\n')

            pops_row = ["*"] * len(filtered_snps[0].all_headers) + [sample_def.population for sample_def in filtered_samples]
            samples_row = list(filtered_snps[0].all_headers.keys()) + [sample_def.id for sample_def in filtered_samples]

            writer.writerow(pops_row)
            writer.writerow(samples_row)

            for snp_idx, snp_def in enumerate(filtered_snps):
                snp_headers = list(snp_def.all_headers.values())

                snp_headers[list(snp_def.all_headers.keys()).index("AlleleID")] = snp_def.allele_id
                snp_headers[list(snp_def.all_headers.keys()).index("CloneID")] = snp_def.clone_id

                writer.writerow(snp_headers + numpy_matrix[0][snp_idx].tolist())

                # Replace the first rows sequence (ref seq.) with the second rows sequence
                if "AlleleSequence" in snp_def.all_headers:
                    snp_headers[list(snp_def.all_headers.keys()).index("AlleleSequence")] = snp_def.sequence

                writer.writerow(snp_headers + numpy_matrix[1][snp_idx].tolist())

            file_out.flush()

        # Write the read counts file
        file_path = os.path.join(folder, filter_name + "_" + self.get_name() + "_read_counts.csv")

        log.info("Outputting DArT Read Counts")

        filtered_counts, filtered_snps, filtered_samples = dataset.get_filtered_counts()
        numpy_matrix = numpy.asarray([filtered_counts[snp.allele_id] for snp in filtered_snps])
        numpy_matrix = numpy.dstack(numpy_matrix)  # [SNPs][samples][calls] -> [samples][calls][SNPs]
        numpy_matrix = numpy.dstack(numpy_matrix)  # [SNPs][calls][samples] -> [calls][SNPs][samples]

        del filtered_counts

        with open(file_path, "w") as file_out:
            writer = csv.writer(file_out, lineterminator='\n')

            a_header = filtered_snps[0].counts_headers if hasattr(filtered_snps[0], "counts_headers") else filtered_snps[0].all_headers

            # Add a row with * up to the first data column to mimic the dart format & so this output can be used as an input
            writer.writerow(["*"] * len(a_header) + [""] * len(filtered_samples))

            samples_row = list(a_header.keys()) + [sample_def.id for sample_def in filtered_samples]
            writer.writerow(samples_row)

            for snp_idx, snp_def in enumerate(filtered_snps):
                headers = snp_def.counts_headers if hasattr(snp_def, "counts_headers")  else snp_def.all_headers

                snp_headers = list(headers.values())

                snp_headers[list(headers.keys()).index("AlleleID")] = snp_def.allele_id
                snp_headers[list(headers.keys()).index("CloneID")] = snp_def.clone_id

                writer.writerow(snp_headers + numpy_matrix[0][snp_idx].tolist())

                # Replace the first rows sequence (ref seq.) with the second rows sequence
                if "AlleleSequence" in headers:
                    snp_headers[list(headers.keys()).index("AlleleSequence")] = snp_def.sequence

                writer.writerow(snp_headers + numpy_matrix[1][snp_idx].tolist())

            file_out.flush()

CSVOutput()
