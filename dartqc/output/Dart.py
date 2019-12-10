import csv
import os

import logging

import numpy
import copy
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

        if encoding is not None and encoding != "11" and encoding != "undefined":
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

        # # --------------Update to un-collapse replicates------------
        # rep_counts = copy.copy(dataset.replicate_counts)
        # Get a copy of all replicate counts (with only unfiltered SNPs)
        rep_counts = numpy.asarray([copy.copy(dataset.replicate_counts[snp]) for snp in dataset.replicate_counts if snp not in dataset.filtered.snps])

        # Delete any replicates for samples that have been filtered out
        del_rep_idxs = [idx for idx, rep in enumerate(dataset.replicates) if rep in dataset.filtered.samples]
        rep_counts = numpy.delete(rep_counts, del_rep_idxs, 1)

        rep_samples = [rep for rep in dataset.replicates if rep not in dataset.filtered.samples]

        # Get the list of sample IDs as per the count of replicates
        rep_sample_headers = []
        for rep_idx, rep_values in enumerate(rep_counts[0]):
            rep_id = rep_samples[rep_idx]

            for counts in rep_values:
                rep_sample_headers.append(rep_id)


        rep_sample_idxs = [idx for idx, sample_def in enumerate(filtered_samples) if sample_def.id in rep_samples]

        # delete existing sample defs
        filtered_samples = [sample_def for sample_def in filtered_samples if sample_def.id not in rep_samples]

        # Delete read counts from data (then add the seperated replicate counts at the end)
        for allele_id, counts in filtered_counts.items():
            filtered_counts[allele_id] = numpy.delete(counts, rep_sample_idxs, 0)

        # Now add them back in with all replicates
        numpy_matrix = numpy.asarray([filtered_counts[snp.allele_id] for snp in filtered_snps])

        # Collapse replicates array (so it matches the read_counts array)
        all_reps_array = []
        for snp_idx, snp in enumerate(filtered_snps):
            reps_array = numpy.asarray([counts for rep_counts in rep_counts[snp_idx] for counts in rep_counts])

            all_reps_array.append(reps_array)

        all_reps_array = numpy.asarray(all_reps_array)

        # Append the replicate array read counts to the end of the existing read counts.
        numpy_matrix = numpy.concatenate((numpy_matrix, all_reps_array), axis=1)

        # Re-ordder the matrix to seperate out the two calls (quicker/easier to write to CSV)
        numpy_matrix = numpy.dstack(numpy_matrix)  # [SNPs][samples][calls] -> [samples][calls][SNPs]
        numpy_matrix = numpy.dstack(numpy_matrix)  # [SNPs][calls][samples] -> [calls][SNPs][samples]

        del filtered_counts

        with open(file_path, "w") as file_out:
            writer = csv.writer(file_out, lineterminator='\n')

            a_header = filtered_snps[0].counts_headers if hasattr(filtered_snps[0], "counts_headers") else filtered_snps[0].all_headers

            # Add a row with * up to the first data column to mimic the dart format & so this output can be used as an input
            writer.writerow(["*"] * len(a_header) + [""] * (len(filtered_samples) + len(rep_sample_headers)))

            samples_row = list(a_header.keys()) + [sample_def.id for sample_def in filtered_samples] + rep_sample_headers
            writer.writerow(samples_row)

            for snp_idx, snp_def in enumerate(filtered_snps):
                headers = snp_def.counts_headers if hasattr(snp_def, "counts_headers") else snp_def.all_headers

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
