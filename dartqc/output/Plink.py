import csv
import os

import logging
import numpy

from dartqc.Dataset import Dataset
from dartqc.PipelineOptions import Output

log = logging.getLogger(__file__)


class PlinkOutput(Output):
    def get_name(self) -> str:
        return "plink"

    def get_description(self) -> str:
        return "Output as PED and MAP files"

    def write(self, filter_name: str, folder: str, encoding:str, dataset: Dataset, unknown_args: [], **kwargs) -> None:
        file_path = os.path.join(folder, filter_name + "_" + self.get_name())

        ped_file = file_path + '.ped'
        map_file = file_path + '.map'

        filtered_calls, filtered_snps, filtered_samples = dataset.get_filtered_calls()

        log.info("Writing MAP file")

        # Write .map file - seems like its just a listing of allele ID's?
        map_data = [["0", snp_def, "0", "0"] for snp_def in filtered_snps]
        with open(map_file, 'w') as map_out:
            ped_writer = csv.writer(map_out, delimiter="\t")
            ped_writer.writerows(map_data)

        log.info("Writing PED file")

        # Get calls as matrix - swap SNP/sample axes & trim back to only first allele's call (much quicker processing)
        numpy_matrix = numpy.asarray([filtered_calls[snp.allele_id] for snp in filtered_snps])
        numpy_matrix = numpy.dstack(numpy_matrix)  # [SNPs][samples][calls] -> [samples][calls][SNPs]
        numpy_matrix = numpy.dstack(numpy_matrix)  # [SNPs][calls][samples] -> [calls][SNPs][samples]
        numpy_matrix = numpy_matrix  # Only get first allele calls (if "-" -> missing)

        del filtered_calls

        # Convert from (0,1) type tuples to requested encoding
        for snp_idx, snp_def in enumerate(filtered_snps):
            if encoding == "ACTG":
                # Find the two possible letters for this SNP
                snp_vals = snp_def.snp.split(":")[1].split(">")

                # Replace all 0's & 1's with ACTG's - missing replaced with 0
                numpy.put(numpy_matrix[0][snp_idx], numpy.where(numpy_matrix[0][snp_idx] == "0"), snp_vals[1])
                numpy.put(numpy_matrix[0][snp_idx], numpy.where(numpy_matrix[0][snp_idx] == "1"), snp_vals[0])
                numpy.put(numpy_matrix[1][snp_idx], numpy.where(numpy_matrix[1][snp_idx] == "0"), snp_vals[0])
                numpy.put(numpy_matrix[1][snp_idx], numpy.where(numpy_matrix[1][snp_idx] == "1"), snp_vals[1])

                numpy.put(numpy_matrix[0][snp_idx], numpy.where(numpy_matrix[0][snp_idx] == "-"), "0")
                numpy.put(numpy_matrix[1][snp_idx], numpy.where(numpy_matrix[1][snp_idx] == "-"), "0")
            elif encoding == "012":
                # Find indexes for 0's and 1's in first allele.
                major = numpy.where(numpy_matrix[0][snp_idx] == "1")
                minor = numpy.where(numpy_matrix[1][snp_idx] == "1")

                # Identify het's as where the second allele has a 1 in the same location as the first allele
                het = numpy.intersect1d(numpy.where(numpy_matrix[1][snp_idx] == "1"), major)

                # Remove het's from major and minor homozygous's
                major = numpy.setdiff1d(major, het, True)
                minor = numpy.setdiff1d(minor, het, True)

                # Replace the first allele values (-,0,1) with new encoding (-,0,1,2)
                numpy.put(numpy_matrix[0][snp_idx], het, "0")
                numpy.put(numpy_matrix[0][snp_idx], minor, "1")
                numpy.put(numpy_matrix[0][snp_idx], major, "2")

            elif encoding == "AB":
                snp_vals = ["A", "B"]

                # Replace all 0's & 1's with A's and B's - replace missing with 0
                numpy.put(numpy_matrix[0][snp_idx], numpy.where(numpy_matrix[0][snp_idx] == "0"), snp_vals[1])
                numpy.put(numpy_matrix[0][snp_idx], numpy.where(numpy_matrix[0][snp_idx] == "1"), snp_vals[0])
                numpy.put(numpy_matrix[1][snp_idx], numpy.where(numpy_matrix[1][snp_idx] == "0"), snp_vals[0])
                numpy.put(numpy_matrix[1][snp_idx], numpy.where(numpy_matrix[1][snp_idx] == "1"), snp_vals[1])

                numpy.put(numpy_matrix[0][snp_idx], numpy.where(numpy_matrix[0][snp_idx] == "-"), "0")
                numpy.put(numpy_matrix[1][snp_idx], numpy.where(numpy_matrix[1][snp_idx] == "-"), "0")

            elif encoding == "IUPAC":
                # Find the two possible letters for this SNP
                snp_vals = snp_def.snp.split(":")[1].split(">")

                both_snps = snp_vals[0] + snp_vals[1]

                het_val = "R" if both_snps == "AG" or both_snps == "GA" \
                    else "Y" if both_snps == "CT" or both_snps == "TC" \
                    else "S" if both_snps == "GC" or both_snps == "CG" \
                    else "W" if both_snps == "AT" or both_snps == "TA" \
                    else "K" if both_snps == "GT" or both_snps == "TG" \
                    else "M" if both_snps == "AC" or both_snps == "CA" \
                    else "???"

                # Find indexes for 0's and 1's in first allele.
                major = numpy.where(numpy_matrix[0][snp_idx] == "1")
                minor = numpy.where(numpy_matrix[1][snp_idx] == "1")

                # Identify het's as where the second allele has a 1 in the same location as the first allele
                het = numpy.intersect1d(numpy.where(numpy_matrix[1][snp_idx] == "1"), major)

                # Remove het's from major and minor homo's
                major = numpy.setdiff1d(major, het, True)
                minor = numpy.setdiff1d(minor, het, True)

                # Replace the first allele values (-,0,1) with new encoding (-,0,1,2)
                numpy.put(numpy_matrix[0][snp_idx], het, het_val)
                numpy.put(numpy_matrix[0][snp_idx], minor, snp_vals[1])
                numpy.put(numpy_matrix[0][snp_idx], major, snp_vals[0])

        if encoding == "012" or encoding == "IUPAC":
            # All data is only in the first allele - so grab it now and drop the second row
            numpy_matrix = numpy_matrix[0]
        else:
            # Reshape the matrix so that there are 2 times the SNPs ([calls][SNPs][Samples] -> [2*SNPs][samples])
            # (eg. basically the double row format dart has)
            numpy_matrix = numpy.reshape(numpy_matrix, (1, numpy_matrix.shape[1] * 2, numpy_matrix.shape[2]), order='F')[0]

        # Swap the axis directions ([SNPs][samples] -> [samples][SNPs])
        # Eg. For each sample there is 1 col per value (012 encoding only has 1 value but ACTG/AB has 2 per call)
        numpy_matrix = numpy.dstack(numpy_matrix)[0]  # [val][SNPs][samples] -> [samples][SNPs][val]

        with open(ped_file, 'w') as ped_out:
            for idx, sample_def in enumerate(filtered_samples):
                sample_details = [sample_def.population.replace(" ", "_"), sample_def.id.replace(" ", "_"), "0", "0", "0", "-9"]

                ped_out.write("\t".join(sample_details + numpy_matrix[idx].tolist()) + "\r\n")
                ped_out.flush()

                # MAP Formatting


PlinkOutput()
