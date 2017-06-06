import os
import csv
import numpy
import json

from dartqc.DartUtils import stamp


class DartWriter:

    def __init__(self, data, attributes):

        self.data = data
        self.attributes = attributes

        self.decoding_scheme = dict()

        self.set_encoding()

    def set_encoding(self, homozygous_major=("A", "A"), homozygous_minor=("B", "B"), heterozygous=("A", "B"),
                     missing=("0", "0")):

        self.decoding_scheme = {"-": missing, "0": heterozygous, "1": homozygous_minor, "2": homozygous_major}

    def write_plink(self, file_name, sep="\t", remove_space=False):

        snp_order = sorted(self.data.keys())

        stamp("Decoding calls...")
        snp_rows = [[self.decoding_scheme[snp] for snp in self.data[snp_id]["calls"]] for snp_id in snp_order]

        stamp("Transposing calls...")
        snps_by_sample = numpy.asarray(snp_rows).transpose(1, 0, 2)

        genotypes = [sample.flatten().tolist() for sample in snps_by_sample]

        names = self.attributes["sample_names"]
        pops = [self.attributes["pops"][sample] for sample in names]

        if remove_space:
            names = ["_".join(name.split()) for name in names]
            pops = ["_".join(pop.split()) for pop in pops]

        ped_file = os.path.join(self.attributes["out_path"], file_name + '.ped')
        map_file = os.path.join(self.attributes["out_path"], file_name + '.map')

        paternal = ["0"] * len(names)
        maternal = ["0"] * len(names)
        sex = ["0"] * len(names)
        phenotype = ["-9"] * len(names)

        plink = zip(pops, names, paternal, maternal, sex, phenotype, genotypes)

        stamp("Formatting calls...")

        ped_data = []
        for row in plink:
            new_row = list(row[:6])
            for geno in row[6]:
                new_row.append(geno)
            ped_data.append(new_row)

        stamp("Writing PLINK")
        stamp("PED file:", ped_file)
        stamp("MAP file:", map_file)

        with open(ped_file, 'w') as ped_out:
            ped_writer = csv.writer(ped_out, delimiter=sep)
            ped_writer.writerows(ped_data)

        # MAP Formatting

        map_data = [["0", snp_id, "0", "0"] for snp_id in snp_order]

        with open(map_file, 'w') as map_out:
            ped_writer = csv.writer(map_out, delimiter=sep)
            ped_writer.writerows(map_data)

    def write_json(self, file_name, data_indent=0, attribute_indent=4):

        data_file = os.path.abspath(os.path.join(self.attributes["out_path"], file_name + "_data.json"))
        attribute_file = os.path.abspath(os.path.join(self.attributes["out_path"], file_name + "_attr.json"))

        stamp("Writing data to JSON")
        stamp("Data file:", data_file)
        stamp("Attribute file:", attribute_file)

        with open(data_file, "w") as data_out:
            json.dump(self.data, data_out, indent=data_indent)

        with open(attribute_file, "w") as attr_out:
            json.dump(self.attributes, attr_out, indent=attribute_indent)
