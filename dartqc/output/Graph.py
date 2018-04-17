import csv
import os

import logging
import copy

import numpy

from dartqc.Dataset import Dataset
from dartqc.PipelineOptions import Output
from dartqc.Graphs import GraphTypes
from dartqc.filters.MinorAlleleFreq import MAFFilter

log = logging.getLogger(__file__)


class GraphOutput(Output):
    GRAPH_IMG_TYPE = ".png"

    def get_name(self) -> str:
        return "graph"

    def get_description(self) -> str:
        return "Generate a set of graphs as a summary of the dataset filtering changes"

    def write(self, filter_name: str, folder: str, encoding: str, dataset: Dataset, unknown_args: [], **kwargs) -> None:
        log.info("Creating Graphs")

        filtered_data, filtered_snps, filtered_samples = dataset.get_filtered_calls()

        maf = MAFFilter.calculate_maf(dataset, False)["pop"]

        filtered_read_counts = copy.deepcopy(dataset.read_counts)
        for allele_id in dataset.filtered.snps:
            del filtered_read_counts[allele_id]

        GraphTypes.read_counts_per_individ([filtered_read_counts], colors=["b"]).to_file(os.path.join(folder, filter_name +"_IndividReadCounts" + self.GRAPH_IMG_TYPE))
        GraphTypes.avg_reads_per_snp([filtered_read_counts], colors=["b"]).to_file(os.path.join(folder,filter_name + "_AvgReadsPerSNP" + self.GRAPH_IMG_TYPE))
        GraphTypes.avg_rep_across_snp([filtered_snps], colors=["b"]).to_file(os.path.join(folder,filter_name + "_AvgRepAcrossSNP" + self.GRAPH_IMG_TYPE))
        GraphTypes.het_across_snp([filtered_data], colors=["b"]).to_file(os.path.join(folder,filter_name + "_HetAcrossSNP" + self.GRAPH_IMG_TYPE))
        GraphTypes.call_rates_across_snp([filtered_data], colors=["b"]).to_file(os.path.join(folder,filter_name + "_CallRatesAcrossSNP" + self.GRAPH_IMG_TYPE))
        GraphTypes.call_rate_across_individ([filtered_data], colors=["b"]).to_file(os.path.join(folder,filter_name + "_CallRatesAcrossIndivid" + self.GRAPH_IMG_TYPE))
        GraphTypes.maf_across_snp([maf], colors=["b"]).to_file(os.path.join(folder,filter_name + "_MAFAcrossSNP" + self.GRAPH_IMG_TYPE))
        GraphTypes.maf_to_read_count([maf], [filtered_read_counts], colors=["b"]).to_file(os.path.join(folder,filter_name + "_MAFToReadCount" + self.GRAPH_IMG_TYPE))
        GraphTypes.call_rate_to_maf([filtered_data], [maf], colors=["b"]).to_file(os.path.join(folder, filter_name + "_CallRateToMAF" + self.GRAPH_IMG_TYPE))

GraphOutput()
