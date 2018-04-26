import csv
import os

import logging
import copy

import numpy

from dartqc.Dataset import Dataset
from dartqc.PipelineOptions import Output
from dartqc.Graphs import GraphTypes
from dartqc.filters.MinorAlleleFreq import MAFFilter
from dartqc.filters.RepAvg import RepAvgFilter
from dartqc.filters.HardyWeinberg import HWEFilter

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

        hwe_disequilibrium = HWEFilter.calculate_hwe_disequilibrium(dataset, False)["pop"]

        filtered_read_counts = dataset.get_filtered_counts()[0]

        rep_avgs = RepAvgFilter.calc_rep_avg(dataset)

        GraphTypes.read_counts_per_individ([filtered_read_counts], colors=["b"]).to_file(os.path.join(folder, filter_name +"_IndividReadCounts" + self.GRAPH_IMG_TYPE))
        GraphTypes.avg_reads_per_snp([filtered_read_counts], colors=["b"]).to_file(os.path.join(folder,filter_name + "_AvgReadsPerSNP" + self.GRAPH_IMG_TYPE))
        GraphTypes.avg_rep_across_snp([rep_avgs], colors=["b"]).to_file(os.path.join(folder,filter_name + "_AvgRepAcrossSNP" + self.GRAPH_IMG_TYPE))
        GraphTypes.het_across_snp([filtered_data], colors=["b"]).to_file(os.path.join(folder,filter_name + "_HetAcrossSNP" + self.GRAPH_IMG_TYPE))
        GraphTypes.call_rates_across_snp([filtered_data], colors=["b"]).to_file(os.path.join(folder,filter_name + "_CallRatesAcrossSNP" + self.GRAPH_IMG_TYPE))
        GraphTypes.call_rate_across_individ([filtered_data], colors=["b"]).to_file(os.path.join(folder,filter_name + "_CallRatesAcrossIndivid" + self.GRAPH_IMG_TYPE))
        GraphTypes.maf_across_snp([maf], colors=["b"]).to_file(os.path.join(folder,filter_name + "_MAFAcrossSNP" + self.GRAPH_IMG_TYPE))
        GraphTypes.maf_to_read_count([maf], [filtered_read_counts], colors=["b"]).to_file(os.path.join(folder,filter_name + "_MAFToReadCount" + self.GRAPH_IMG_TYPE))
        GraphTypes.call_rate_to_maf([filtered_data], [maf], colors=["b"]).to_file(os.path.join(folder, filter_name + "_CallRateToMAF" + self.GRAPH_IMG_TYPE))
        GraphTypes.maf_to_hwe([maf], [hwe_disequilibrium], [filtered_read_counts], colors=["b"]).to_file(os.path.join(folder, filter_name + "_MAFToHWE" + self.GRAPH_IMG_TYPE))
        # GraphTypes.read_depth_outliers([filtered_read_counts], colors=["b"]).to_file(os.path.join(folder, filter_name + "_ReadDepthOutliers" + self.GRAPH_IMG_TYPE), is_3d=False)

GraphOutput()
