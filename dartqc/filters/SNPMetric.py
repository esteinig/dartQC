import logging
import numpy
import re

from dartqc import PipelineOptions
from dartqc.Dataset import Dataset
from dartqc.FilterResult import FilterResult

log = logging.getLogger(__file__)


class SNPMetricFilter(PipelineOptions.Filter):
    def get_name(self) -> str:
        return "SNP Metric"

    def get_cmd_names(self):
        return ["--metric", "--rep_average"]

    def get_cmd_help(self) -> str:
        return "Filter based on a quality metric provided in the data (dataset.all_headers).  " \
               "Pattern: [<col_name><comparison><value>] such as [RepAvg>0.2,RepAvg>0.3].  " \
               "Comparisons includ <, > and ="

    def get_description(self) -> str:
        return "Filter based on a quality metric provided in the data (make sure the header matches exactly!)"

    def get_cmd_type(self):
        return lambda s: [(re.split(r'[<=>]', item)[0], re.findall(r'[<=>]', item)[0], float(re.split(r'[<=>]', item)[1]))
                          if len(item.strip()) > 0 else None for item in re.sub(r'(^[(\[]{2})|([)\]]{2}$)', r"\1", s).split(',')]

    def get_order(self) -> int:
        return 5

    def filter(self, dataset: Dataset, threshold, unknown_args: [], **kwargs) -> FilterResult:
        silenced = FilterResult()

        for snp_def in dataset.snps:
            if snp_def.allele_id in dataset.filtered.snps:
                continue

            if threshold[1] == ">" and float(snp_def.all_headers[threshold[0]]) < threshold[2]:
                silenced.silenced_snp(snp_def.allele_id)
            elif threshold[1] == "<" and float(
                    snp_def.all_headers[threshold[0]]) > threshold[2]:
                silenced.silenced_snp(snp_def.allele_id)
            elif threshold[1] == "=" and float(snp_def.all_headers[threshold[0]]) != threshold[2]:
                silenced.silenced_snp(snp_def.allele_id)
            elif threshold[1] != ">" and threshold[1] != "=" and threshold[1] != "<":
                log.error("Invalid comparison type: {} only < = or > are supported".format(threshold[1]))

        return silenced

SNPMetricFilter()
