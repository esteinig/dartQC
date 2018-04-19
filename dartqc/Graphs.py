import random
import time
from math import ceil

import logging
import numpy
import os
from matplotlib import pyplot
from matplotlib import patches
from matplotlib import colors

from dartqc.Dataset import SNPDef, Dataset

GRAPH_IMG_TYPE = ".jpg"

log = logging.getLogger(__file__)


class Graph:
    SCATTERPLOT = 0
    LINE_GRAPH = 1
    BAR_GRAPH = 2

    def __init__(self, title: str, graph_type: int, x_label: str = None, y_label: str = None, x_ticks: [] = None, x_tick_labels: [] = None,
                 y_ticks: [] = None, y_tick_labels: [] = None, width: int = 15, height: int = 9):
        self.title = title

        self.x_label = x_label
        self.x_ticks = x_ticks
        self.x_tick_labels = x_tick_labels

        self.y_label = y_label
        self.y_tick_labels = y_tick_labels
        self.y_ticks = y_ticks

        self.graph_type = graph_type

        self.width = width
        self.height = height

        self.sets = []

    def add_set(self, x_data: [], y_data: [], legend: str = None, color=None):
        if color is None:
            color = random.choice(list(colors.get_named_colors_mapping().keys()))

        if legend is None:
            legend = "Set {}".format(len(self.sets) + 1)

        if x_data is None and self.graph_type == self.BAR_GRAPH:
            x_data = [i for i in range(len(y_data))]

        self.sets.append(GraphSet(x_data, y_data, legend, color))

    def to_file(self, file_path):
        fig, ax = pyplot.subplots()
        fig.set_size_inches(self.width, self.height)

        legend_handles = []

        for set_idx, set in enumerate(self.sets):
            if self.graph_type == self.SCATTERPLOT:
                # for idx, y_val in enumerate(set.y_data):
                ax.scatter(set.x_data, set.y_data, .5, color=set.color)
            elif self.graph_type == self.BAR_GRAPH:
                full_width = .8
                width = full_width / len(self.sets)
                x = numpy.asarray(set.x_data) + ((set_idx * width) + (width / 2) - full_width / 2)

                # for idx, y_val in enumerate(set.y_data):
                ax.bar(x, set.y_data, width, color=set.color)
            else:
                raise Exception("Invalid graph type: {}".format(self.graph_type))

            patch = patches.Patch(color=set.color, label=set.legend)
            legend_handles.append(patch)

        pyplot.legend(handles=legend_handles)

        # Make the labels look nicer.
        ax.set_title(self.title, fontdict={'weight': 'bold', 'size': 'large'})

        # X labels/ticks
        if self.x_label is not None:
            ax.set_xlabel(self.x_label, fontdict={'weight': 'bold'})

        if self.x_ticks is not None:
            ax.set_xticks(self.x_ticks)

        if self.x_tick_labels is not None:
            ax.set_xticklabels(self.x_tick_labels, fontdict={'style': 'italic'}, rotation=45)

        # Make the x-axis labels align correctly
        for tick in ax.xaxis.get_majorticklabels():
            tick.set_horizontalalignment("right")

        # Y labels/ticks
        if self.y_label is not None:
            ax.set_ylabel(self.y_label, fontdict={'weight': 'bold'})

        if self.y_ticks is not None:
            ax.set_yticks(self.y_ticks)

        if self.y_tick_labels is not None:
            ax.set_yticklabels(self.y_tick_labels, fontdict={'style': 'italic'})

        # Save the graph to file
        # log.info("Generating graph: " + file_path)
        pyplot.savefig(file_path, bbox_inches='tight', transparent=True)

        pyplot.close(fig)


class GraphSet:
    def __init__(self, x_data: [], y_data: [], legend: str = None, color: str = None):
        self.x_data = x_data
        self.y_data = y_data
        self.legend = legend
        self.color = color


class GraphTypes:
    @staticmethod
    def _find_categories(y_data: [[]], num_categories: int = 10):
        biggest = 0
        for set_data in y_data:
            for val in set_data:
                if val > biggest:
                    biggest = val

        mag = 2.0
        while pow(10, mag) < biggest:
            mag += 1

        mag = pow(10, mag)
        if biggest / mag < .5:
            mag /= 10

        biggest = ceil(biggest / mag) * mag

        label_divisor = 1000000 if (mag / num_categories) >= 1000000 else 1000 if (mag / num_categories) >= 1000 else 1
        label_unit = "M" if label_divisor == 1000000 else "k" if label_divisor == 1000 else ""

        categories = [round((x / (num_categories + 1)) * biggest) for x in range(1, num_categories + 1)]
        x_tick_labels = [("Less than " if idx == 0
                          else str(round(categories[idx - 1] / label_divisor)) + " to ") + str(round(num / label_divisor)) + label_unit if idx < num_categories - 1
                         else "Greater than " + str(round(categories[len(categories) - 2] / label_divisor)) + label_unit for idx, num in enumerate(categories)]

        return categories, x_tick_labels, [i for i in range(len(categories))]

    @staticmethod
    def _to_category_counts(values: [], categories: []):
        results = [0 for i in range(len(categories))]

        for val in values:
            found = False
            for idx, num in enumerate(categories):
                if val < num:
                    results[idx] += 1
                    found = True
                    break

            if not found:
                results[len(results) - 1] += 1

        return results

    # (Graph 1) Distribution of total read count per individual
    @staticmethod
    def read_counts_per_individ(read_counts: [[]], colors: [] = None, legends: [] = None):
        start = time.time()

        graph = Graph(
            title="Distribution of total read counts per individual",
            graph_type=Graph.BAR_GRAPH,
            y_label="Number of individuals",
            x_label="Total read count per individual"
        )

        all_indiv_read_cnts = []
        for set_counts in read_counts:
            if len(set_counts) == 0:
                continue

            # swap the matrix [snp][sample][allele] -> [sample][snp][allele]
            individ_counts = numpy.asarray([counts for allele_id, counts in set_counts.items()])
            individ_counts = numpy.stack(individ_counts, axis=1)  # [SNPs][samples][calls] -> [SNPs][calls][samples]

            # Add all counts per individual allele pair together
            individ_counts = numpy.sum(individ_counts, axis=(1, 2))

            all_indiv_read_cnts.append(individ_counts)

        categories, graph.x_tick_labels, graph.x_ticks = GraphTypes._find_categories(all_indiv_read_cnts, 10)

        # Count each value based on which categoruy it falls into
        for idx, individ_read_counts in enumerate(all_indiv_read_cnts):
            set_y_data = GraphTypes._to_category_counts(individ_read_counts, categories)

            graph.add_set(x_data=None,
                          y_data=set_y_data,
                          legend=legends[idx] if legends is not None and len(legends) > idx else None,
                          color=colors[idx] if colors is not None and len(colors) > idx else None)

        log.info(graph.title + ": " + str(round((time.time() - start), 2)) + "s")

        return graph

    # (Graph 2) Distribution of average reads per SNP
    @staticmethod
    def avg_reads_per_snp(read_counts: [[]], colors: [] = None, legends: [] = None):
        start = time.time()

        #  Distribution of total read counts per individual
        graph = Graph(
            title="Distribution of average reads per SNP",
            graph_type=Graph.BAR_GRAPH,
            y_label="Number of SNPs",
            x_label="Average number of reads per SNP",
            x_tick_labels=["<5", "5 to 10", "10 to 15", "15 to 20", "20 to 25", "25 to 30", "30 to 35", ">35"]
        )

        categories = [5, 10, 15, 20, 25, 30, 35, 1000]
        graph.x_ticks = [i for i in range(len(categories))]

        for idx, set_counts in enumerate(read_counts):
            if len(set_counts) == 0:
                continue

            set_counts = numpy.asarray([counts for allele_id, counts in set_counts.items()])
            reads_per_snp = numpy.sum(set_counts, axis=(1, 2))  # / 2.0
            num_samples = numpy.shape(set_counts)[1]
            avg_reads_per_snp = (reads_per_snp / num_samples).tolist()

            set_y_data = GraphTypes._to_category_counts(avg_reads_per_snp, categories)
            graph.add_set(x_data=None,
                          y_data=set_y_data,
                          legend=legends[idx] if legends is not None and len(legends) > idx else None,
                          color=colors[idx] if colors is not None and len(colors) > idx else None)

        log.info(graph.title + ": " + str(round((time.time() - start), 2)) + "s")

        return graph

    # (Graph 3) Distribution of call rates across SNPs
    @staticmethod
    def call_rates_across_snp(calls, colors: [] = None, legends: [] = None):
        start = time.time()

        #  Distribution of total read counts per individual
        graph = Graph(
            title="Distribution of call rates across SNPs",
            graph_type=Graph.BAR_GRAPH,
            y_label="Number of SNPs",
            x_label="Call rate per SNP",
            x_tick_labels=["<50%", "50 to 55%", "55 to 60%", "60 to 65%", "65 to 70%", "70 to 75%", "75 to 80%", "80 to 85%",
                           "85 to 90%", "90 to 95%", "95 to 100%"]
        )

        categories = [.5, .55, .60, .65, .70, .75, .80, .85, .90, .95, 100]
        graph.x_ticks = [i for i in range(len(categories))]

        for idx, graph_data in enumerate(calls):
            if len(graph_data) == 0:
                continue

            # blanks = [(calls == "-").sum() for calls in callsArray]
            # positive_calls = [numpy.asarray([(call[0] == "1" and call[1] == "0") or (call[1] == "1" and call[0] == "0") for call in call_list]).sum() for call_list in graph_data.values()]
            # tot_calls = [numpy.asarray([call[0] != "-" for call in call_list]).sum() for call_list in graph_data.values()]
            #
            # call_rate_per_snp = [0 if tot == 0 else pos / tot for pos, tot in zip(positive_calls, tot_calls)]

            numpy_matrix = numpy.asarray([calls for calls in graph_data.values()])
            numpy_matrix = numpy.dstack(numpy_matrix)  # [SNPs][samples][calls] -> [samples][calls][SNPs]
            numpy_matrix = numpy.dstack(numpy_matrix)  # [samples][calls][SNPs] -> [calls][SNPs][samples]
            numpy_matrix = numpy_matrix[0]  # Only get first allele calls (if "-" -> missing)

            called_counts = numpy.asarray([len(numpy.where(calls != "-")[0]) for calls in numpy_matrix])
            call_rate_per_snp = called_counts / numpy.shape(numpy_matrix)[1]

            set_y_data = GraphTypes._to_category_counts(call_rate_per_snp, categories)
            graph.add_set(x_data=None,
                          y_data=set_y_data,
                          legend=legends[idx] if legends is not None and len(legends) > idx else None,
                          color=colors[idx] if colors is not None and len(colors) > idx else None)

        log.info(graph.title + ": " + str(round((time.time() - start), 2)) + "s")
        return graph

    # (Graph 4) Distribution of call rates across individuals
    @staticmethod
    def call_rate_across_individ(calls, colors: [] = None, legends: [] = None):
        start = time.time()

        graph = Graph(
            title="Distribution of call rates across individuals",
            graph_type=Graph.BAR_GRAPH,
            y_label="Number of individuals",
            x_label="Call rate per individual",
            x_tick_labels=["<50%", "50 to 55%", "55 to 60%", "60 to 65%", "65 to 70%", "70 to 75%", "75 to 80%", "80 to 85%",
                           "85 to 90%", "90 to 95%", "95 to 100%"]
        )
        categories = [.50, .55, .60, .65, .70, .75, .80, .85, .90, .95, .100]
        graph.x_ticks = [i for i in range(len(categories))]

        for idx, graph_data in enumerate(calls):
            if len(graph_data) == 0:
                continue

            # blanks = [(calls == "-").sum() for calls in callsArray]
            # positive_calls = [((calls == "1") | (calls == "2")).sum() for calls in calls_array]
            # tot_calls = [(calls != "-").sum() for calls in calls_array]
            # positive_calls = [numpy.asarray([call == "1" or call == "2" for call in call_list]).sum() for call_list in graph_data]
            # tot_calls = [numpy.asarray([call != "-" for call in call_list]).sum() for call_list in graph_data]
            #
            # call_rate_per_snp = [0 if tot == 0 else pos / tot for pos, tot in zip(positive_calls, tot_calls)]

            # Get calls as matrix - swap SNP/sample axes & trim back to only first allele's call (much quicker processing)
            numpy_matrix = numpy.asarray([calls for calls in graph_data.values()])
            numpy_matrix = numpy.stack(numpy_matrix, axis=2)  # [SNPs][samples][calls] -> [SNPs][calls][samples]
            numpy_matrix = numpy.stack(numpy_matrix, axis=1)  # [SNPs][calls][samples] -> [calls][samples][SNPs]
            numpy_matrix = numpy_matrix[0]  # Only get first allele calls (if "-" -> missing)

            called_counts = numpy.asarray([len(numpy.where(calls != "-")[0]) for calls in numpy_matrix])
            call_rate_per_sample = called_counts / numpy.shape(numpy_matrix)[1]

            set_y_data = GraphTypes._to_category_counts(call_rate_per_sample, categories)
            graph.add_set(x_data=None,
                          y_data=set_y_data,
                          legend=legends[idx] if legends is not None and len(legends) > idx else None,
                          color=colors[idx] if colors is not None and len(colors) > idx else None)

        log.info(graph.title + ": " + str(round((time.time() - start), 2)) + "s")

        return graph

    # (Graph 5) Distribution of MAF across SNPs
    @staticmethod
    def maf_across_snp(snp_maf, colors: [] = None, legends: [] = None):
        start = time.time()

        graph = Graph(
            title="Distribution of MAF across SNPs",
            graph_type=Graph.BAR_GRAPH,
            y_label="Number of SNPs",
            x_label="Minor Allele Frequency (MAF) of SNPs",
            x_tick_labels=["less than 0.01", "0.01 to 0.02", "0.02 to 0.03", "0.03 to 0.04", "0.04 to 0.05", "0.05 to 0.06",
                           "0.06 to 0.07", "0.07 to 0.08", "0.08 to 0.09", "0.09 to 0.1", "0.1 to 0.15", "0.15 to 0.2",
                           "0.2 to 0.3", "0.3 to 0.4", "0.4 to 0.5", "greater than 0.5"]
        )
        categories = [.01, .02, .03, .04, .05, .06, .07, .08, .09, .1, .15, .2, .3, .4, .5, 1]
        graph.x_ticks = [i for i in range(len(categories))]

        for idx, set_maf in enumerate(snp_maf):
            if len(set_maf) == 0:
                continue

            set_y_data = GraphTypes._to_category_counts(list(set_maf.values()), categories)
            graph.add_set(x_data=None,
                          y_data=set_y_data,
                          legend=legends[idx] if legends is not None and len(legends) > idx else None,
                          color=colors[idx] if colors is not None and len(colors) > idx else None)

        log.info(graph.title + ": " + str(round((time.time() - start), 2)) + "s")

        return graph

    # (Graph 6) Distribution of average repeatability (repAvg) across SNPs
    @staticmethod
    def avg_rep_across_snp(snp_defs: [[SNPDef]], colors: [] = None, legends: [] = None):
        start = time.time()

        graph = Graph(
            title="Distribution of average repeatability (RepAvg) across SNPs",
            graph_type=Graph.BAR_GRAPH,
            y_label="Number of SNPs",
            x_label="Average repeatability of SNPs",
            x_tick_labels=["<.9", "<.92", "<.94", "<.96", "<.98", "<.99", ">.99"]
        )
        categories = [.9, .92, .94, .96, .98, .99, 1]
        graph.x_ticks = [i for i in range(len(categories))]

        for idx, snp_defs in enumerate(snp_defs):
            if len(snp_defs) == 0:
                continue

            # Get the repeatability averages
            # TODO: Own calc of rep avg.
            rep_avgs = [float(snp_def.rep_average) for snp_def in snp_defs]

            set_y_data = GraphTypes._to_category_counts(rep_avgs, categories)
            graph.add_set(x_data=None,
                          y_data=set_y_data,
                          legend=legends[idx] if legends is not None and len(legends) > idx else None,
                          color=colors[idx] if colors is not None and len(colors) > idx else None)

        log.info(graph.title + ": " + str(round((time.time() - start), 2)) + "s")

        return graph

    # (Graph 7) Distribution of Heterozygosity across SNPs
    @staticmethod
    def het_across_snp(calls: {str:numpy.array}, colors: [] = None, legends: [] = None):
        start = time.time()

        graph = Graph(
            title="Distribution of Heterozygosity across SNPs",
            graph_type=Graph.BAR_GRAPH,
            y_label="Number of SNPs",
            x_label="Heterozygosity per SNP",
            x_tick_labels=["0 to 10%", "10 to 20%", "20 to 30%", "30 to 40%", "40 to 50%", "50 to 60%", "60 to 70%",
                           "70 to 80%", "80 to 90%", "90 to 100%"]
        )
        categories = [val / 100 for val in range(10, 101, 10)]
        graph.x_ticks = [i for i in range(len(categories))]

        for idx, set_calls in enumerate(calls):
            if len(set_calls) == 0:
                continue

            freq_hetz = []
            for allele_id, calls in set_calls.items():
                swapped_calls = numpy.dstack(calls)[0]  # [samples][calls][SNPs] -> [calls][SNPs][samples]

                num_calls = len(numpy.where(swapped_calls[0] != "-")[0])
                allele_1_idxs = numpy.where(swapped_calls[0] == "1")[0]
                allele_2_idxs = numpy.where(swapped_calls[1] == "1")[0]
                het_idxs = numpy.intersect1d(allele_1_idxs, allele_2_idxs)

                freq_hetz.append(len(het_idxs) / num_calls)

            # freq_hetz = [float(snp_def.all_headers["FreqHomSnp"]) for snp_def in set_calls]
            set_y_data = GraphTypes._to_category_counts(freq_hetz, categories)
            graph.add_set(x_data=None,
                          y_data=set_y_data,
                          legend=legends[idx] if legends is not None and len(legends) > idx else None,
                          color=colors[idx] if colors is not None and len(colors) > idx else None)

        log.info(graph.title + ": " + str(round((time.time() - start), 2)) + "s")

        return graph

    # (Graph 8) Relationship between MAF and Read Count
    @staticmethod
    def maf_to_read_count(snp_maf: [[]], read_counts: [numpy.array], colors: [] = None, legends: [] = None):
        start = time.time()

        #  Make sure that the read_counts and snp_maf SNP order/count match (eg. if MAF calc'ed on a filtered dataset, read counts will need SNP's removed)

        graph = Graph(
            title="Relationship between MAF and Read Count",
            graph_type=Graph.SCATTERPLOT,
            y_label="MAF by SNP",
            x_label="Read Count by SNP",
        )

        for (idx, maf_values) in enumerate(snp_maf):
            if len(maf_values) == 0:
                continue

            counts = numpy.asarray([counts for allele_id, counts in read_counts[idx].items()])
            reads_per_snp = (numpy.sum(counts, axis=(1, 2))).tolist() #/ 2.0

            # log.info("{} - {}".format(len(counts), len(reads_per_snp)))

            maf_vals = list(maf_values.values())
            reads_per_snp = reads_per_snp

            graph.add_set(x_data=reads_per_snp,
                          y_data=maf_vals,
                          legend=legends[idx] if legends is not None and len(legends) > idx else None,
                          color=colors[idx] if colors is not None and len(colors) > idx else None)

        log.info(graph.title + ": " + str(round((time.time() - start), 2)) + "s")

        return graph

    # (Graph 9) Relationship between Call rate and MAF
    @staticmethod
    def call_rate_to_maf(calls: [numpy.array], snp_maf, colors:[]=None, legends:[]=None):
        start = time.time()

        graph = Graph(
            title="Relationship between Call rate and MAF",
            graph_type=Graph.SCATTERPLOT,
            y_label="MAF by SNP",
            x_label="Call rate by SNP",
            x_tick_labels=["0%", "20%", "40%", "60%", "80%", "100%"],
            x_ticks=[0, .2, .4, .6, .8, 1]
        )

        for (idx, set_calls) in enumerate(calls):
            if len(set_calls) == 0:
                continue

            # positive_calls = [numpy.asarray([call == "1" or call == "2" for call in call_list]).sum() for call_list in set_calls]
            # tot_calls = [numpy.asarray([call != "-" for call in call_list]).sum() for call_list in set_calls]
            #
            # call_rate_per_snp = [0 if tot == 0 else pos / tot for pos, tot in zip(positive_calls, tot_calls)]

            numpy_matrix = numpy.asarray([calls for calls in set_calls.values()])
            numpy_matrix = numpy.dstack(numpy_matrix)  # [SNPs][samples][calls] -> [samples][calls][SNPs]
            numpy_matrix = numpy.dstack(numpy_matrix)  # [samples][calls][SNPs] -> [calls][SNPs][samples]
            numpy_matrix = numpy_matrix[0]  # Only get first allele calls (if "-" -> missing)

            called_counts = numpy.asarray([len(numpy.where(calls != "-")[0]) for calls in numpy_matrix])
            call_rate_per_snp = called_counts / numpy.shape(numpy_matrix)[1]

            # Find MAF for a SNP
            graph_y_data = list(snp_maf[idx].values())

            graph.add_set(x_data=call_rate_per_snp,
                          y_data=graph_y_data,
                          legend=legends[idx] if legends is not None and len(legends) > idx else None,
                          color=colors[idx] if colors is not None and len(colors) > idx else None)

        log.info(graph.title + ": " + str(round((time.time() - start), 2)) + "s")

        return graph
