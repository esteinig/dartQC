import random
import time

import numpy
from matplotlib import pyplot
from matplotlib import patches
from matplotlib import colors

from DartModules import MarkerModule


class Graphs:
    # Some graphs won't be changed based on filtering (eg. read counts, repAvg & freq Hetz)
    @staticmethod
    def create_static_plots(data, read_data, output_dir, project):
        print("Creating static plots")

        Graphs.read_counts_per_individ(read_data, output_dir + project + "_IndividReadCounts" + ".jpg", color='blue')
        Graphs.avg_reads_per_snp(read_data, output_dir + project + "_AvgReadsPerSNP" + ".jpg", color='blue')
        Graphs.avg_rep_across_snp(data, output_dir + project + "_AvgRepAcrossSNP" + ".jpg", color='blue')
        Graphs.het_across_snp(data, output_dir + project + "_HetAcrossSNP" + ".jpg", color='blue')

    @staticmethod
    def create_plots(data, read_data, attrs, name, output_dir, project, color=None, legend=None):
        print("Creating " + name + " plots")

        if not isinstance(data, list):
            data = [data]
        if not isinstance(attrs, list):
            attrs = [attrs]

        snpMAF = []
        for (index, data_item) in enumerate(data):
            mm = MarkerModule(data=data_item, attributes=attrs[index])
            snpMAF.append({k: v["maf"] for (k, v) in data_item.items()})

        Graphs.call_rates_across_snp(data=data, outfile=output_dir + project + "_" + name + "_CallRatesAcrossSNP" + ".jpg", color=color, legend=legend)
        Graphs.call_rate_across_individ(data=data, outfile=output_dir + project + "_" + name + "_CallRatesAcrossIndivid" + ".jpg", color=color, legend=legend)
        Graphs.maf_across_snp(snp_maf=snpMAF, outfile=output_dir + project + "_" + name + "_MAFAcrossSNP" + ".jpg", color=color, legend=legend)
        Graphs.maf_to_read_count(snp_maf=snpMAF, read_data=read_data, outfile=output_dir + project + "_" + name + "_MAFToReadCount" + ".jpg", color=color, legend=legend)
        Graphs.call_rate_to_maf(data=data, snp_maf=snpMAF, outfile=output_dir + project + "_" + name + "_CallRateToMAF" + ".jpg", color=color, legend=legend)


    @staticmethod
    def create_bar_graph(y, title, x_tick_labels, x_label, y_label, outfile, color=None, legend=None):
        # ax = pyplot.axes([0.1, 0.1, 0.8, 0.8])
        fig, ax = pyplot.subplots()
        fig.set_size_inches(6, 3)

        if not isinstance(y, list):
            y = [y]

        if color is None:
            color = []
            for index in range(len(y)):
                color.append(random.choice(list(colors.get_named_colors_mapping().keys())))

        elif not isinstance(color, list):
            color = [color]

        if legend is not None and not isinstance(legend, list):
            legend = [legend]

        if len(y) == 0:
            return

        # x = range(len(y[0]))
        x = numpy.arange(len(y[0]))
        width = .8 / len(y)

        if color is None:
            color = ["blue"]

        legend_handles = []

        # print(color)

        # Add the graph data
        for (index, val) in enumerate(y):

            ax.bar(x + ((index * width) + (width / 2) - 0.4), val, width, color=color[index])

            if legend is not None:
                patch = patches.Patch(color=color[index], label=legend[index])
                legend_handles.append(patch)

        if len(legend_handles) > 0:
            pyplot.legend(handles=legend_handles)

        # Make the labels look nicer.
        ax.set_title(title, fontdict={'weight': 'bold', 'size': 'large'})
        ax.set_xticklabels(x_tick_labels, fontdict={'style': 'italic'}, rotation=45)
        ax.set_xlabel(x_label, fontdict={'weight': 'bold'})
        ax.set_ylabel(y_label, fontdict={'weight': 'bold'})
        ax.set_xticks(x)

        # Make the x-axis labels light up correctly
        for tick in ax.xaxis.get_majorticklabels():
            tick.set_horizontalalignment("right")

        # Save the graph to file
        pyplot.savefig(outfile, bbox_inches='tight')

        pyplot.close(fig)

    @staticmethod
    def create_scatter_plot(x, y, title, x_tick_labels, x_label, y_label, outfile, color="blue", legend=None):
        # ax = pyplot.axes([0.1, 0.1, 0.8, 0.8])
        fig, ax = pyplot.subplots()
        fig.set_size_inches(6, 3)

        width = .5

        if not isinstance(y, list):
            y = [y]
        if not isinstance(x, list):
            x = [x]

        if color is None:
            color = []
            for index in range(len(y)):
                color.append(random.choice(list(colors.get_named_colors_mapping().keys())))

        elif not isinstance(color, list):
            color = [color]

        if legend is not None and not isinstance(legend, list):
            legend = [legend]

        if len(y) == 0:
            return

        legend_handles = []

        # Add the graph data
        for (index, val) in enumerate(y):
            ax.scatter(x[index], val, width, color=color[index])

            if legend is not None:
                patch = patches.Patch(color=color[index], label=legend[index])
                legend_handles.append(patch)

        if len(legend_handles) > 0:
            pyplot.legend(handles=legend_handles)

        # Make the labels look nicer.
        ax.set_title(title, fontdict={'weight': 'bold', 'size': 'large'})
        ax.set_xlabel(x_label, fontdict={'weight': 'bold'})
        ax.set_ylabel(y_label, fontdict={'weight': 'bold'})
        # ax.set_xticks(x)

        # Set the x axis labels
        if x_tick_labels is not None:
            ax.set_xticklabels(x_tick_labels, fontdict={'style': 'italic'})

        # for tick in ax.xaxis.get_majorticklabels():
        #     tick.set_horizontalalignment("right")

        # Save the graph to file
        pyplot.savefig(outfile, bbox_inches='tight')

        pyplot.close(fig)

    # (Graph 1) Distribution of total read count per individual
    @staticmethod
    def read_counts_per_individ(read_data, outfile, color=None, legend=None):
        start = time.time()

        #  Distribution of total read counts per individual
        title = "Distribution of total read counts per individual"
        y_label = "Number of individuals"
        x_label = "Total read count per individual"
        x_tick_labels = ["less than 1 million", "1 million to 1.25 million", "1.5 million to 1.75 million",
                         "1.75 million to 2 million", "2 million to 2.25 million", "2.25 million to 2.5 million",
                         "greater than 2.5 million"]

        if not isinstance(read_data, list):
            read_data = [read_data]

        y_data = []

        for graph_data in read_data:
            if len(graph_data) == 0:
                continue

            counts = {k: v["calls"] for (k, v) in graph_data.items()}

            individ_counts = list(zip(*counts.values()))
            counts_array = numpy.asarray(individ_counts)

            individ_read_counts = (numpy.sum(counts_array, axis=(1, 2)) / 2.0).tolist()

            # individReadCounts = [sum(sum(tuple) for tuple in cnts) for cnts in individCounts]

            graph_y_data = [0, 0, 0, 0, 0, 0, 0, 0]
            for count in individ_read_counts:
                if count <= 1000000:
                    graph_y_data[0] += 1
                elif count <= 1250000:
                    graph_y_data[1] += 1
                elif count <= 1500000:
                    graph_y_data[2] += 1
                elif count <= 1750000:
                    graph_y_data[3] += 1
                elif count <= 2000000:
                    graph_y_data[4] += 1
                elif count <= 2250000:
                    graph_y_data[5] += 1
                elif count <= 2500000:
                    graph_y_data[6] += 1
                else:
                    graph_y_data[7] += 1

            y_data.append(graph_y_data)

        Graphs.create_bar_graph(y_data, title, x_tick_labels, x_label, y_label, outfile, color=color, legend=legend)

        print(title + ": " + str(round((time.time() - start), 2)) + "s")

    # (Graph 2) Distribution of average reads per SNP
    @staticmethod
    def avg_reads_per_snp(read_data, outfile, color=None, legend=None):
        start = time.time()

        #  Distribution of total read counts per individual
        title = "Distribution of average reads per SNP"
        y_label = "Number of SNPs"
        x_label = "Average number of reads per SNP"
        x_tick_labels = ["<5", "5 to 10", "10 to 15", "15 to 20", "20 to 25", "25 to 30", "30 to 35", ">35"]

        if not isinstance(read_data, list):
            read_data = [read_data]

        y_data = []

        for graph_data in read_data:
            if len(graph_data) == 0:
                continue

            # Find the average reads per SNP
            counts = {k: v["calls"] for (k, v) in graph_data.items()}

            counts_array = numpy.asarray(list(counts.values()))
            reads_per_snp = numpy.sum(counts_array, axis=(1, 2)) / 2.0
            num_samples = numpy.shape(counts_array)[1]
            avg_reads_per_snp = (reads_per_snp / num_samples).tolist()

            # Convert actual values to a count of SNP's in each percentage range
            graph_y_data = [0, 0, 0, 0, 0, 0, 0, 0]
            for count in avg_reads_per_snp:
                if count <= 5:
                    graph_y_data[0] += 1
                elif count <= 10:
                    graph_y_data[1] += 1
                elif count <= 15:
                    graph_y_data[2] += 1
                elif count <= 20:
                    graph_y_data[3] += 1
                elif count <= 25:
                    graph_y_data[4] += 1
                elif count <= 30:
                    graph_y_data[5] += 1
                elif count <= 35:
                    graph_y_data[6] += 1
                else:
                    graph_y_data[7] += 1

            y_data.append(graph_y_data)

        Graphs.create_bar_graph(y=y_data, title=title, x_tick_labels=x_tick_labels, x_label=x_label, y_label=y_label,
                                outfile=outfile, color=color, legend=legend)

        print(title + ": " + str(round((time.time() - start), 2)) + "s")

    # (Graph 3) Distribution of call rates across SNPs
    @staticmethod
    def call_rates_across_snp(data, outfile, color=None, legend=None):
        start = time.time()

        #  Distribution of total read counts per individual
        title = "Distribution of call rates across SNPs"
        y_label = "Number of SNPs"
        x_label = "Call rate per SNP"
        x_tick_labels = ["<50%", "50 to 55%", "55 to 60%", "60 to 65%", "70 to 75%", "75 to 80%", "80 to 85%",
                         "85 to 90%", "90 to 95%", "95 to 100%"]

        if not isinstance(data, list):
            data = [data]

        y_data = []

        for graph_data in data:
            if len(graph_data) == 0:
                continue

            # Find the call rates per SNP
            calls = [v["calls"] for (k, v) in graph_data.items()]
            calls_array = numpy.asarray(list(calls))

            # blanks = [(calls == "-").sum() for calls in callsArray]
            positive_calls = [((calls == "1") | (calls == "2")).sum() for calls in calls_array]
            tot_calls = [(calls != "-").sum() for calls in calls_array]

            call_rate_per_snp = [pos / tot for pos, tot in zip(positive_calls, tot_calls)]

            # Convert actual call rates into counts in each percentage range.
            graph_y_data = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            for count in call_rate_per_snp:
                if count <= .50:
                    graph_y_data[0] += 1
                elif count <= .55:
                    graph_y_data[1] += 1
                elif count <= .60:
                    graph_y_data[2] += 1
                elif count <= .65:
                    graph_y_data[3] += 1
                elif count <= .70:
                    graph_y_data[4] += 1
                elif count <= .75:
                    graph_y_data[5] += 1
                elif count <= .80:
                    graph_y_data[6] += 1
                elif count <= .85:
                    graph_y_data[7] += 1
                elif count <= .90:
                    graph_y_data[8] += 1
                elif count <= .95:
                    graph_y_data[9] += 1
                else:
                    graph_y_data[10] += 1

            y_data.append(graph_y_data)

        Graphs.create_bar_graph(y=y_data, title=title, x_tick_labels=x_tick_labels, x_label=x_label, y_label=y_label,
                                outfile=outfile, color=color, legend=legend)

        print(title + ": " + str(round((time.time() - start), 2)) + "s")

    # (Graph 4) Distribution of call rates across individuals
    @staticmethod
    def call_rate_across_individ(data, outfile, color=None, legend=None):
        start = time.time()

        title = "Distribution of call rates across individuals"
        y_label = "Number of individuals"
        x_label = "Call rate per individual"
        x_tick_labels = ["<50%", "50 to 55%", "55 to 60%", "60 to 65%", "70 to 75%", "75 to 80%", "80 to 85%",
                         "85 to 90%", "90 to 95%", "95 to 100%"]

        if not isinstance(data, list):
            data = [data]

        y_data = []

        for graph_data in data:
            if len(graph_data) == 0:
                continue

            # Find the call rates per individual
            calls = [v["calls"] for (k, v) in graph_data.items()]
            individ_calls = list(zip(*calls))
            calls_array = numpy.asarray(list(individ_calls))

            # blanks = [(calls == "-").sum() for calls in callsArray]
            positive_calls = [((calls == "1") | (calls == "2")).sum() for calls in calls_array]
            tot_calls = [(calls != "-").sum() for calls in calls_array]

            call_rate_per_snp = [pos / tot for pos, tot in zip(positive_calls, tot_calls)]

            # Convert to counts within eac call rate percentage range.
            graph_y_data = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            for count in call_rate_per_snp:
                if count <= .50:
                    graph_y_data[0] += 1
                elif count <= .55:
                    graph_y_data[1] += 1
                elif count <= .60:
                    graph_y_data[2] += 1
                elif count <= .65:
                    graph_y_data[3] += 1
                elif count <= .70:
                    graph_y_data[4] += 1
                elif count <= .75:
                    graph_y_data[5] += 1
                elif count <= .80:
                    graph_y_data[6] += 1
                elif count <= .85:
                    graph_y_data[7] += 1
                elif count <= .90:
                    graph_y_data[8] += 1
                elif count <= .95:
                    graph_y_data[9] += 1
                else:
                    graph_y_data[10] += 1

            y_data.append(graph_y_data)

        Graphs.create_bar_graph(y=y_data, title=title, x_tick_labels=x_tick_labels, x_label=x_label, y_label=y_label,
                                outfile=outfile, color=color, legend=legend)

        print(title + ": " + str(round((time.time() - start), 2)) + "s")

    # (Graph 5) Distribution of MAF across SNPs
    @staticmethod
    def maf_across_snp(snp_maf, outfile, color=None, legend=None):
        start = time.time()

        title = "Distribution of MAF across SNPs"
        y_label = "Number of SNPs"
        x_label = "Minor Allele Frequency (MAF) of SNPs"
        x_tick_labels = ["less than 0.01", "0.01 to 0.02", "0.02 to 0.03", "0.03 to 0.04", "0.04 to 0.05",
                         "0.05 to 0.06",
                         "0.06 to 0.07", "0.07 to 0.08", "0.08 to 0.09", "0.09 to 0.1", "0.1 to 0.15", "0.15 to 0.2",
                         "0.2 to 0.3", "0.3 to 0.4", "0.4 to 0.5", "greater than 0.5"]

        if not isinstance(snp_maf, list):
            snp_maf = [snp_maf]

        y_data = []

        for graph_data in snp_maf:
            if len(graph_data) == 0:
                continue

            # Conver the actual MAF data into counts within the graphs ranges
            graph_y_data = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            for (k, maf) in graph_data.items():
                if maf <= .01:
                    graph_y_data[0] += 1
                elif maf <= .02:
                    graph_y_data[1] += 1
                elif maf <= .03:
                    graph_y_data[2] += 1
                elif maf <= .04:
                    graph_y_data[3] += 1
                elif maf <= .05:
                    graph_y_data[4] += 1
                elif maf <= .06:
                    graph_y_data[5] += 1
                elif maf <= .07:
                    graph_y_data[6] += 1
                elif maf <= .08:
                    graph_y_data[7] += 1
                elif maf <= .09:
                    graph_y_data[8] += 1
                elif maf <= .1:
                    graph_y_data[9] += 1
                elif maf <= .15:
                    graph_y_data[10] += 1
                elif maf <= .2:
                    graph_y_data[11] += 1
                elif maf <= .3:
                    graph_y_data[12] += 1
                elif maf <= .4:
                    graph_y_data[13] += 1
                elif maf <= .5:
                    graph_y_data[14] += 1
                else:
                    graph_y_data[15] += 1

            y_data.append(graph_y_data)

        Graphs.create_bar_graph(y=y_data, title=title, x_tick_labels=x_tick_labels, x_label=x_label, y_label=y_label,
                                outfile=outfile, color=color, legend=legend)

        print(title + ": " + str(round((time.time() - start), 2)) + "s")

    # (Graph 6) Distribution of average repeatability (repAvg) across SNPs
    @staticmethod
    def avg_rep_across_snp(data, outfile, color=None, legend=None):
        start = time.time()

        title = "Distribution of average repeatability (RepAvg) across SNPs"
        y_label = "Number of SNPs"
        x_label = "Average repeatability of SNPs"
        x_tick_labels = ["<.9", "<.92", "<.94", "<.96", "<.98", "<.99", ">.99"]

        if not isinstance(data, list):
            data = [data]

        y_data = []

        for graph_data in data:
            if len(graph_data) == 0:
                continue

            # Get the repeatability averages
            rep_avgs = [v["rep_average"] for (k, v) in graph_data.items()]

            # Convert to counts within each graph range
            graph_y_data = [0, 0, 0, 0, 0, 0, 0]
            for count in rep_avgs:
                if count <= .9:
                    graph_y_data[0] += 1
                elif count <= .92:
                    graph_y_data[1] += 1
                elif count <= .94:
                    graph_y_data[2] += 1
                elif count <= .96:
                    graph_y_data[3] += 1
                elif count <= .98:
                    graph_y_data[4] += 1
                elif count <= .99:
                    graph_y_data[5] += 1
                else:
                    graph_y_data[6] += 1

            y_data.append(graph_y_data)

        Graphs.create_bar_graph(y=y_data, title=title, x_tick_labels=x_tick_labels, x_label=x_label, y_label=y_label,
                                outfile=outfile, color=color, legend=legend)

        print(title + ": " + str(round((time.time() - start), 2)) + "s")

    # (Graph 7) Distribution of Heterozygosity across SNPs
    @staticmethod
    def het_across_snp(data, outfile, color=None, legend=None):
        start = time.time()

        title = "Distribution of Heterozygosity across SNPs"
        y_label = "Number of SNPs"
        x_label = "Heterozygosity per SNP"
        x_tick_labels = ["0 to 10%", "10 to 20%", "20 to 30%", "30 to 40%", "40 to 50%", "50 to 60%", "60 to 70%",
                         "70 to 80%", "80 to 90%", "90 to 100%"]

        if not isinstance(data, list):
            data = [data]

        y_data = []

        for graph_data in data:
            if len(graph_data) == 0:
                continue

            freq_hetz = [v["freq_heterozygous"] for (k, v) in graph_data.items()]

            # Convert into counts within each graph range
            graph_y_data = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            for count in freq_hetz:
                if float(count) < .10:
                    graph_y_data[0] += 1
                elif float(count) <= .20:
                    graph_y_data[1] += 1
                elif float(count) <= .30:
                    graph_y_data[2] += 1
                elif float(count) <= .40:
                    graph_y_data[3] += 1
                elif float(count) <= .50:
                    graph_y_data[4] += 1
                elif float(count) <= .60:
                    graph_y_data[5] += 1
                elif float(count) <= .70:
                    graph_y_data[6] += 1
                elif float(count) <= .80:
                    graph_y_data[7] += 1
                elif float(count) <= .90:
                    graph_y_data[8] += 1
                else:
                    graph_y_data[9] += 1

            y_data.append(graph_y_data)

        Graphs.create_bar_graph(y=y_data, title=title, x_tick_labels=x_tick_labels, x_label=x_label, y_label=y_label,
                                outfile=outfile, color=color, legend=legend)

        print(title + ": " + str(round((time.time() - start), 2)) + "s")

    # (Graph 8) Relationship between MAF and Read Count
    @staticmethod
    def maf_to_read_count(snp_maf, read_data, outfile, color=None, legend=None):
        start = time.time()

        title = "Relationship between MAF and Read Count"
        y_label = "MAF by SNP"
        x_label = "Read Count by SNP"
        x_tick_labels = None

        if not isinstance(snp_maf, list):
            snp_maf = [snp_maf]

        y_data = []
        x_data = []

        for (index, graph_data) in enumerate(snp_maf):
            if len(graph_data) == 0:
                continue

            # Get the read count values in the same order as the MAF data
            counts = [read_data[k]["calls"] for (k, v) in graph_data.items() if k in read_data]

            counts_array = numpy.asarray(counts)

            reads_per_snp = (numpy.sum(counts_array, axis=(1, 2)) / 2.0).tolist()
            x_data.append(reads_per_snp)

            # Find MAF for a SNP
            graph_y_data = [v for (k, v) in graph_data.items() if k in read_data]
            y_data.append(graph_y_data)

        Graphs.create_scatter_plot(x=x_data, y=y_data, title=title, x_tick_labels=x_tick_labels, x_label=x_label,
                                   y_label=y_label, outfile=outfile, color=color, legend=legend)

        print(title + ": " + str(round((time.time() - start), 2)) + "s")

    # (Graph 9) Relationship between Call rate and MAF
    @staticmethod
    def call_rate_to_maf(data, snp_maf, outfile, color=None, legend=None):
        start = time.time()

        title = "Relationship between Call rate and MAF"
        y_label = "MAF by SNP"
        x_label = "Call rate by SNP"
        x_tick_labels = ["50%", "60%", "70%", "80%", "90%", "100%", "110%"]

        if not isinstance(snp_maf, list):
            snp_maf = [snp_maf]
        if not isinstance(data, list):
            data = [data]

        y_data = []
        x_data = []

        for (index, graph_data) in enumerate(data):
            if len(graph_data) == 0:
                continue

            # Get the call rate by SNP
            calls = [v["calls"] for (k, v) in graph_data.items()]
            calls_array = numpy.asarray(list(calls))

            # blanks = [(calls == "-").sum() for calls in callsArray]
            positive_calls = [((calls == "1") | (calls == "2")).sum() for calls in calls_array]
            tot_calls = [(calls != "-").sum() for calls in calls_array]

            call_rate_per_snp = [pos / tot for pos, tot in zip(positive_calls, tot_calls)]
            x_data.append(call_rate_per_snp)

            # Find MAF for a SNP
            graph_y_data = list(snp_maf[index].values())
            y_data.append(graph_y_data)

        Graphs.create_scatter_plot(x=x_data, y=y_data, title=title, x_tick_labels=x_tick_labels, x_label=x_label,
                                   y_label=y_label, outfile=outfile, color=color, legend=legend)

        print(title + ": " + str(round((time.time() - start), 2)) + "s")
