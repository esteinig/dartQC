import time

import numpy
from matplotlib import pyplot

from DartModules import MarkerModule


class Graphs:
    @staticmethod
    def create_bar_graph(y, title, xTickLabels, xLabel, yLabel, outFile):
        # ax = pyplot.axes([0.1, 0.1, 0.8, 0.8])
        fig, ax = pyplot.subplots()
        fig.set_size_inches(6, 3)

        x = range(len(y))
        width = .5
        rects1 = ax.bar(x, y, width, color="blue")

        ax.set_title(title, fontdict={'weight': 'bold', 'size': 'large'})
        ax.set_xticklabels(xTickLabels, fontdict={'style': 'italic'}, rotation=45)
        ax.set_xlabel(xLabel, fontdict={'weight': 'bold'})
        ax.set_ylabel(yLabel, fontdict={'weight': 'bold'})
        ax.set_xticks(x)

        for tick in ax.xaxis.get_majorticklabels():
            tick.set_horizontalalignment("right")

        pyplot.savefig(outFile, bbox_inches='tight')

    @staticmethod
    def create_scatter_plot(x, y, title, xTickLabels, xLabel, yLabel, outFile):
        # ax = pyplot.axes([0.1, 0.1, 0.8, 0.8])
        fig, ax = pyplot.subplots()
        fig.set_size_inches(6, 3)

        width = .5
        rects1 = ax.scatter(x, y, width, color="blue")

        ax.set_title(title, fontdict={'weight': 'bold', 'size': 'large'})
        ax.set_xlabel(xLabel, fontdict={'weight': 'bold'})
        ax.set_ylabel(yLabel, fontdict={'weight': 'bold'})
        # ax.set_xticks(x)

        if (xTickLabels is not None):
            ax.set_xticklabels(xTickLabels, fontdict={'style': 'italic'})

        # for tick in ax.xaxis.get_majorticklabels():
        #     tick.set_horizontalalignment("right")

        pyplot.savefig(outFile, bbox_inches='tight')

    # (Graph 1) Distribution of total read count per individual
    @staticmethod
    def read_counts_per_individ(read_data, outfile):
        #  Distribution of total read counts per individual
        title = "Distribution of total read counts per individual"
        yLabel = "Number of individuals"
        xLabel = "Total read count per individual"
        xTickLabels = ["less than 1 million", "1 million to 1.25 million", "1.5 million to 1.75 million",
                       "1.75 million to 2 million", "2 million to 2.25 million", "2.25 million to 2.5 million",
                       "greater than 2.5 million"]

        start = time.time()

        counts = {k: v["calls"] for (k, v) in read_data.items()}

        individCounts = list(zip(*counts.values()))
        countsArray = numpy.asarray(individCounts)

        individReadCounts = (numpy.sum(countsArray, axis=(1, 2)) / 2.0).tolist()

        # individReadCounts = [sum(sum(tuple) for tuple in cnts) for cnts in individCounts]

        yData = [0, 0, 0, 0, 0, 0, 0, 0]
        for count in individReadCounts:
            if count <= 1000000:
                yData[0] += 1
            elif count <= 1250000:
                yData[1] += 1
            elif count <= 1500000:
                yData[2] += 1
            elif count <= 1750000:
                yData[3] += 1
            elif count <= 2000000:
                yData[4] += 1
            elif count <= 2250000:
                yData[5] += 1
            elif count <= 2500000:
                yData[6] += 1
            else:
                yData[7] += 1
        #
        Graphs.create_bar_graph(yData, title, xTickLabels, xLabel, yLabel, outfile)

        print(title + ": " + str(round((time.time() - start), 2)) + "s")

    # (Graph 2) Distribution of average reads per SNP
    @staticmethod
    def avg_reads_per_snp(read_data, outfile):
        #  Distribution of total read counts per individual
        title = "Distribution of average reads per SNP"
        yLabel = "Number of SNPs"
        xLabel = "Average number of reads per SNP"
        xTickLabels = ["<5", "5 to 10", "10 to 15", "15 to 20", "20 to 25", "25 to 30", "30 to 35", ">35"]

        start = time.time()

        counts = {k: v["calls"] for (k, v) in read_data.items()}

        countsArray = numpy.asarray(list(counts.values()))
        readsPerSNP = numpy.sum(countsArray, axis=(1, 2)) / 2.0
        numSamples = numpy.shape(countsArray)[1]
        avgReadsPerSNP = (readsPerSNP / numSamples).tolist()

        yData = [0, 0, 0, 0, 0, 0, 0, 0]
        for count in avgReadsPerSNP:
            if count <= 5:
                yData[0] += 1
            elif count <= 10:
                yData[1] += 1
            elif count <= 15:
                yData[2] += 1
            elif count <= 20:
                yData[3] += 1
            elif count <= 25:
                yData[4] += 1
            elif count <= 30:
                yData[5] += 1
            elif count <= 35:
                yData[6] += 1
            else:
                yData[7] += 1

        Graphs.create_bar_graph(yData, title, xTickLabels, xLabel, yLabel, outfile)

        print(title + ": " + str(round((time.time() - start), 2)) + "s")

    # (Graph 3) Distribution of call rates across SNPs
    @staticmethod
    def call_rates_across_snp(data, read_data, outfile):
        #  Distribution of total read counts per individual
        title = "Distribution of call rates across SNPs"
        yLabel = "Number of SNPs"
        xLabel = "Call rate per SNP"
        xTickLabels = ["<50%", "50 to 55%", "55 to 60%", "60 to 65%", "70 to 75%", "75 to 80%", "80 to 85%",
                       "85 to 90%", "90 to 95%", "95 to 100%"]

        start = time.time()

        calls = [v["calls"] for (k, v) in data.items()]
        callsArray = numpy.asarray(list(calls))

        # blanks = [(calls == "-").sum() for calls in callsArray]
        positiveCalls = [((calls == "1") | (calls == "2")).sum() for calls in callsArray]
        totCalls = [(calls != "-").sum() for calls in callsArray]

        callRatePerSNP = [pos / tot for pos, tot in zip(positiveCalls, totCalls)]

        yData = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        for count in callRatePerSNP:
            if count <= .50:
                yData[0] += 1
            elif count <= .55:
                yData[1] += 1
            elif count <= .60:
                yData[2] += 1
            elif count <= .65:
                yData[3] += 1
            elif count <= .70:
                yData[4] += 1
            elif count <= .75:
                yData[5] += 1
            elif count <= .80:
                yData[6] += 1
            elif count <= .85:
                yData[7] += 1
            elif count <= .90:
                yData[8] += 1
            elif count <= .95:
                yData[9] += 1
            else:
                yData[10] += 1

        Graphs.create_bar_graph(yData, title, xTickLabels, xLabel, yLabel, outfile)

        print(title + ": " + str(round((time.time() - start), 2)) + "s")

    # (Graph 4) Distribution of call rates across individuals
    @staticmethod
    def call_rate_across_individ(data, read_data, outfile):
        title = "Distribution of call rates across individuals"
        yLabel = "Number of individuals"
        xLabel = "Call rate per individual"
        xTickLabels = ["<50%", "50 to 55%", "55 to 60%", "60 to 65%", "70 to 75%", "75 to 80%", "80 to 85%",
                       "85 to 90%", "90 to 95%", "95 to 100%"]

        start = time.time()

        calls = [v["calls"] for (k, v) in data.items()]
        individCalls = list(zip(*calls))
        callsArray = numpy.asarray(list(individCalls))

        # blanks = [(calls == "-").sum() for calls in callsArray]
        positiveCalls = [((calls == "1") | (calls == "2")).sum() for calls in callsArray]
        totCalls = [(calls != "-").sum() for calls in callsArray]

        callRatePerSNP = [pos / tot for pos, tot in zip(positiveCalls, totCalls)]

        yData = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        for count in callRatePerSNP:
            if count <= .50:
                yData[0] += 1
            elif count <= .55:
                yData[1] += 1
            elif count <= .60:
                yData[2] += 1
            elif count <= .65:
                yData[3] += 1
            elif count <= .70:
                yData[4] += 1
            elif count <= .75:
                yData[5] += 1
            elif count <= .80:
                yData[6] += 1
            elif count <= .85:
                yData[7] += 1
            elif count <= .90:
                yData[8] += 1
            elif count <= .95:
                yData[9] += 1
            else:
                yData[10] += 1

        Graphs.create_bar_graph(yData, title, xTickLabels, xLabel, yLabel, outfile)

        print(title + ": " + str(round((time.time() - start), 2)) + "s")

    # (Graph 5) Distribution of MAF across SNPs
    @staticmethod
    def maf_across_snp(snpMAF, outfile):
        title = "Distribution of MAF across SNPs"
        yLabel = "Number of SNPs"
        xLabel = "Minor Allele Frequency (MAF) of SNPs"
        xTickLabels = ["less than 0.01", "0.01 to 0.02", "0.02 to 0.03", "0.03 to 0.04", "0.04 to 0.05", "0.05 to 0.06",
                       "0.06 to 0.07", "0.07 to 0.08", "0.08 to 0.09", "0.09 to 0.1", "0.1 to 0.15", "0.15 to 0.2",
                       "0.2 to 0.3", "0.3 to 0.4", "0.4 to 0.5", "greater than 0.5"]
        start = time.time()

        yData = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        for (k, maf) in snpMAF.items():
            if maf <= .01:
                yData[0] += 1
            elif maf <= .02:
                yData[1] += 1
            elif maf <= .03:
                yData[2] += 1
            elif maf <= .04:
                yData[3] += 1
            elif maf <= .05:
                yData[4] += 1
            elif maf <= .06:
                yData[5] += 1
            elif maf <= .07:
                yData[6] += 1
            elif maf <= .08:
                yData[7] += 1
            elif maf <= .09:
                yData[8] += 1
            elif maf <= .1:
                yData[9] += 1
            elif maf <= .15:
                yData[10] += 1
            elif maf <= .2:
                yData[11] += 1
            elif maf <= .3:
                yData[12] += 1
            elif maf <= .4:
                yData[13] += 1
            elif maf <= .5:
                yData[14] += 1
            else:
                yData[15] += 1

        Graphs.create_bar_graph(yData, title, xTickLabels, xLabel, yLabel, outfile)

        print(title + ": " + str(round((time.time() - start), 2)) + "s")

    # (Graph 6) Distribution of average repeatability (repAvg) across SNPs
    @staticmethod
    def avg_rep_across_snp(data, outfile):
        title = "Distribution of average repeatability (RepAvg) across SNPs"
        yLabel = "Number of SNPs"
        xLabel = "Average repeatability of SNPs"
        xTickLabels = ["<.9", "<.92", "<.94", "<.96", "<.98", "<.99", ">.99"]

        start = time.time()

        repAvgs = [v["rep_average"] for (k, v) in data.items()]

        yData = [0, 0, 0, 0, 0, 0, 0]
        for count in repAvgs:
            if count <= .9:
                yData[0] += 1
            elif count <= .92:
                yData[1] += 1
            elif count <= .94:
                yData[2] += 1
            elif count <= .96:
                yData[3] += 1
            elif count <= .98:
                yData[4] += 1
            elif count <= .99:
                yData[5] += 1
            else:
                yData[6] += 1

        Graphs.create_bar_graph(yData, title, xTickLabels, xLabel, yLabel, outfile)

        print(title + ": " + str(round((time.time() - start), 2)) + "s")

    # (Graph 7) Distribution of Heterozygosity across SNPs
    @staticmethod
    def het_across_snp(data, outfile):
        title = "Distribution of Heterozygosity across SNPs"
        yLabel = "Number of SNPs"
        xLabel = "Heterozygosity per SNP"
        xTickLabels = ["0 to 10%", "10 to 20%", "20 to 30%", "30 to 40%", "40 to 50%", "50 to 60%", "60 to 70%",
                       "70 to 80%", "80 to 90%", "90 to 100%"]

        start = time.time()

        freqHetz = [v["freq_heterozygous"] for (k, v) in data.items()]

        yData = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        for count in freqHetz:
            if float(count) < .10:
                yData[0] += 1
            elif float(count) <= .20:
                yData[1] += 1
            elif float(count) <= .30:
                yData[2] += 1
            elif float(count) <= .40:
                yData[3] += 1
            elif float(count) <= .50:
                yData[4] += 1
            elif float(count) <= .60:
                yData[5] += 1
            elif float(count) <= .70:
                yData[6] += 1
            elif float(count) <= .80:
                yData[7] += 1
            elif float(count) <= .90:
                yData[8] += 1
            else:
                yData[9] += 1

        Graphs.create_bar_graph(yData, title, xTickLabels, xLabel, yLabel, outfile)

        print(title + ": " + str(round((time.time() - start), 2)) + "s")

    # (Graph 8) Relationship between MAF and Read Count
    @staticmethod
    def maf_to_read_count(snpMAF, read_data, outfile):
        title = "Relationship between MAF and Read Count"
        yLabel = "MAF by SNP"
        xLabel = "Read Count by SNP"
        xTickLabels = None

        start = time.time()

        # Get the read count values in the same order as the MAF data
        counts = [read_data[k]["calls"] for (k, v) in snpMAF.items() if k in read_data]

        countsArray = numpy.asarray(counts)
        readsPerSNP = (numpy.sum(countsArray, axis=(1, 2)) / 2.0).tolist()

        # Find MAF for a SNP
        yData = [v for (k, v) in snpMAF.items() if k in read_data]

        Graphs.create_scatter_plot(readsPerSNP, yData, title, xTickLabels, xLabel, yLabel, outfile)

        print(title + ": " + str(round((time.time() - start), 2)) + "s")

    # (Graph 9) Relationship between Call rate and MAF
    @staticmethod
    def call_rate_to_maf(data, snpMAF, outfile):
        title = "Relationship between Call rate and MAF"
        yLabel = "MAF by SNP"
        xLabel = "Call rate by SNP"
        xTickLabels = ["50%", "60%", "70%", "80%", "90%", "100%", "110%"]

        start = time.time()

        # Get the call rate by SNP
        calls = [v["calls"] for (k, v) in data.items()]
        callsArray = numpy.asarray(list(calls))

        # blanks = [(calls == "-").sum() for calls in callsArray]
        positiveCalls = [((calls == "1") | (calls == "2")).sum() for calls in callsArray]
        totCalls = [(calls != "-").sum() for calls in callsArray]

        callRatePerSNP = [pos / tot for pos, tot in zip(positiveCalls, totCalls)]

        # Find MAF for a SNP
        yData = list(snpMAF.values())

        Graphs.create_scatter_plot(callRatePerSNP, yData, title, xTickLabels, xLabel, yLabel, outfile)

        print(title + ": " + str(round((time.time() - start), 2)) + "s")