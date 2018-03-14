import csv
import operator
import os
import shutil
from subprocess import call

import logging
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from install import cdhit_config

from dartqc.PipelineOptions import Filter
from dartqc.Dataset import Dataset
from dartqc.FilterResult import FilterResult
from dartqc.filters.MinorAlleleFreq import MAFFilter

log = logging.getLogger(__file__)


# TODO: Not sure this is completely finished...
class ClusterFilter(Filter):
    def get_order(self) -> int:
        return 100

    def get_cmd_type(self):
        return lambda s: [float(item.strip()) if len(item.strip()) > 0 else None for item in s[1:-1].split(',')]

    def get_name(self) -> str:
        return "cluster"

    def get_cmd_help(self) -> str:
        return "remove SNPs in similar sequence clusters (cd-hit-est, value is clustering cloeseness/distance/identity)"

    def filter(self, dataset: Dataset, identity: float, unknown_args: [], **kwargs) -> FilterResult:
        silenced = FilterResult()

        log.info("Identifying clusters using cd-hit-est")

        clusters = ClusterFilter._find_clusters(dataset, identity)

        log.info("{} clusters identified - silencing least important SNPs based on MAF".format(len(clusters)))

        removed, retained = ClusterFilter._select_by_maf(dataset, clusters)

        for snp in removed:
            silenced.silenced_snp(snp)

        return silenced

    @staticmethod
    def _find_clusters(dataset, identity=0.95, word_size=10, description_length=0):
        """
        Clusters the reference allele sequences with CDHIT-EST and parses the clusters for selecting and
        retaining best sequences.

        CD-HIT returns slightly different cluster configurations for each run due to greedy incremental algorithm,
        but little variation observed between runs in the data for P. monodon. Know thyself!

        """

        tmp_path = os.path.join(dataset.working_dir, "tmp")
        os.makedirs(tmp_path, exist_ok=True)
        fasta_path = os.path.join(tmp_path, "snp_seqs.fasta")

        ClusterFilter._write_fasta(dataset, fasta_path)

        out_path = dataset.batch_id + "_IdentityClusters_" + str(identity) + ".clstr"
        cluster_path = ClusterFilter._run_cdhit(fasta_path=fasta_path, identity=identity, word_size=word_size,
                                                description_length=description_length, out_path=out_path)

        clusters = ClusterFilter._parse_cdhit(cluster_path)

        shutil.rmtree(tmp_path, ignore_errors=True)

        return clusters

    @staticmethod
    def _write_fasta(dataset, path):
        """ Write fasta file of sequences with SNP IDs for CD-HIT. """

        snp_defs = [snp_def for snp_def in dataset.snps if snp_def.allele_id not in dataset.filtered.snps]

        seqs = [SeqRecord(Seq(snp.sequence_ref, IUPAC.unambiguous_dna), id=snp.allele_id, name="", description="") for
                snp in snp_defs]

        with open(path, "w") as fasta_file:
            SeqIO.write(seqs, fasta_file, "fasta")

    @staticmethod
    def _run_cdhit(fasta_path, identity=0.95, word_size=5, description_length=0, out_path=None):
        """ Run CDHIT-EST for sequences, install with sudo apt install cd-hit on Ubuntu """

        cdhit_path = None
        with open(cdhit_config, "r") as settings:
            cdhit_path = settings.readline()

        assert cdhit_path is not None, "cdhit path isn't configured.  Configure by running install.py with --cdhit_path <path>"

        if out_path.endswith(".clstr"):
            out_path = out_path[:-6]

        cluster_path = out_path + '.clstr'

        print("Calling cd-hit-est: " + cdhit_path + " -i " + fasta_path + " -o " + out_path + " -c " + str(identity) +
              " -n " + str(word_size) + " -d " + str(description_length))

        with open(os.devnull, "w") as devnull:
            call([cdhit_path, "-i", fasta_path, "-o", out_path, "-c", str(identity), "-n", str(word_size),
                  "-d", str(description_length)], stdout=devnull)

        return cluster_path

    @staticmethod
    def _parse_cdhit(file):
        """
        Parses the CDHIT cluster output and retains only cluster with more than one sequence for selection and
        removal from total SNPs in dictionary: {cluster_id: snp_ids}

        """

        identity_clusters = {}

        with open(file, "r") as clstr_file:
            nrows = len(list(clstr_file))

        with open(file, "r") as clstr_file:
            clstr_reader = csv.reader(clstr_file, delimiter=' ')

            cluster_id = 0
            ids = []

            row_index = 1
            for row in clstr_reader:
                if row[0] == ">Cluster":
                    if len(ids) > 1:
                        identity_clusters[cluster_id] = ids
                    ids = []
                    cluster_id += 1
                    row_index += 1
                else:
                    allele_id = ClusterFilter._find_between(row[1], ">", "...")
                    ids.append(allele_id)
                    # EOF
                    if row_index == nrows:
                        if len(ids) > 1:
                            identity_clusters[cluster_id] = ids
                    row_index += 1

        return identity_clusters

    @staticmethod
    def _find_between(s, first, last):
        """Finds the substring between two strings."""

        try:
            start = s.index(first) + len(first)
            end = s.index(last, start)
            return s[start:end]
        except ValueError:
            return ""

    @staticmethod
    def _select_by_maf(dataset, clusters):
        """ Select best markers from clusters by selector. """
        retained_sequences = []
        removed_sequences = []

        clustered_alleles = [allele_id for cluster_allele_ids in clusters.values() for allele_id in cluster_allele_ids]
        maf_values = MAFFilter.calculate_maf(dataset, False, None, clustered_alleles)
        maf_no_pops = maf_values[list(maf_values.keys())[0]]

        for cluster, allele_ids in clusters.items():
            # best_sequence = ClusterFilter._compare_entries(dataset, cluster_members)
            entries_stats = [[allele_id, maf_no_pops[allele_id]] for allele_id in allele_ids]
            entries_ranked = sorted(entries_stats, key=operator.itemgetter(1), reverse=True)

            best_sequence = entries_ranked[0][0]

            retained_sequences.append(best_sequence)
            removed_sequences += [member for member in allele_ids if member != best_sequence]

        return retained_sequences, removed_sequences

ClusterFilter()
