# Task: Filter

Main task to filter SNPs in DartQC.

Inputs are either the call data and scheme files with `--calls` and `--call_scheme` or the directory containing the project's (global option `--project`) pre-processed files with `--preprocessed`.

The following filters remove SNPs:
- `--maf` <= minor allele frequence, default is None
- `--call_rate` <= call rate of SNP, default is None
- `--hwe` >= p-value of Hardy-Weinberg Equilibrium, default is None
- `--rep` <= replicatation average provided by DArT, default is None

The following filters apply redundancy tests after filtering SNPs:
- `--duplicates`: remove SNPs wit duplicate `CloneID`, default False
- `--clusters`: cluster SNP `AlleleSequence` with nucleotide CD-HIT and pick one SNP per cluster by highest MAF, default False
- `--identity`: nucleotide identity by which to cluster with CD-HIT, default is 0.95 (95%)

Output are: `project_filtered.ped`, `project_filtered.map`, `project_filtered_data.json`,  `project_filtered_attr.json`

---

Filter call file in working directory by MAF and Call Rate, remove duplicate `CloneID`, cluster `AlleleSequence`:

`dartqc.py filter --calls example_calls.csv --call_scheme example_calls_scheme.json --maf 0.02 --call_rate 0.7 --duplicates --clusters`

---

Filter pre-processed files from project `preprocess` in current working directory:

`dartqc.py --project preprocess filter --processed . --maf 0.02 --call_rate 0.7

---
