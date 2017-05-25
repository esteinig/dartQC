# Task: Filter

```
dartqc filter [--help] [--processed PROCESSED_PATH] [--calls CALL_FILE]
              [--call_scheme CALL_SCHEME] [--maf MAF] [--hwe HWE]
              [--call_rate CALL_RATE] [--rep REP] [--mind MIND]
              [--mono MONO] [--mono_comparison MONO_COMP]
              [--split_clones SPLIT_CLONES] [--duplicates] [--clusters]
              [--identity IDENTITY]
              
Arguments:

--processed           input path to processed data files (project_data.json, project_attr.json)
--calls, -c           path to called read file
--call_scheme         path to call scheme json file
--maf                 filter snps <= minor allele frequency
--hwe                 filter snps <= p-value of hardy-weinberg test
--call_rate           filter snps <= call rate of snp
--rep                 filter snps <= replication average of snp
--mind                filter samples > missingness per sample
--mono                filter samples monomorphic in <mono> populations ('all', int)
--mono_comparison     filter samples monomorphic in >=, <=, == populations ('==')
--duplicates          remove snps with duplicate clone IDs
--clusters            remove snps in identical sequence clusters
--identity            remove snps in identical sequence clusters
```

Main task to filter SNPs in DartQC.

Inputs are either the call data and scheme files with `--calls` and `--call_scheme` or the directory containing the project's (global option `--project`) pre-processed files with `--preprocessed`.

The following filters remove samples:
- `--mind` > missing data per sample across all SNPs, default is None

Note: `--mind` is really slow at the moment, especially with larger sample sizes.

The following filters remove SNPs:
- `--maf` <= minor allele frequency, default is None
- `--call_rate` <= call rate of SNP, default is None
- `--hwe` <= p-value of Hardy-Weinberg Equilibrium, default is None
- `--rep` <= replicatation average provided by DArT, default is None
- `--mono` == 'all' or int, remove monoorphic snps in mono populations, needs global option --pop

The following filters apply redundancy tests after filtering SNPs:
- `--duplicates`: remove SNPs wit duplicate `CloneID`, default False
- `--clusters`: cluster SNP `AlleleSequence` with nucleotide CD-HIT and pick one SNP per cluster by highest MAF, default False
- `--identity`: nucleotide identity by which to cluster with CD-HIT, default is 0.95 (95%)

Output are: `project_filtered.ped`, `project_filtered.map`, `project_filtered_data.json`,  `project_filtered_attr.json`

---

Filter call file in working directory by MAF and Call Rate, remove duplicate `CloneID`, cluster `AlleleSequence`:

`dartqc filter --calls example_calls.csv --call_scheme example_calls_scheme.json --maf 0.02 --call_rate 0.7 --duplicates --clusters`

---

Filter pre-processed files (`dartqc_data.json`, `dartqc_attr.json`) from project `dartqc` in current working directory:

`dartqc --project preprocess filter --processed . --maf 0.02 --call_rate 0.7`

---
