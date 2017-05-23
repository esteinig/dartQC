## DartQC
### Quality Control Pipeline

Command line pipeline to facilitate quality control of SNP data from Diversity Array Technologies (DArT). This version is written to be user-friendly and executable on the HPC. Currently dartQC integrates:

- preprocessing using raw and called SNPs from DArT
- basic filters:
  - `minor allele frequency`
  - `replication average`
  - `hardy-weinberg`
  - `call rate`
- redundancy filtering:
  - `duplicate clones`
  - `sequence similarity`
- output for PLINK

This is roughly what's going on:

<p align="center">
 <img src="https://github.com/esteinig/dartQC/blob/master/workflow.png">
</p>

#### Dependencies

DartQC is written for local Unix systems or JCU's HPC. It relies on the package and environment manager [Conda]() with a code base in Python.

- Python > 3.5
- Miniconda or Anaconda for Python 3

Briefly:

```
wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

#### How to use DartQC

This section provides a brief guide of how to install and use the pipeline components. This assumes you are using a Bash shell on a local Unix system or the JCU's HPC (Zodiac). There may be some trouble on Zodiac due to the default login shell (TCSH instead of Bash). If you are unsure follow the guide to setting up on JCU's [Zodiac]().

1. [Install DartQC]()
2. [Preparation]()
3. [Preprocessing]()
4. [Filtering]()
5. [DartQC on Zodiac]()

#### Tasks

DartQC has a hierarchical parser structure that allows you to set global options and execute a task (prepare, process, filter) with its own specific arguments:

```
# global options
dartqc --help

# task options
dartqc prepare --help
dartqc process --help
dartqc filter -- help
```

Global arguments are specified *before* the command for a task, like this:

**`dartqc`**`--prefix example --output_path ./example`**`prepare`**`--file example_data.csv`


#### Quick Start

Example workflow without pre-processing:

```
source activate dartqc

dartqc -p example -o ./example prepare --file example.csv
dartqc -p example -o ./example filter --call example.csv --call_scheme example_scheme.json --maf 0.02 --clusters

source deactivate
```

Example workflow with pre-processing:

```
source activate dartqc

dartqc -p example -o ./example prepare --file calls.csv
dartqc -p example -o ./example prepare --file raw.csv

dartqc -p example -o ./example process --calls calls.csv --call_scheme calls_scheme.json --raw raw.csv --raw_scheme raw_scheme.json --read_threshold 7

dartqc -p example -o ./example filter --processed ./example --maf 0.02 --call_rate 0.7 --duplicates --clusters

source deactivate
```

### Contact

If you find any bugs, please submit an issue for this repository on GitHub.



