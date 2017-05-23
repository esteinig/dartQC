## DartQC
### Quality Control Pipeline

Command line pipeline to facilitate quality control of SNP data from Diversity Array Technologies (DArT). This version is a re-qrite of the original scripts aiming to be somewhat more user-friendly and executable on a HPC.

<p align="center">
 <img src="https://github.com/esteinig/dartQC/blob/master/workflow.png">
</p>

#### Dependencies

DartQC is written for local Unix systems or JCU's HPC. It relies on the package and environment manager [Conda]() with a code base in Python

- Miniconda or Anaconda for Python 3

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
dartqc --help

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



