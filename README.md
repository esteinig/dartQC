## DartQC
### Quality Control Pipeline

Command line pipeline to facilitate quality control of SNP data from Diversity Array Technologies (DArT). This version is a re-write of the original scripts aiming to be somewhat more user-friendly and executable on JCU's HPC.

<p align="center">
 <img src="https://github.com/esteinig/dartQC/blob/master/workflow.png" height="1024" width="768">
</p>

#### How to use DartQC

This section provides a brief guide of how to install and use the pipeline components. This assumes you are using a Bash shell on a local Unix system or the JCU's HPC (Zodiac).

1. [Install DartQC](https://github.com/esteinig/dartQC/blob/master/readme/install.md)
2. [Task: prepare](https://github.com/esteinig/dartQC/blob/master/readme/task.prepare.md)
3. [Task: process](https://github.com/esteinig/dartQC/blob/master/readme/task.process.md)
4. [Task: filter](https://github.com/esteinig/dartQC/blob/master/readme/task.filter.md)
5. [DartQC on Zodiac]()

#### Tasks

DartQC has a hierarchical parser structure that allows you to set global options and execute a task (prepare, process, filter) with its own specific arguments:

```
dartqc.py --help

dartqc.py prepare --help
dartqc.py process --help
dartqc.py filter -- help
```

Global arguments are specified *before* the command for a task, like this:

**`dartqc.py`**`--prefix example --output_path ./example`**`prepare`**`--file example_data.csv`


#### Quick Start

Example workflow without pre-processing from Excel or CSV:

```
source activate dartqc

# CSV
dartqc.py prepare --file example.csv
# Excel
dartqc.py prepare --file example.xlsx --sheet double_row_snps

dartqc.py filter --call example.csv --call_scheme example_scheme.json --maf 0.02 --clusters

source deactivate
```

Example workflow with pre-processing:

```
source activate dartqc

dartqc.py prepare --file calls.csv
dartqc.py prepare --file raw.csv

dartqc.py filter -c calls.csv --call_scheme calls_scheme.json -r raw.csv --raw_scheme raw_scheme.json --read_threshold 7

dartqc.py filter --processed ./example --maf 0.02 --call_rate 0.7 --duplicates --clusters

source deactivate
```

### Contact

If you find any bugs, please submit an issue for this repository on GitHub.



