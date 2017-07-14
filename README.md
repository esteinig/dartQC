## DartQC
### Quality Control Pipeline

Command line pipeline to facilitate quality control of SNP data from Diversity Array Technologies (DArT). This version is a re-write of the original scripts aiming to be somewhat more user-friendly and executable on JCU's HPC.

#### Install

Requires conda package manager, e.g. `miniconda` see [Install DartQC](https://github.com/esteinig/dartQC/blob/master/docs/install.md).

```
conda install dartqc -c bioconda -c esteinig
```

#### How to use DartQC

This section provides a brief guide of how to install and use DartQC, assuming a Bash shell and Unix. If you are using the pipeline on JCU's HPC (Zodiac) please read the relevant sections in [Install DartQC](https://github.com/esteinig/dartQC/blob/master/docs/install.md) and [Task: pbs](https://github.com/esteinig/dartQC/blob/master/docs/task.pbs.md)

1. [Install DartQC](https://github.com/esteinig/dartQC/blob/master/docs/install.md)
2. [Task: prepare](https://github.com/esteinig/dartQC/blob/master/docs/task.prepare.md)
3. [Task: validate](https://github.com/esteinig/dartQC/blob/master/docs/task.validate.md)
4. [Task: process](https://github.com/esteinig/dartQC/blob/master/docs/task.process.md)
5. [Task: filter](https://github.com/esteinig/dartQC/blob/master/docs/task.filter.md)
6. [Task: pbs](https://github.com/esteinig/dartQC/blob/master/docs/task.pbs.md)

#### Tasks

DartQC has a hierarchical parser structure that allows you to set global options and execute a task (prepare, process, filter) with its own specific arguments:

```
dartqc [--help] [--project] [--output_path] [--pop] task

Arguments:

--project, -p          output prefix
--output_path, -o      output directory
--populations, --pop   csv file with header: id, population

Tasks:

dartqc prepare --help
dartqc validate --help
dartqc process --help
dartqc filter -- help

Support Tasks:

dartqc install --help
dartqc pbs --help
```

Global arguments are specified before the command for a task, like this:

**`dartqc`**`--project example --output_path ./example`**`prepare`**`--file example_data.csv`


#### Quick Start

Example workflow without pre-processing from Excel or CSV:

```
# CSV
dartqc prepare --file example.csv

# Excel
dartqc prepare --file example.xlsx --sheet double_row_snps

dartqc filter --call example.csv --call_scheme example_scheme.json --maf 0.02 --clusters
```

Example workflow with pre-processing:

```
dartqc prepare --file calls.csv
dartqc prepare --file raw.csv

dartqc process -c calls.csv --call_scheme calls_scheme.json -r raw.csv --raw_scheme raw_scheme.json --read_sum 7

dartqc filter --processed . --maf 0.02 --call_rate 0.7 --duplicates --clusters
```

---

<p align="center">
 <img src="https://github.com/esteinig/dartQC/blob/master/workflow.png" height="768" width="768">
</p>

---

### Contact

If you find any bugs, please submit an issue for this repository on GitHub.



