## DartQC
### Quality Control Pipeline

Command line pipeline to facilitate quality control of SNP data from Diversity Array Technologies (DArT). This version is written to be user-friendly and executable on the HPC. I made the data input less flexible than before, aiming to provide a user-friendly input for basic quality control. If you need to filter on parameters provided by DArT, please contact me. 

Currently dartQC integrates:

- preprocessing using raw and called SNPs from DArT
- basic filters:
  - minor allele frequency
  - replication average
  - hardy-weinberg
  - call rate
- redundancy filtering:
  - duplicate clones
  - sequence similarity
- output for PLINK

This is roughly what's going on:

<p align="center">
 <img src="https://github.com/esteinig/dartQC/blob/master/workflow.png">
</p>

#### Dependencies

- Python 3.5
- Conda

#### How to use DartQC

This section provides a brief guide of how to install and use the pipeline components. This assumes you are using a Bash shell on a local Linux system or the HPC. There may be some trouble on Zodiac (JCU's HPC) which may still be using TCSH as default login shell. If you are unsure follow the provided guide.

#### Installation

Install by cloning this repository and installing the virtual environment that will take care of installing all other library dependencies, as well as CD-HIT.

```
git clone https://github.com/esteinig/dartQC
conda env create -f ./dartQC/env/dartqc.yaml
```

#### Data Preparation

Essentially, you need to make sure the pipeline knows where to find the right data. This is done via a configuration file (JON) which tells the program the indices (non-pythonic, starting from 1) of necessary columns and rows for raw and call data sets. You can either provide a manually curated scheme file (`--call_scheme` or `--raw_scheme`) or attempt to guess the format from the file via the command `prepare`.


