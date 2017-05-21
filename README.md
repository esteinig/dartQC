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

Unfortunately, the data is usually received as an Excel file and we have seen it change formats regularly. This is a nightmare for handling input data in downstream applications. Essentially, you need to make sure the pipeline knows where to find the right data in this mess. This is done via a configuration file (JON) which tells the program the indices (non-pythonic, starting from 1) of necessary columns and rows for raw and call data sets. 

At the moment this is manual (and you only need to do it once or create seperate scheme files for different data sets), but we are looking into adding a component to the pipeline that attempts to guess these files from the data.

For now, have a look at the file for the call data `dartQC/schemes/call.scheme.json`:

```json
{
  "clone_column": 1,
  "allele_column": 2,
  "sequence_column": 3,
  "replication_column": 17,
  "call_column": 18,
  "pop_row": 0,
  "sample_row": 7,
  "data_row": 8,
  "split_clone": false,
  "split_char": "|"
}
```

You can see that for each essential column and row an index is given. Data row is where the data (SNPs) begin, sample and pop row are the rows in which you have your sample (unique) and population identifiers.
