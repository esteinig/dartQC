## DartQC
### Quality Control Pipeline

Command line pipeline to facilitate quality control of SNP data from Diversity Array Technologies (DArT). This version is written to be user-friendly and executable on the HPC. Input complexity reduced for basic quality control. If you need to filter on parameters provided by DArT, please contact me. Currently integrates:

- preprocessing using raw and called data from DArT
- basic filters (minor allele frequency, replication average, hardy-weinberg, call rate)
- redundancy filtering (duplicate clones, sequence similarity clusters)
- output for PLINK

This is roughly what's going on:

<p align="center">
 <img src="https://github.com/esteinig/dartQC/blob/master/workflow.png">
</p>


#### How to use DartQC

...
