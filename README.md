##DartQC
### Quality Control Pipeline

<p align="center">
 <img src="https://github.com/esteinig/dartQC/blob/master/dart_qc.png">
</p>

**Workshop @ JCU**

* [Workshop Tutorial](https://github.com/esteinig/dartQC/blob/master/workshop.md)

**Tutorials**

* [Tutorial 1 - Command Line]()
* [Tutorial 2 - Code DartQC in Python]()
* [Tutorial 3 - Shiny UI for R]()

**Dependencies**

* Python >= v.3.4
* Numpy
* BioPython
* CDHIT-EST

I highly recommend the [Anaconda](https://www.continuum.io/downloads) installer for a clean and accessible manager for Python. 

If you are using this pipeline you are likely doing a bit of Bioinformatics. If you are using Windows, consider switching to Linux or install a Virtual Machine.

Install in Anaconda:

`conda install biopython`

Install in Ubuntu: 

`sudo apt-get install cd-hit`

**Command Line Usage**

`python dart_qc.py --help`

Parameters:

```
Required:

-i, --input           Comma-delimited input file with raw calls (.csv) ['']

Defaults:

-f, --format          Row format, currently only double-row format ['double']
-t, --type            Marker type, currently only SNPs ['snp']
-o, --out             Output file type, one of: 'plink' or 'structure' ['plink']
-p, --pop             File specifying ID and Population (two columns with header, comma-delimited), 
                      see below for specification of populations in input file ['']

--maf                 Filter markers by minor allele frequency <= [0.02]
--call                Filter markers by call rate <= [0.70]
--rep                 Filter markers by average replication statistic <= [0.95]

--sequence-identity   Filter reference allele using CDHIT-EST, identity threshold [0.95]
--identity-selector   Statistic for sequence identity picks: 'call_rate', 'maf' or 'rep' ['maf']
--clone-selector      Statistic for duplicate clone picks: 'call_rate', 'maf' or 'rep' ['maf']

--major               Homozygous major encoding, bi-allelic string ["10"]
--minor               Homozygous minor encoding, bi-allelic string ["01"]
--hetero              Heterozygous encoding, bi-allelic string ["11"]
--missing             Missing encoding, bi-allelic string ["--"]

--data-row            Row number - start of allele calls [7]
--sample-row          Row number - sample names [6]
--pop-row             Row number - population names, zero for generic 'Pop' [0]

--id-col              Column number - Allele IDs [1]
--clone-col           Column number - Clone IDs [2]
--seq-col             Column number - Sequences [3]
--rep-col             Column number - Average Replication Statistic [17]
--call-col            Column number - Start of sample names and allele calls [18]

--project             Project name for writing output ['Data']
--verbose             Print results to screen during QC
--keep                Keep intermediate output files (from CDHIT-EST)
```

**Coding your own QC Pipeline**

```p

# Example QC

def exampleQC():
    
    dart_data = DartReader()
    dart_data.project = "Koala"
    dart_data.set_options(project="Koala", data_start_row=9, sample_row=8,
                          read_count_ref_col=15, read_count_snp_col=16, rep_average_col=18,
                          sample_start_col=19, call_start_col=19)
    
    dart_data.read_data("koalaInput.csv")
    dart_data.read_pops("koalaPopulations.csv")

    dart_qc = DartControl(dart_data)
    
    dart_qc.find_duplicate_clones()
    dart_qc.select_best_clones(selector="call_rate")
    dart_qc.find_identity_clusters()
    dart_qc.select_best_identity_seqs(selector="maf")

    dart_qc.filter_snps(data="total", selector="call_rate", threshold=0.70, comparison="<=")
    dart_qc.filter_snps(data="filtered", selector="maf", threshold=0.02, comparison="<=")
    dart_qc.filter_snps(data="filtered", selector="average_read_count_ref", threshold=50, comparison="<=")
    
    dart_writer = DartWriter(dart_qc)
    dart_writer.write_snps(mode='dart)
    dart_writer.write_snps(mode='plink')
    dart_writer.write_log()
    
```
