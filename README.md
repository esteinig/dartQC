##DartQC
### Quality Control Pipeline for Diversity Array Technology (DArT) 

![](https://github.com/esteinig/dartQC/blob/master/dartQC.png)

**Dependencies**

Install in Ubuntu: 

`sudo apt-get install python-numpy python-biopython cd-hit`

If you are using this pipeline you are likely doing a bit of Bioinformatics. If you are also using Windows, consider switching to Linux or install a Virtual Machine. This program is *not* compatible with Windows.

* Numpy
* BioPython
* CDHIT-EST
* Python >= v.3.4

**Command Line Usage**

`python dart_qc.py --help`

Parameters:

```
Required:

-i, --input           Input file from DArT (.csv) ['']

Defaults:

-f, --format          Row format, currently only double-row ['double']
-t, --type            Marker type, currently only SNPs ['snp']
-o, --out             Output file type, one of: 'plink' or 'structure' ['plink']
-p, --pop             File specifying ID and Population (two columns with header, comma-delimited) ['']

--maf                 Filter markers by minor allele frequency <= [0.02]
--call                Filter markers by call rate <= [0.70]
--rep                 Filter markers by average replication statistic <= [0.95]

--sequence-identity   Filter reference allele using CDHIT-EST, identity threshold [0.95]
--identity-selector   Statistic for sequence identity picks: 'call_rate', 'maf' or 'rep' ['maf']
--clone-selector      Statistic for duplicate clone picks: 'call_rate', 'maf' or 'rep' ['maf']

--major               Homozygous major encoding in DArT, tuple of strings [ ('1', '0') ]
--minor               Homozygous minor encoding in DArT, tuple of strings [ ('0', '1') ]
--hetero              Heterozygous encoding in DArT, tuple of strings [ ('1', '1') ]
--missing             Missing encoding in DArT, tuple of strings [ ('-', '-') ]

--data-row            Row number - start of allele calls in DArT [7]
--sample-row          Row number - sample names in DArT [6]

--id-col              Column number - Allele IDs [1]
--clone-col           Column number - Clone IDs [2]
--seq-col             Column number - Sequences [3]
--rep-col             Column number - Average Replication Statistic [17]
--call-col            Column number - Start of sample names and allele calls [18]

--project             Project name for writing output ['DartData']
--verbose             Print results to screen during QC
--keep                Keep intermediate output files (from CDHIT-EST)
```

**Coding your own QC Pipeline**

```p

# Example QC
# For all options and data dictionary keys, see annotations in the classes for DartQC

def exampleQC():
    
    dart_data = DartReader()
    dart_data.project = "Koala"
    dart_data._data_row = 9
    dart_data._sample_row = 8
    dart_data._read_count_ref_column = 15
    dart_data._read_count_snp_column = 16
    dart_data._replication_count = 18
    dart_data._sample_column = 19
    dart_data._call = 19

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
