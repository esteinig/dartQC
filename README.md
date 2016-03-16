##DartQC
### Quality Control Pipeline for Diversity Array Technology (DArT) 



**Dependencies**

If you are using this pipeline you are likely doing a bit of Bioinformatics. Do not use Windows, switch to Linux / Ubuntu. Seriously.

* Numpy
* BioPython
* CDHIT-EST
* Python >= v.3.4

`sudo apt-get install python-numpy python-biopython cd-hit`

**Usage**

`python dart_qc.py --help`

Parameters:

```
Required:

-i, --input     Input file from DArT (.csv) ['']

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

