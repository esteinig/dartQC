###DartQC Workshop 
---

####Tutorial

This workshop tutorial is meant to walk you through some of the steps of setting up an analysis environment and running *DartQC*. Mainly, this will involve setting up a Virtual Machine running Linux (Ubuntu), installation of dependencies and the (inevitable) pitfalls of running the program. Please do not attempt to install *DartQC* under Windows, unless you are rather experienced. It is much easier, faster and more worthwhile in the long run to set up a VM.

---
####1. Virtual Machine running Ubuntu
---

Alright let's get started! Since you are likely using Windows as a primary OS, we need to install a virtual computer that runs a Linux OS. Ubuntu is one of many flavours of Linux, but has very goofd support, functionality and is closest to what you are used to under Windows.

1. Download and install [VirtualBox](https://www.virtualbox.org/)
2. Download the [Ubuntu 14.04 LTS (Long Term Support) ISO](http://www.ubuntu.com/download/desktop)

Now, there are many excellent resources out there that will show you how to set up the VM. I recommend working through the first Tutorial, since it also discusses how to share files and adjust your display sizes in the VM.

1. [Tutorial NHI](http://linus.nci.nih.gov/bdge/installUbuntu.html)
2. [Tutorial Video](https://www.youtube.com/watch?v=CkDd6jClqEE&nohtml5=False)

Now, just a couple of things to look out for when you set it up for the first time:

1. Make sure you assign enough processors to your VM (at least 2)
2. Make sure you assign enough disk space to your VM (depending on your system, but try to get around 100 GB)
3. Usually on Intel processors the Virtualization is disabled in the BIOS, if this is the case and you can't get the machine running, look at [this explanation](http://www.howtogeek.com/213795/how-to-enable-intel-vt-x-in-your-computers-bios-or-uefi-firmware/?PageSpeed=noscript).

Good job, now that you have the machine running, fire it up and log into your account, we will get to the good stuff now! 

---
####2. Installing DartQC and Dependencies
---

#####DartQC

Let's get started by opening the command line or Terminal with `Ctrl` + `Alt` + `T`.

GitHub is an excellent way for developers of programs, but also for users since you can very easily stay up-to-date with changes in a repository. This is best achieved with the command-line functionality of `git`, so we wil linstall it now with:

`sudo apt-get install git`

Now, make sure you are in your home directory and clone this repository - this will download all files onto your local system:

`git clone https://github.com/esteinig/dartQC.git`

If you want to keep up-to-date, all you now need is to navigate to the program's directory in your home `dartQC` and execute the pull command, which pulls the latest version of dartQC:

`git pull`

#####Dependencies

Alright, we now have the program on our system. It is written in Python 3 and since we are working with sequence data and clustering, we need some additional packages and install CD-HIT. 

Instead of using the standard installation of Python we will use the manager Anaconda. It keeps a separate installation of Python and has some very beautiful installation tools that make life so much easier! [Download Anaconda](https://www.continuum.io/downloads) for Python >= 3.4 and navigate to the directory of your download in the Terminal, then install:

```
cd /home/esteinig/Downloads
bash Anaconda3-4.0.0-Linux-x86_64.sh
```

The cool thing now is, you can use the `conda` installer to add packages to Python, we already have `Numpy` but need `BioPython`:

`conda install biopython`

The last thing we need is the sequence clustering package [CD-HIT](http://bioinformatics.org/cd-hit/), which is also simple to install through the Terminal:

`sudo apt-get install cd-hit`

---
#####Dependency Checks and Troubleshoots
---

It's good practice to check if everything is in order, let's first have a look if we are using the right Python installation in the OS:

`python --version`

If Anaconda has not added the new version automatically, this will likely return the pre-installed Python 2.7. If we want to run the correct version, we need to either use the complete path to the program or add it to the system variable `$PATH`. The latter tells your command line where to look for things when you execute them from the Terminal.

We first need to go to the home directory and open a file called `.bashrc` with a text editor - this file essentially loads settings for the command line, so that whenever you open it up, you will have your specified `$PATH`:

```
cd ~
gedit .bashrc
```

At the end of the file, we now add the path (replace `esteinig` with your user name, you can get the current directory path with `pwd`) to Anaconda's Python to the beginning and attach the original `$PATH` to the end:

`export PATH="/home/esteinig/Anaconda3/bin:$PATH"`

You can add additional paths by seperating them with a semicolon. Save and exit the file, don't forget to reload the file or re-open the Terminal:

`source .bashrc`

Now check the path and version and run Python 3. Within Python, try importing Bio (BioPython):

```
which python
python --version
python
import Bio
quit()
```

If no errors are returned, you are good to go after one last check for CDHIT-EST:

`cdhit-est`

---
#####Adding DartQC to `$PATH`
---

Now that you know how it works, you can also add the script itself to your `$PATH`, simply:

```
cd ~
gedit .bashrc
```

Add with your username replaced to the document:

```
export PATH="$PATH:/home/esteinig/dartQC"
```

Load the document or re-open the Terminal, then make the script executable:

```
source .bashrc
chmod 755 /home/esteinig/dartQC/dart_qc.py 
```

Now you can call the script from any directory in your Terminal:

`dart_qc.py --help`

Woop, woop.

*If you update your local repository with `pull`, you need to make the script executable again!*

---
####3. Running DartQC from the Command Line
---

Alright, let's get started - the first thing is to convert your Excel sheet to a comma-delimited file (`.csv`). At the moment, the input is limited to the double-row format and SNPs.

Get an overview of the options on the main page of this repository or type:

`python dart_qc.py --help`

#####Usage and Pitfalls

The main thing to look out for is how your data file is structured - since we wil likely see changes in the files or people want to use previous versions, the important parameters for calculations can be specified by the user, with defaults indicated in square brackets. The counting is non-pythonic, i.e. first line in the data file is line 1. **Make sure you have the correct column and row designations, error checks are not possible for variable format input!** If you do run into errors, check the error messages first, most of the time they should indicate what is wrong and usually it is related to the input format or specifications of options.

```
--data-row            Row number - start of allele calls [7]
--sample-row          Row number - sample names [6]
--pop-row             Row number - population names, zero for generic 'Pop' [0]

--id-col              Column number - Allele IDs [1]
--clone-col           Column number - Clone IDs [2]
--seq-col             Column number - Sequences [3]
--rep-col             Column number - Average Replication Statistic [17]
--call-col            Column number - Start of sample names and allele calls [18]
```

AlleleIDs must be unique, and sometimes are already given by the clone ID including the SNP, such as:

`15900384|F|0--14:T>G`

The clone IDs can (and should be) non-unique, in this case simply:

`15900384`

Currently, the double row format has a weird encoding of the alleles, where:

```
1 0   =   Homozygous Major
1 1   =   Heterozygous
0 1   =   Homozygous Minor
```

You can change these depending on your file format version with:

```
--major               Homozygous major encoding, bi-allelic string ["10"]
--minor               Homozygous minor encoding, bi-allelic string ["01"]
--hetero              Heterozygous encoding, bi-allelic string ["11"]
--missing             Missing encoding, bi-allelic string ["--"]
```

If you have populations encoded in the file, you can give the row number to extract, otherwise generic name `Pop` is assigned to the ouput files - you can also specify a separate population file, which contains two columns with headers: `ID` and `Population`, where your `ID` must match the ones present in the data file.

CD-HIT will remove duplicate clones or identical sequences and retain one sequence for each duplicate or cluster, according to the best statistic - you can change these with:

```
--identity-selector   Statistic for sequence identity picks: 'call_rate', 'maf' or 'rep' ['maf']
--clone-selector      Statistic for duplicate clone picks: 'call_rate', 'maf' or 'rep' ['maf']
```

The filters run in the order of minor allele frequency, call rate and replication average - if you want to deactivate a filter, you need to set it to -1.

```
--maf                 Filter markers by minor allele frequency <= [0.02]
--call                Filter markers by call rate <= [0.70]
--rep                 Filter markers by average replication statistic <= [0.95]
```

Finally, let's look at some examples:

#####Examples for DartQC

```
# Default:

dart_qc.py -i monodon.csv --project Monodon

# Set different filters, deactivate Average Replication Filter:

dart_qc.py -i monodon.csv --project Monodon --maf 0.5 --call 0.2 --rep -1

# Set different selectors for duplicate picks:

dart_qc.py -i monodon.csv --project Monodon --maf 0.5 --call 0.2 --rep -1 --clone-selector call_rate --identity-selector maf

# Input file with different format:

dart_qc.py -i koala.csv --project Koala --data-row 5 --sample-row 4 --pop-row 3 --id-col 2 --clone-col 1
```

Well done, this should cover the basics - you will usually notice in your output files if something has gone majorly wrong, it is always good to check the output for any program and see if it makes sense! I will add some more information on how to code your own pipeline next week.




