###DartQC Workshop 
---

####Tutorial

This more extensive workshop tutorial is meant to walk you through some of the steps of setting up an analysis environment and running *DartQC*. Mainly, this will involve setting up a Virtual Machine running Linux (Ubuntu), installation of dependencies and the (inevitable) pitfalls of running the program. Please do not attempt to install *DartQC* under Windows, unless you are rather experienced. It is much easier, faster and more worthwhile in the long run to set up a VM.

####1. Virtual Machine running Ubuntu

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

####2. Installing DartQC and Dependencies

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

#####Dependency Checks

It's good practice to check if everything is in order, let's first have a look if we are using the right Python installation in the OS:

`python --version`

This will likely return the pre-installed Python 2.7 - if we want to run the correct version, we need to either add it to 


