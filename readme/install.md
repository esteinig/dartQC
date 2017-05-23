## Install DartQC

At the moment, installation is a manual process but a stable version will eventually be on PyPI or Conda. 

Clone this repository and setup the Conda virtual environment that contains Python dependencies and CD-HIT from BioConda:

```
cd ~
git clone https://github.com/esteinig/dartQC
conda env create -f ~/dartQC/env/dartqc.yaml
```

#### How do I run the script?

```
source activate dartqc
~/dartQC/dartqc.py --help
```


#### How do I add a directory to the system's `PATH`?

It's good practice to create a `bin` directory containing sym-links to executables and scripts, making them available in your `PATH`. You can do the following in a Bash shell:

```
cd ~
mkdir bin
echo 'export PATH="$PATH:~/bin"' >> ~/.bashrc
source ~/.bashrc
```

#### How do I execute the script from anywhere on the system?

You can then sym-link the script to `bin` and use it from anywhere on your system. For example if `~/bin` is in `PATH`:

```
ln -s ~/dartQC/dartqc.py ~/bin

source activate dartqc
dartqc.py --help
```

Please note that activating the conda environment requires a Bash shell instead of the default shell on HPC (TCSH).
