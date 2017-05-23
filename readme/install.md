## Install DartQC

Please note that activating the conda environment requires a Bash shell instead of the default shell on HPC (TCSH). You can check what shell you are using and if return is `-tcsh`:

```
echo $0      # if -tcsh
/bin/bash    # enter bash shell
 ```
 
#### How do I install the dependencies?

Dependencies (Python, CD-HIT) are handled by the package and environment manager Conda. You don't have to do anything except installing `miniconda`:

```
wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

#### How do I install DartQC?

Clone this repository and setup the virtual environment that contains Python dependencies and CD-HIT from BioConda:

```
cd ~
git clone https://github.com/esteinig/dartQC
conda env create -f ~/dartQC/env/dartqc.yaml
```

#### How do I run the script?

Before you run the script, you need to activate the virtual environment. It can be deactivated after you are done with the command `source deactivate`.

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

You can then sym-link the script to `bin` and use it from anywhere on your system without typing the full path to the script. For example if `~/bin` is in `PATH`:

```
ln -s ~/dartQC/dartqc.py ~/bin

source activate dartqc
dartqc.py --help
```
