## Unix system with Python 3.5+

With `conda` (see below) look for `cd-hit` in channel `bioconda` and `dartqc` in `esteinig`:

```
conda install -c bioconda -c esteinig dartqc
```

Python Package Index (PyPI). This is a cumbersome way to install the program since the environment requires CD-HIT from the BioConda channel. The preferred way is through `conda` as described above.

You can install with pip and use the native installer task [`install`](https://github.com/esteinig/dartQC/blob/master/docs/task.install.md) in the program to install `miniconda` and the virtual environment. The environment needs to be then activated before using `dartqc`, i.e. `source activate dartqc`. When you want to update `dartqc` to a more recent version as the program is developed, you need to call `pip install dartqc --upgrade` from within the virtual environment. 

```
pip install dartqc
dartqc install
source activate dartqc
```
 
#### How do I install the Conda and Dependencies?

Dependencies (Python, CD-HIT) are handled by the package and environment manager Conda. You don't have to do anything except installing `miniconda`. In the final step of its installation, answer `yes` and let the installer append `miniconda` to you `PATH`.

```
wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

## JCU's HPC (Zodiac)

Please note that installing via conda may require a Bash shell instead of the default shell on HPC (TCSH). You can check what shell you are using and activate Bash by:

```
echo $0      # if -tcsh
bash         # enter bash shell
 ```
