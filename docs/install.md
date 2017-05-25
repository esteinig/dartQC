## Unix system with Python 3

Requires Python Packacke Index (PyPI) installer `pip`, bundled with the latest versions for Python 3.

```
pip install dartqc
dartqc install
```

## JCU's HPC (Zodiac)

JCU's HPC needs a manual installation of `conda`, since the system-wide default Python 2.6 does not come with an installer for PyPI.

Please note that activating the conda environment requires a Bash shell instead of the default shell on HPC (TCSH). You can check what shell you are using and activate Bash by:

```
echo $0      # if -tcsh
/bin/bash    # enter bash shell
 ```
 
#### How do I install the dependencies?

Dependencies (Python, CD-HIT) are handled by the package and environment manager Conda. You don't have to do anything except installing `miniconda`. In the final step of its installation, answer `yes` and let the installer append `miniconda` to you `PATH`.

```
wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh --prefix=$HOME/miniconda
source ~/.bashrc
```

#### How do I install DartQC?

Simply with `pip` which is now available through `conda` followed by the environment installer in `dartqc`:

```
pip install dartqc
dartqc install
```
