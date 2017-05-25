## Task: install

Once you have installed `dartqc` through `pip`, run `dartqc install` to check for a `conda` environment manager. If not present, the install task will download `miniconda` and install to `$HOME/miniconda`. The install task then installs the virtual environment `dartqc.yaml` for DartQC. This environment contains the dependencies for the program and needs to be activated before using DartQC:

`source activate dartqc`

The environment includes:
* CD-HIT (BioConda)
* BioPython
* SciPy
* NumPy
* Pandas
