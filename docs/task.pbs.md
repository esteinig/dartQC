## Task: pbs

The only purpose of this task is to create a job submission script for PBS/Torque scheduling systems. It is good practice and behaviour  to run the processing and filtering on compute nodes instead of the login node. The command requires an e-mail address for user identification and you can set the epected walltime (default is 02:00:00). This generates `dartqc.pbs`:

`dartqc pbs --email esteinig@jcu.edu.au --walltime 02:00:00`

The script activates the environment before execution on the nodes and you can run any command like this:

`qsub -v -command="dartqc filter --calls call_file.csv --call_scheme call_scheme.json --maf 0.02" ./dartqc.pbs`

The script can be found in this repository at `pbs/dartqc.pbs`.
