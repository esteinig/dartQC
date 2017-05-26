## Task: pbs

```
dartqc pbs [--help] [--email] [--walltime]

Arguments:

--email, -e         email address [user@hpc.jcu.edu.au]
--walltime, -w      expected run time in hh:mm:ss [02:00:00]
--memory, -m        memory for job job
--processors, -p    number of processors for job
--pypi              add call to activate environment if installed via pypi
```

The only purpose of this task is to create a job submission script for PBS/Torque scheduling systems. It is good practice and behaviour  to run the processing and filtering on compute nodes instead of the login node. The command requires an e-mail address for user identification and you can set the expected walltime, default is `02:00:00`. This generates `dartqc.pbs`:

`dartqc pbs --email esteinig@jcu.edu.au --walltime 02:00:00`

You can then run any command like this:

`qsub -v -command="dartqc filter --calls call_file.csv --call_scheme call_scheme.json --maf 0.02" ./dartqc.pbs`


If `--pypi` is switched on, the script activates the `dartqc` environment before execution on the node. You need to have run `dartqc install` for this to work. 

The script can also be found in this repository at `pbs/dartqc.pbs`.
