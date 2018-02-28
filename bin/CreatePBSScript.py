import argparse
import sys

# Creates a PBS job installation script for scheduling PBS jobs such as on a HPC
import textwrap


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--email", "-e", type=str, default="user@hpc.jcu.edu.au", required=False,
                            dest="email", help="user email for pbs")

    parser.add_argument("--walltime", "-w", type=str, default="02:00:00", required=False,
                            dest="walltime", help="expected walltime for pbs")
    parser.add_argument("--memory", "-m", type=str, default="1gb", required=False,
                            dest="memory", help="memory for pbs")
    parser.add_argument("--cpus", "-c", type=str, default=1, required=False,
                            dest="cpus", help="cpus for pbs")
    parser.add_argument("--pypi", "-p", action="store_true", required=False,
                            dest="pypi", help="package installed with pypi, load dartqc env before command")

    args = vars(parser.parse_args())

    PBS(email=args["email"], walltime=args["walltime"], processors=args["cpus"], memory=args["memory"],
        pypi_install=args["pypi"])


class PBS:
    def __init__(self, walltime="02:00:00", memory="1gb", processors=1, email="user@hpc.jcu.edu.au",
                 pypi_install=False):

        self.walltime = walltime
        self.memory = memory
        self.processors = processors
        self.pypi_install = pypi_install

        self.email = email

        self._write_pbs()

    def _write_pbs(self):

        if self.pypi_install:
            pypi = "true"
        else:
            pypi = "false"

        pbs = textwrap.dedent("""
        #PBS -S /bin/sh
        #PBS -M {email}
        #PBS -m n
        #PBS -j oe
        #PBS -l ncpus={cpus}
        #PBS -l pmem={mem}
        #PBS -l walltime={walltime}

        pypi_install={pypi}

        if [ "$pypi_install" = true ] ; then
            source activate dartqc
        fi

        source ~/.bashrc
        cd $PBS_O_WORKDIR

        echo $PATH

        if [ ! "$command" ]
            then
                echo -e "Error: script must be supplied with a command"
                exit 1
            else
                echo "$command"
                eval "$command"
        fi

        exit 0

        """.format(walltime=self.walltime, email=self.email, cpus=self.processors, mem=self.memory, pypi=pypi))

        with open("dartqc.pbs", "w") as outfile:
            outfile.write(pbs)


if __name__ == "__main__":
    main()