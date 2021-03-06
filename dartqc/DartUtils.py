import os
import time
import argparse
import textwrap

from subprocess import call, check_output, CalledProcessError
import logging


class Installer:
    def __init__(self, miniconda=True):

        self.miniconda = miniconda
        self.miniconda_url = "https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh"
        self.miniconda_install = "Miniconda3-latest-Linux-x86_64.sh"

        self.base_path = os.path.dirname(os.path.realpath(__file__))

        self.env = os.path.join(self.base_path, "env", "dartqc.yaml")

        if not self._check() and self.miniconda:
            self._install_miniconda()

        self._install_env()

    @staticmethod
    def _check():

        try:
            check_output(["conda", "--version"])
            stamp("Conda manager detected, use: conda --version")
            return True
        except CalledProcessError:
            stamp("Could not detect conda manager, please see README.")
            return False

    def _install_miniconda(self):

        stamp("Could not detect Conda. Installing Miniconda for Python 3.")

        try:
            stamp("Downloading and installing...")
            with open("install.log", "w") as err_file:
                check_output(["wget", self.miniconda_url, "&&", "bash", self.miniconda_install,
                              "-b", "-p", "$HOME/miniconda"], stderr=err_file)
                stamp("Success. Removing installer for Miniconda...")
                os.remove(os.path.join(os.getcwd(), self.miniconda_install))
                if os.path.exists("install.log"):
                    os.remove("install.log")
        except CalledProcessError:
            stamp("Could not install Miniconda, please see install.log and README.")

        try:
            stamp("Adding $HOME/miniconda to PATH.")
            with open("install.log", "w") as err_file:
                call("""echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> $HOME/.bashrc && source .bashrc""",
                     stderr=err_file)
                stamp("Done. Testing...")
                self._check()
                if os.path.exists("install.log"):
                    os.remove("install.log")
        except CalledProcessError:
            stamp("Could not install Miniconda, please see install.log and README.")

    def _install_env(self):

        try:
            with open("install.log", "w") as err_file:
                stamp("Installing environment for DartQC...")
                check_output(["conda", "env", "create", "--name", "dartqc", "--file", self.env], stderr=err_file)
            stamp("Installed environment, activate with: source activate dartqc")
            if os.path.exists("install.log"):
                os.remove("install.log")
        except CalledProcessError:
            stamp("Could not install environment, please see install.log and README.")


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


class CommandLine:
    def __init__(self):
        parser = argparse.ArgumentParser()

        parser.add_argument("--version", "-v", dest="version", action="store_true", help="print version and exit")

        parser.add_argument("--output_path", "-o", type=lambda p: os.path.abspath(p),
                            default=os.getcwd(), required=False, dest="out_path", help="output path")

        parser.add_argument("--project", "-p", type=str, default="dartqc", required=False,
                            dest="project", help="project name")

        parser.add_argument("--populations", "--pop", type=lambda p: os.path.abspath(p), default=None, required=False,
                            dest="pop_file", help="CSV file with header columns ID and Population")

        parser.add_argument("--log_path", type=lambda p: os.path.abspath(p), required=False, dest="log_path",
                            help="log path")

        subparsers = parser.add_subparsers(help='Command-line interface for DartQC')

        install_parser = subparsers.add_parser("install")

        install_parser.set_defaults(subparser='install')

        pbs_parser = subparsers.add_parser("pbs")

        pbs_parser.add_argument("--email", "-e", type=str, default="user@hpc.jcu.edu.au", required=False,
                                dest="email", help="user email for pbs")

        pbs_parser.add_argument("--walltime", "-w", type=str, default="02:00:00", required=False,
                                dest="walltime", help="expected walltime for pbs")
        pbs_parser.add_argument("--memory", "-m", type=str, default="1gb", required=False,
                                dest="memory", help="memory for pbs")
        pbs_parser.add_argument("--cpus", "-c", type=str, default=1, required=False,
                                dest="cpus", help="cpus for pbs")
        pbs_parser.add_argument("--pypi", "-p", action="store_true", required=False,
                                dest="pypi", help="package installed with pypi, load dartqc env before command")

        pbs_parser.set_defaults(subparser='pbs')

        prepare_parser = subparsers.add_parser("prepare")

        prepare_parser.add_argument("--file", "-f", type=lambda p: os.path.abspath(p), required=True,
                                    dest="file", help="path to input file")
        prepare_parser.add_argument("--name", "-n", type=str, default=None, required=False,
                                    dest="output_name", help="name of output scheme file")
        prepare_parser.add_argument("--sheet", "-s", type=str, default=None, required=False,
                                    dest="sheet", help="if file is excel: name of sheet that contains double row data")

        prepare_parser.set_defaults(subparser='prepare')

        validate_parser = subparsers.add_parser("validate")

        validate_parser.add_argument("--raw", "-r", type=lambda p: os.path.abspath(p), required=True,
                                     dest="raw_file", help="path to raw read file")

        validate_parser.add_argument("--raw_scheme", default="raw_scheme.json",
                                     type=lambda p: os.path.abspath(p), required=True,
                                     dest="raw_scheme", help="path to raw scheme json file")

        validate_parser.add_argument("--calls", "-c", default="calls.csv", type=lambda p: os.path.abspath(p),
                                     required=True, dest="call_file", help="path to called read file")

        validate_parser.add_argument("--call_scheme", default="call_scheme.json",
                                     type=lambda p: os.path.abspath(p), required=True,
                                     dest="call_scheme", help="path to call scheme json file")

        validate_parser.add_argument("--id_list", "-i", default="id_list.csv", type=lambda p: os.path.abspath(p),
                                     required=True, dest="id_list",
                                     help="path to CSV file with list of official clone IDs that should be used (eg. to fix ID's that Dart outputs wrong)")

        validate_parser.add_argument("--cdhit_path", type=lambda p: os.path.abspath(p), required=False,
                                     dest="cdhit_path", default="cd-hit-est",
                                     help="Path to the cdhit 2d executable (required if cd-hit-2d doesn't work on cmd line)")

        validate_parser.set_defaults(subparser='validate')

        process_parser = subparsers.add_parser("process")

        process_parser.add_argument("--raw", "-r", type=lambda p: os.path.abspath(p), required=True,
                                    dest="raw_file", help="path to raw read file")

        process_parser.add_argument("--raw_scheme", default="raw_scheme.json",
                                    type=lambda p: os.path.abspath(p), required=False,
                                    dest="raw_scheme", help="path to raw scheme json file")

        process_parser.add_argument("--read_sum", "--reads", default=[10],
                                    type=lambda s: [int(item) for item in s.split(',')], required=False,
                                    dest="raw_read_threshold",
                                    help="silence call if ref and snp allele raw read sum < threshold")

        process_parser.add_argument("--calls", "-c", default="calls.csv", type=lambda p: os.path.abspath(p),
                                    required=True, dest="call_file", help="path to called read file")

        process_parser.add_argument("--call_scheme", default="call_scheme.json",
                                    type=lambda p: os.path.abspath(p), required=False,
                                    dest="call_scheme", help="path to call scheme json file")

        process_parser.add_argument("--graph", "-g", default=False, type=bool, required=False,
                                    dest="graph", help="Create graphs")

        process_parser.set_defaults(subparser='process')

        filter_parser = subparsers.add_parser("filter")

        filter_parser.add_argument("--processed", "--pp", type=lambda p: os.path.abspath(p), required=False,
                                   dest="processed_path", default=None,
                                   help="input path to processed data files (project_data.json, project_attr.json)")

        filter_parser.add_argument("--calls", "-c", default="calls.csv", type=lambda p: os.path.abspath(p),
                                   required=False, dest="call_file", help="path to called read file")
        filter_parser.add_argument("--call_scheme", default="call_scheme.json",
                                   type=lambda p: os.path.abspath(p), required=False,
                                   dest="call_scheme", help="path to call scheme json file")

        # Input values must follow [<val>, <val>, ...] - missing values should be empty as in [,,,]
        filter_parser.add_argument("--maf", default=[], type=lambda s: [float(item.strip()) if len(item.strip()) > 0 else None for item in s[1:-1].split(',')],
                                   dest="maf", help="filter snps <= minor allele frequency")
        filter_parser.add_argument("--hwe", default=[], type=lambda s: [float(item.strip()) if len(item.strip()) > 0 else None for item in s[1:-1].split(',')],
                                   dest="hwe", help="filter snps <= p-value of hardy-weinberg test")
        filter_parser.add_argument("--call_rate", default=[], type=lambda s: [float(item.strip()) if len(item.strip()) > 0 else None for item in s[1:-1].split(',')],
                                   dest="call_rate", help="filter snps <= call rate of snp")
        filter_parser.add_argument("--rep", default=[], type=lambda s: [float(item.strip()) if len(item.strip()) > 0 else None for item in s[1:-1].split(',')],
                                   dest="rep", help="filter snps <= replication average of snp")

        filter_parser.add_argument("--mind", default=[], type=lambda s: [float(item.strip()) if len(item.strip()) > 0 else None for item in s[1:-1].split(',')],
                                   dest="mind", help="filter samples > missingness per sample")
        filter_parser.add_argument("--mono", default=None,
                                   dest="mono", help="filter samples monomorphic in <mono> populations ('all', int)")
        filter_parser.add_argument("--mono_comparison", default="==",
                                   dest="mono_comp", help="filter samples monomorphic in >=, <=, == populations ('==')")

        filter_parser.add_argument("--split_clones", default="", type=str, dest="split_clones",
                                   help="split clone ids on this character")

        filter_parser.add_argument("--duplicates", default=False, action="store_true",
                                   dest="remove_duplicates", help="remove snps with duplicate clone IDs")
        filter_parser.add_argument("--clusters", default=False, action="store_true",
                                   dest="remove_clusters", help="remove snps in identical sequence clusters")
        filter_parser.add_argument("--identity", default=0.95, type=float,
                                   dest="identity", help="remove snps in identical sequence clusters")
        filter_parser.add_argument("--cdhit_path", type=lambda p: os.path.abspath(p), required=False,
                                   dest="cdhit_path", default="cd-hit-est",
                                   help="Path to the cdhit executable (required if cd-hit-est doesn't work on cmd line)")

        filter_parser.add_argument("--graph", "-g", default=False, type=bool, required=False,
                                   dest="graph", help="Create graphs")

        filter_parser.add_argument("--raw", "-r", type=lambda p: os.path.abspath(p), required=False,
                                   dest="raw_file", help="path to raw read file (only needed if graphing)")

        filter_parser.add_argument("--raw_scheme", default="raw_scheme.json",
                                   type=lambda p: os.path.abspath(p), required=False,
                                   dest="raw_scheme", help="path to raw scheme json file (only needed if graphing)")

        filter_parser.set_defaults(subparser='filter')

        self.args = parser.parse_args()


def stamp(*args):
    message = str(time.strftime("[%H:%M:%S]")) + " " + " ".join([str(arg) for arg in args])
    print(message)
    logging.info(message)
