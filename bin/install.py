# Installation script to help with dependency management
import argparse
import inspect
import os
from subprocess import call, check_output, CalledProcessError

import time

import sys

cdhit_config = os.path.join(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))), "cdhit.txt")


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--miniconda", "-m", required=False, type=bool,default=True,
                        dest="miniconda", help="Use miniconda package manager")

    parser.add_argument("--cdhit_path", type=lambda p: os.path.abspath(p), required=False,
                        dest="cdhit_path", default="cd-hit-est",
                        help="Path to the cdhit 2d executable (required if cd-hit-2d doesn't work on cmd line)")

    args = vars(parser.parse_args())

    Installer(args["miniconda"], args["cdhit_path"])


class Installer:
    def __init__(self, miniconda=True, cdhit_path=None):

        self.miniconda = miniconda
        self.miniconda_url = "https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh"
        self.miniconda_install = "Miniconda3-latest-Linux-x86_64.sh"

        self.base_path = os.path.dirname(os.path.realpath(__file__))

        self.env = os.path.join(self.base_path, "env", "dartqc.yaml")

        self.cdhit_path = cdhit_path

        if not self._check() and self.miniconda:
            self._install_miniconda()

        self._install_env()
        self._install_cdhit()

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

    def _install_cdhit(self):
        try:
            # Check that this cdhit path works!
            with open("install.log", "w") as err_file:
                check_output([self.cdhit_path, "-h"], stderr=err_file)
        except CalledProcessError:
            output = sys.exc_info()[1].output

            if str(output).find("CD-HIT version") != -1:
                with open(cdhit_config, "w") as config:
                    config.write(self.cdhit_path)
                    config.flush()
            else:
                stamp("Could not install cdhit path, please see install.log and README.")


def stamp(*args):
    message = str(time.strftime("[%H:%M:%S]")) + " " + " ".join([str(arg) for arg in args])
    print(message)


if __name__ == "__main__":
    main()
