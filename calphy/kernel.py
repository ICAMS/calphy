"""
calphy: a Python library and command line interface for automated free
energy calculations.

Copyright 2021  (c) Sarath Menon^1, Yury Lysogorskiy^2, Ralf Drautz^2
^1: Max Planck Institut für Eisenforschung, Dusseldorf, Germany
^2: Ruhr-University Bochum, Bochum, Germany

calphy is published and distributed under the Academic Software License v1.0 (ASL).
calphy is distributed in the hope that it will be useful for non-commercial academic research,
but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
calphy API is published and distributed under the BSD 3-Clause "New" or "Revised" License
See the LICENSE FILE for more details.

More information about the program can be found in:
Menon, Sarath, Yury Lysogorskiy, Jutta Rogal, and Ralf Drautz.
“Automated Free Energy Calculation from Atomistic Simulations.” Physical Review Materials 5(10), 2021
DOI: 10.1103/PhysRevMaterials.5.103801

For more information contact:
sarath.menon@ruhr-uni-bochum.de/yury.lysogorskiy@icams.rub.de
"""

import os
import numpy as np
import time
import yaml
import warnings

from calphy.input import read_inputfile  # , create_identifier
import calphy.scheduler as pq
import argparse as ap
from calphy import __version__ as version


def _getattr_safe(obj, attr):
    """Safely get attribute whether obj is Pydantic model or dict."""
    if isinstance(obj, dict):
        return obj.get(attr)
    return getattr(obj, attr)


def _to_dict(obj):
    """Convert to dict whether obj is Pydantic model or dict."""
    if isinstance(obj, dict):
        return obj
    return obj.__dict__


def run_jobs(inputfile, validate=False):
    """
    Spawn jobs which are submitted to cluster

    Parameters
    ----------
    options : dict
        dict containing input options
    validate : bool
        if True, perform full validation during parsing (slower).
        if False, uses fast parsing with dict-like access.

    Returns
    -------
    None
    """

    # the jobs are well set up in calculations dict now
    # Step 1 - loop over calcs dict
    # Step 2 - check input structure of calc and create lattice if needed
    # Step 3 - Submit job

    # Fast parsing by default - use helper functions for attribute access
    calculations = read_inputfile(inputfile, validate=validate)
    print("Total number of %d calculations found" % len(calculations))

    for count, calc in enumerate(calculations):
        script_mode = _getattr_safe(calc, "script_mode")
        queue = _getattr_safe(calc, "queue")

        if script_mode:
            identistring = calc.create_identifier()
            scriptpath = os.path.join(os.getcwd(), ".".join([identistring, "sub"]))
            errfile = os.path.join(os.getcwd(), ".".join([identistring, "err"]))

            scheduler_type = _getattr_safe(queue, "scheduler")
            queue_dict = _to_dict(queue)
            cores = _getattr_safe(queue, "cores")

            if scheduler_type == "local":
                scheduler = pq.Local(queue_dict, cores=cores)
            elif scheduler_type == "slurm":
                scheduler = pq.SLURM(queue_dict, cores=cores)
            elif scheduler_type == "sge":
                scheduler = pq.SGE(queue_dict, cores=cores)
            else:
                raise ValueError("Unknown scheduler")

            # add run commands
            command = f"calphy_run_averaging -i {inputfile} -k {count}"
            scheduler.queueoptions["commands"].append(command)

            command = f"cp input.yaml {os.path.join(os.getcwd(), identistring)}"
            scheduler.queueoptions["commands"].append(command)

            command = f"cd {os.path.join(os.getcwd(), identistring)}"
            scheduler.queueoptions["commands"].append(command)

            mpi_exec = _getattr_safe(calc, "mpi_executable")
            lammps_exec = _getattr_safe(calc, "lammps_executable")

            if cores > 1:
                # here turn on mpi
                command = f"{mpi_exec} -np {cores} {lammps_exec} -in averaging.lmp"
            else:
                command = f"{lammps_exec} -in averaging.lmp"
            scheduler.queueoptions["commands"].append(command)
            scheduler.queueoptions["commands"].append("cd ..")

            command = f"calphy_process_averaging -i {inputfile} -k {count}"
            scheduler.queueoptions["commands"].append(command)

            command = f"calphy_run_integration -i {inputfile} -k {count}"
            scheduler.queueoptions["commands"].append(command)

            command = f"cd {os.path.join(os.getcwd(), identistring)}"
            scheduler.queueoptions["commands"].append(command)
            if cores > 1:
                # here turn on mpi
                command = f"{mpi_exec} -np {cores} {lammps_exec} -in integration.lmp"
            else:
                command = f"{lammps_exec} -in integration.lmp"
            scheduler.queueoptions["commands"].append(command)
            scheduler.queueoptions["commands"].append("cd ..")

            command = f"calphy_process_integration -i {inputfile} -k {count}"
            scheduler.maincommand = command
            scheduler.write_script(scriptpath)
            # _ = scheduler.submit()

        else:

            identistring = calc.create_identifier()
            scriptpath = os.path.join(os.getcwd(), ".".join([identistring, "sub"]))
            errfile = os.path.join(os.getcwd(), ".".join([identistring, "err"]))

            # the below part assigns the schedulers
            # now we have to write the submission scripts for the job
            # parse Queue and import module

            scheduler_type = _getattr_safe(queue, "scheduler")
            queue_dict = _to_dict(queue)
            cores = _getattr_safe(queue, "cores")

            if scheduler_type == "local":
                scheduler = pq.Local(queue_dict, cores=cores)
            elif scheduler_type == "slurm":
                scheduler = pq.SLURM(queue_dict, cores=cores)
            elif scheduler_type == "sge":
                scheduler = pq.SGE(queue_dict, cores=cores)
            else:
                raise ValueError("Unknown scheduler")

            scheduler.queueoptions["jobname"] = "".join(
                e for e in identistring if e.isalnum()
            )

            # for lattice just provide the number of position
            scheduler.maincommand = "calphy_kernel -i %s -k %d" % (inputfile, count)
            scheduler.write_script(scriptpath)
            _ = scheduler.submit()


def main():
    """
    Main method to parse arguments and run jobs

    Paramaters
    ----------
    None

    Returns
    -------
    None
    """
    arg = ap.ArgumentParser()

    # argument name of input file
    arg.add_argument(
        "-i", "--input", required=False, type=str, help="name of the input file"
    )

    arg.add_argument(
        "-v", "--version", action="store_true", help="name of the input file"
    )

    arg.add_argument(
        "--validate",
        action="store_true",
        default=False,
        help="perform full validation during input parsing (slower)",
    )

    # parse args
    args = vars(arg.parse_args())

    if args["version"]:
        print(version)
    else:
        # spawn job
        if args["input"]:
            run_jobs(args["input"], validate=args["validate"])
