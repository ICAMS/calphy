"""
calphy: a Python library and command line interface for automated free
energy calculations.

Copyright 2021  (c) Sarath Menon^1, Yury Lysogorskiy^1, Ralf Drautz^1
^1: Ruhr-University Bochum, Bochum, Germany

More information about the program can be found in:
Menon, Sarath, Yury Lysogorskiy, Jutta Rogal, and Ralf Drautz. 
“Automated Free Energy Calculation from Atomistic Simulations.” 
ArXiv:2107.08980 [Cond-Mat], July 19, 2021. 
http://arxiv.org/abs/2107.08980.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

See the LICENSE file.

For more information contact:
sarath.menon@ruhr-uni-bochum.de
"""

import os
import numpy as np
import shutil
import argparse as ap
import subprocess
import yaml

from calphy.input import read_yamlfile, create_identifier
from calphy.liquid import Liquid
from calphy.solid import Solid
from calphy.alchemy import Alchemy


def routine_fe(job):
    """
    Perform an FE calculation routine
    """
    job.run_averaging()
    #now run integration loops
    for i in range(job.nsims):
        job.run_integration(iteration=(i+1))

    job.thermodynamic_integration()
    job.submit_report()

def routine_ts(job):
    """
    Perform ts routine
    """
    routine_fe(job)

    #now do rev scale steps
    for i in range(job.nsims):
        job.reversible_scaling(iteration=(i+1))
    
    job.integrate_reversible_scaling(scale_energy=True)


def routine_only_ts(job):
    """
    Perform sweep without free energy calculation
    """
    job.run_averaging()
    for i in range(job.nsims):
        job.reversible_scaling(iteration=(i+1))


def routine_alchemy(job):
    """
    Perform an FE calculation routine
    """
    job.run_averaging()
    #now run integration loops
    for i in range(job.nsims):
        job.run_integration(iteration=(i+1))

    job.thermodynamic_integration()
    job.submit_report()    

def main():
    arg = ap.ArgumentParser()
    
    #argument name of input file
    arg.add_argument("-i", "--input", required=True, type=str,
    help="name of the input file")

    arg.add_argument("-k", "--kernel", required=True, type=int, 
    help="kernel number of the calculation to be run.")


    #parse input
    #parse arguments
    args = vars(arg.parse_args())
    kernel = args["kernel"]
    options = read_yamlfile(args["input"])

    calc = options["calculations"][kernel]
    
    #format and parse the arguments
    #thigh is for now hardcoded    
    identistring = create_identifier(calc)
    simfolder = os.path.join(os.getcwd(), identistring)

    #if folder exists, delete it -> then create
    if os.path.exists(simfolder):
        shutil.rmtree(simfolder)
    os.mkdir(simfolder)

    #now we need to modify the routines
    if calc["mode"] == "alchemy":
        job = Alchemy(options=options, kernel=kernel, simfolder=simfolder)
    else:
        if calc["state"] == "liquid":
            job = Liquid(options=options, kernel=kernel, simfolder=simfolder)
        else:
            job = Solid(options=options, kernel=kernel, simfolder=simfolder)

    #integration routine
    os.chdir(simfolder)

    if calc["mode"] == "fe":
        routine_fe(job)
    elif calc["mode"] == "ts":
        routine_ts(job)
    elif calc["mode"] == "mts":
        routine_only_ts(job)
    elif calc["mode"] == "alchemy":
        routine_alchemy(job)
    else:
        raise ValueError("Mode should be either fe/ts/mts/alchemy")
