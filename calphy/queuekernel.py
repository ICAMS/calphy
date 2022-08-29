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
import shutil
import argparse as ap
import subprocess
import yaml
import time
import datetime

from calphy.input import read_inputfile
from calphy.liquid import Liquid
from calphy.solid import Solid
from calphy.alchemy import Alchemy
from calphy.routines import MeltingTemp, routine_fe, routine_ts, routine_only_ts, routine_pscale, routine_tscale, routine_alchemy, routine_composition_scaling


def setup_calculation(calc):
    """
    Set up a calculation

    Parameters
    ----------
    options: dict
        options object

    kernel: int
        index of the calculation to be run

    Returns
    -------
    job: Phase class
        job class
    """
    #now we need to modify the routines
    if calc.mode == "melting_temperature":
        simfolder = None
        job = MeltingTemp(calculation=calc, simfolder=simfolder)
    elif calc.mode == "alchemy" or calc.mode == "composition_scaling":
        simfolder = calc.create_folders()
        job = Alchemy(calculation=calc, simfolder=simfolder)
    else:
        simfolder = calc.create_folders()
        if calc.reference_phase == "liquid":
            job = Liquid(calculation=calc, simfolder=simfolder)
        else:
            job = Solid(calculation=calc, simfolder=simfolder)

    return job

def run_calculation(job):
    """
    Run calphy calculation

    Parameters
    ----------
    job: Phase class

    Returns
    -------
    job : Phase class
    """
    if job.calc.mode == "fe":
        job = routine_fe(job)
    elif job.calc.mode == "ts":
        job = routine_ts(job)
    elif job.calc.mode == "mts":
        job = routine_only_ts(job)
    elif job.calc.mode == "alchemy":
        job = routine_alchemy(job)
    elif job.calc.mode == "melting_temperature":
        job.calculate_tm()
    elif job.calc.mode == "tscale":
        job = routine_tscale(job)
    elif job.calc.mode == "pscale":
        job = routine_pscale(job)
    elif job.calc.mode == "composition_scaling":
        job = routine_composition_scaling(job)
    else:
        raise ValueError("Mode should be either fe/ts/mts/alchemy/melting_temperature/tscale/pscale/composition_scaling")
    return job

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
    calculations = read_inputfile(args["input"])

    calc = calculations[kernel]
    
    #format and parse the arguments
    identistring = calc.create_identifier()
    simfolder = os.path.join(os.getcwd(), identistring)

    #if folder exists, delete it -> then create
    if os.path.exists(simfolder):
        shutil.rmtree(simfolder)
    os.mkdir(simfolder)

    if calc.mode == "melting_temperature":
        os.rmdir(simfolder)
        simfolder = None
        job = MeltingTemp(calculation=calc, simfolder=simfolder)
    elif calc.mode == "alchemy" or calc.mode == "composition_scaling":
        job = Alchemy(calculation=calc, simfolder=simfolder)
        os.chdir(simfolder)
    else:
        if calc.reference_phase == "liquid":
            job = Liquid(calculation=calc, simfolder=simfolder)
        else:
            job = Solid(calculation=calc, simfolder=simfolder)
        os.chdir(simfolder)

    if job.calc.mode == "fe":
        _ = routine_fe(job)
    elif job.calc.mode == "ts":
        _ = routine_ts(job)
    elif job.calc.mode == "mts":
        _ = routine_only_ts(job)
    elif job.calc.mode == "alchemy":
        _ = routine_alchemy(job)
    elif job.calc.mode == "melting_temperature":
        job.calculate_tm()
    elif job.calc.mode == "tscale":
        _ = routine_tscale(job)
    elif job.calc.mode == "pscale":
        _ = routine_pscale(job)
    elif job.calc.mode == "composition_scaling":
        _ = routine_composition_scaling(job)
    else:
        raise ValueError("Mode should be either fe/ts/mts/alchemy/melting_temperature/tscale/pscale/composition_scaling")
