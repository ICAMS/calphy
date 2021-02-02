"""
Queue kernel module contains the workflow that is run on single thread
"""

import os
import numpy as np
import shutil
import argparse as ap
import subprocess
import yaml

from pytint.input import read_yamlfile, create_identifier
from pytint.liquid import Liquid
from pytint.solid import Solid
import pytint.lattice as pl

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
    
    #process lattice
    lattice = calc["lattice"].upper()
    if lattice in ["BCC", "FCC", "HCP", "DIA", "SC", "LQD"]:
        #process lattice
        alat, apc, l = pl.get_lattice(element, lattice)
    elif os.path.exists(calc["lattice"]):
        #its a file - do something
        l = "file"
    else:
        raise ValueError("Unknown lattice found. Allowed options are BCC, FCC, HCP, DIA, SC or LQD; or an input file.")


    #format and parse the arguments
    #thigh is for now hardcoded    
    identistring = create_identifier(calc)
    simfolder = os.path.join(os.getcwd(), identistring)

    #if folder exists, delete it -> then create
    if os.path.exists(simfolder):
        shutil.rmtree(simfolder)
    os.mkdir(simfolder)

    #time to set up the job
    #create a lattice object
    #just tweak for diamond!
    if l == "dia":
        l = "diamond"

    #now we need to modify the routines
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
    else:
        raise ValueError("Mode should be either fe or ts")
