"""
Queue kernel module contains the workflow that is run on single thread
"""

import os
import numpy as np
import shutil
import argparse as ap
import subprocess
import yaml

from pytint.input import read_yamlfile
from pytint.liquid import Liquid
from pytint.solid import Solid

def main():
    arg = ap.ArgumentParser()
    
    #argument name of input file
    arg.add_argument("-i", "--input", required=True, type=str,
    help="name of the input file")

    arg.add_argument("-j", "--job", required=False, type=str, 
    choices=["integrate", "rs"],
    default="integrate", help="mode of simulation - integrate or rs")

    arg.add_argument("-t", "--temperature", required=True, type=float,
    help="temperature for the simulation")

    arg.add_argument("-p", "--pressure", required=True, type=float,
    help="pressure for the simulation")

    arg.add_argument("-l", "--lattice", required=True, type=str,
    help="lattice for the simulation")

    arg.add_argument("-apc", "--atomspercell", required=True, type=int,
    help="Number of atoms per cell")

    arg.add_argument("-a", "--latticeconstant", required=True, type=float,
    help="lattice constant for the simulation")

    arg.add_argument("-c", "--concentration", required=True, type=float,
    help="concentration for the simulation")

    arg.add_argument("-m", "--mainlattice", required=True, type=str,
    help="lammps lattice for the simulation")

    arg.add_argument("-fe", "--freenergy", required=False, type=float,
    help="Initial free energy for reversible scaling method")

    #parse arguments
    args = vars(arg.parse_args())
    options = read_yamlfile(args["input"])

    #we have some housekeepin to do now
    #format and parse the arguments
    #thigh is for now hardcoded    
    l       = args["lattice"]
    ml      = args["mainlattice"]
    t       = int(args["temperature"])
    p       = "%.02f"%args["pressure"]
    c       = "%.02f"%args["concentration"]
    thigh   = max(options["main"]["temperature"])*1.5

    #check the kind of job we have at hand
    if args["job"] == "integrate":
        integrate = True
        skey = "in"
    elif args["job"] == "rs":
        integrate = False
        skey = "rs"

    #create an string which should be unique for the job in hand
    #the job should have an extra argument to indicate job time        
    identistring = "-".join([skey, l, str(t), p, c])
    simfolder = os.path.join(os.getcwd(), identistring)

    #if folder exists, delete it -> then create
    if os.path.exists(simfolder):
        shutil.rmtree(simfolder)
    os.mkdir(simfolder)

    #time to set up the job
    #create a lattice object
    if args["lattice"] == "LQD":
        job = Liquid(t = args["temperature"], p = args["pressure"],
                    l = ml, apc = args["atomspercell"],
                    alat = args["latticeconstant"],
                    c = args["concentration"], options=options,
                    simfolder = simfolder, thigh = thigh)
    else:
        job = Solid(t = args["temperature"], p = args["pressure"],
                    l = ml, apc = args["atomspercell"],
                    alat = args["latticeconstant"],
                    c = args["concentration"], options=options,
                    simfolder = simfolder)

    #integration routine
    os.chdir(simfolder)

    if integrate:
        job.run_averaging()
        #now run integration loops
        for i in range(options["main"]["nsims"]):
            job.run_integration(iteration=(i+1))

        job.thermodynamic_integration()
        job.submit_report()
    
    #reversible scaling routine
    else:
        #the rs routine
        job.run_averaging()

        #now find fe for one temp
        for i in range(options["main"]["nsims"]):
            job.run_integration(iteration=(i+1))

        #now do rev scale steps
        for i in range(options["main"]["nsims"]):
            job.reversible_scaling(iteration=(i+1))
        #we do not integrate rev scaling!
        #do it manually from a jupyter notebook or so
        job.integrate_reversible_scaling()
