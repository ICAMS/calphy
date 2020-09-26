import os
import numpy as np
import shutil
import argparse as ap
import subprocess
import yaml

from pytint.input import read_yamlfile
from pytint.liquid import Liquid
from pytint.solid import Solid

def submit_job(options, simfolder, scriptpath):
    """
    Submit a job
    """
    #create run command
    cores = options["queue"]["cores"]
    if cores > 1:
        cmd = ["mpirun", "--oversubscribe", "-n", str(cores), options["main"]["lammps_exec"], 
        "-in", scriptpath, "-screen", "lammscreen.log"]
    else:
        cmd = [options["main"]["lammps_exec"], "-in", scriptpath, "-screen", "lammscreen.log"]

    #now launch subprocess
    process = subprocess.Popen(
        cmd,
        stdout=None,
        stderr=subprocess.PIPE,
        stdin=None,
        cwd = simfolder,
    )

    #now wait until process is done
    out, err = process.communicate()
    err = err.decode("utf-8")
    if len(err) > 0:
        raise RuntimeError(err)    


def main():
    arg = ap.ArgumentParser()
    
    #argument name of input file
    arg.add_argument("-i", "--input", required=True, type=str,
    help="name of the input file")

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


    args = vars(arg.parse_args())
    options = read_yamlfile(args["input"])

    #we have some housekeepin to do now    
    l = args["lattice"]
    ml = args["mainlattice"]
    t = int(args["temperature"])
    p = "%.02f"%args["pressure"]
    c = "%.02f"%args["concentration"]
    thigh = max(options["main"]["temperature"])*1.5

    identistring = "-".join([l, str(t), p, c])
    simfolder = os.path.join(os.getcwd(), identistring)
    
    #if folder exists, delete it -> then create
    if os.path.exists(simfolder):
        shutil.rmtree(simfolder)
    os.mkdir(simfolder)

    #time to set up the job file
    if args["lattice"] == "liquid":
        job = Liquid(t = args["temperature"], p = args["pressure"],
                    l = ml, apc = args["atomspercell"],
                    alat = args["latticeconstant"],
                    c = args["concentration"], options=options,
                    simfolder = simfolder, thigh = thigh)
    else:
        job = Solid(t = args["temperature"], p = args["pressure"],
                    l = args["lattice"], apc = args["atomspercell"],
                    alat = args["latticeconstant"],
                    c = args["concentration"], options=options,
                    simfolder = simfolder)

    #initial routine
    scriptpath = os.path.join(simfolder, "mdavg.in")
    #switch folder
    os.chdir(simfolder)
    job.write_average_script(scriptpath)
    submit_job(options, simfolder, scriptpath)
    
    #gather data and submit second routine
    job.gather_average_data()
    scriptpath = os.path.join(simfolder, "mdint.in")
    job.write_integrate_script(scriptpath)
    submit_job(options, simfolder, scriptpath)

    #integrate and write report
    job.thermodynamic_integration()
    job.submit_report()