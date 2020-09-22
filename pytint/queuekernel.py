import os
import numpy as np
import shutil

import pytint.average_scripts as pavg
from pytint.input import read_yamlfile
import argparse as ap

def main():
    arg = ap.ArgumentParser()
    
    #argument name of input file
    arg.add_argument("-i", "--input", required=True, type=str,
    help="name of the input file")

    arg.add_argument("-t", "--temperature", required=True, type=float,
    help="temperature for the simulation")

    arg.add_argument("-c", "--concentration", required=False, type=float,
    default=0.00, help="concentration for the simulation")

    arg.add_argument("-s", "--solid", required=True, type=bool,
    help="is the system in solid state?")

    args = vars(arg.parse_args())
    options = read_yamlfile(args["input"])
    #now write the inputfile for first part
    if args["solid"]:
        basename = os.path.join(os.getcwd(), "solid")

    simfolder = ".".join([basename, str(args["temperature"]), str(args["concentration"])])
    
    if os.path.exists(simfolder):
        shutil.rmtree(simfolder)

    os.mkdir(simfolder)

    #we need to copy the potential file?
    #no no pair style should provide complete path
    scriptpath = os.path.join(simfolder, "mdavg.in")

    if args["solid"]:
        pavg.write_script_solid(scriptpath, args["temperature"],
            options["md"]["pressure"], options)
    else:
        thigh = options["main"]["tm"][1]*1.5
        pavg.write_script_liquid(scriptpath, args["temperature"], args["temperature"],
            options["md"]["pressure"], options)


      