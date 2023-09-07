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

from calphy.input import read_inputfile #, create_identifier
import calphy.scheduler as pq
import argparse as ap
from calphy import __version__ as version

def run_jobs(inputfile):
    """
    Spawn jobs which are submitted to cluster

    Parameters
    ----------
    options : dict
        dict containing input options
    
    Returns
    -------
    None
    """
    
    #the jobs are well set up in calculations dict now
    #Step 1 - loop over calcs dict
    #Step 2 - check input structure of calc and create lattice if needed
    #Step 3 - Submit job
    
    #read the input file
    calculations = read_inputfile(inputfile)
    print("Total number of %d calculations found" % len(calculations))

    if calculations[0].script_mode:
        for count, calc in enumerate(calculations):
            identistring = calc.create_identifier()
            scriptpath = os.path.join(os.getcwd(), ".".join([identistring, "sub"]))
            errfile = os.path.join(os.getcwd(), ".".join([identistring, "err"]))
            if calc.queue.scheduler == "local":
                scheduler = pq.Local(calc.queue.__dict__, cores=calc.queue.cores)
            elif calc.queue.scheduler == "slurm":
                scheduler = pq.SLURM(calc.queue.__dict__, cores=calc.queue.cores)
            elif calc.queue.scheduler == "sge":
                scheduler = pq.SGE(calc.queue.__dict__, cores=calc.queue.cores)
            else:
                raise ValueError("Unknown scheduler")

            #add run commands
            command = f'calphy_run_averaging -i {inputfile} -k {count}'
            scheduler.queueoptions['commands'].append(command)

            command = f'cp input.yaml {os.path.join(os.getcwd(), identistring)}'
            scheduler.queueoptions['commands'].append(command)

            command = f'cd {os.path.join(os.getcwd(), identistring)}'
            scheduler.queueoptions['commands'].append(command)
            if calc.queue.cores > 1:
                #here turn on mpi
                command = f'{calc.mpi_executable} -np {calc.queue.cores} {calc.lammps_executable} -in averaging.lmp'
            else:
                command = f'{calc.lammps_executable} -in averaging.lmp'
            scheduler.queueoptions['commands'].append(command)
            scheduler.queueoptions['commands'].append('cd ..')

            command = f'calphy_process_averaging -i {inputfile} -k {count}'
            scheduler.queueoptions['commands'].append(command)

            command = f'calphy_run_integration -i {inputfile} -k {count}'
            scheduler.queueoptions['commands'].append(command)

            command = f'cd {os.path.join(os.getcwd(), identistring)}'
            scheduler.queueoptions['commands'].append(command)
            if calc.queue.cores > 1:
                #here turn on mpi
                command = f'{calc.mpi_executable} -np {calc.queue.cores} {calc.lammps_executable} -in integration.lmp'
            else:
                command = f'{calc.lammps_executable} -in integration.lmp'
            scheduler.queueoptions['commands'].append(command)
            scheduler.queueoptions['commands'].append('cd ..')

            command = f'calphy_process_integration -i {inputfile} -k {count}'
            scheduler.maincommand = command
            scheduler.write_script(scriptpath)
            #_ = scheduler.submit()



    else:
        for count, calc in enumerate(calculations):

            identistring = calc.create_identifier()
            scriptpath = os.path.join(os.getcwd(), ".".join([identistring, "sub"]))
            errfile = os.path.join(os.getcwd(), ".".join([identistring, "err"]))

            #the below part assigns the schedulers
            #now we have to write the submission scripts for the job
            #parse Queue and import module
            if calc.queue.scheduler == "local":
                scheduler = pq.Local(calc.queue.__dict__, cores=calc.queue.cores)
            elif calc.queue.scheduler == "slurm":
                scheduler = pq.SLURM(calc.queue.__dict__, cores=calc.queue.cores)
            elif calc.queue.scheduler == "sge":
                scheduler = pq.SGE(calc.queue.__dict__, cores=calc.queue.cores)
            else:
                raise ValueError("Unknown scheduler")

            #for lattice just provide the number of position
            scheduler.maincommand = "calphy_kernel -i %s -k %d"%(inputfile, 
                count)
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
    
    #argument name of input file
    arg.add_argument("-i", "--input", required=False, type=str,
    help="name of the input file")

    arg.add_argument("-v", "--version", action='store_true',
    help="name of the input file")
    
    #parse args
    args = vars(arg.parse_args())

    if args["version"]:
        print(version)
    else:
        #spawn job
        if args["input"]:
            run_jobs(args["input"])
