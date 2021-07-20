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
import time
import yaml
import warnings

from calphy.input import read_yamlfile, create_identifier
import calphy.scheduler as pq
import argparse as ap


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
    options = read_yamlfile(inputfile)
    print("Total number of %d calculations found" % len(options["calculations"]))

    for count, calc in enumerate(options["calculations"]):

        identistring = create_identifier(calc)
        scriptpath = os.path.join(os.getcwd(), ".".join([identistring, "sub"]))
        errfile = os.path.join(os.getcwd(), ".".join([identistring, "err"]))

        #the below part assigns the schedulers
        #now we have to write the submission scripts for the job
        #parse Queue and import module
        if options["queue"]["scheduler"] == "local":
            scheduler = pq.Local(options["queue"], cores=options["queue"]["cores"])
        elif options["queue"]["scheduler"] == "slurm":
            scheduler = pq.SLURM(options["queue"], cores=options["queue"]["cores"])
        elif options["queue"]["scheduler"] == "sge":
            scheduler = pq.SGE(options["queue"], cores=options["queue"]["cores"])
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
    arg.add_argument("-i", "--input", required=True, type=str,
    help="name of the input file")
    
    #parse args
    args = vars(arg.parse_args())

    #spawn job
    run_jobs(args["input"])
