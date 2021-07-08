"""
Main kernel methods for pytint

There are two modes of operation:

- The basic one in which the temperature will be split on regular intervals
  and an integration calculation is done at each point.

- Reversible scaling mode in which one integration calculation will be done
  and one reversible scaling calculation will be done.

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
