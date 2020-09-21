import os
import numpy as np

from pytint.input import read_yamlfile
import pytint.queue as pq

def spawn_jobs(inputfile):
    options = read_yamlfile(inputfile)

    #check the reqd temps
    if not len(options["main"]["tm"]) == 2:
        raise valueError("Length of input temperature should be 2")

    if not len(options["main"]["tsims"]) > 1:
        raise valueError("Required sims should be atleast 2")

    temparray = np.linspace(options["main"]["tm"][0], options["main"]["tm"][1], options["main"]["tsims"], endpoint=True)

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


    #now we have to create a list of commands for the scheduler
    #which is to run queuekernel - which will then write everything
    