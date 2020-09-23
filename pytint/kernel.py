import os
import numpy as np
import time
import yaml

from pytint.input import read_yamlfile
import pytint.queue as pq
import argparse as ap

def spawn_jobs(inputfile, monitor=False):
    
    options = read_yamlfile(inputfile)
    #gather job arrays
    temp = options["main"]["temperature"]
    press = options["main"]["pressure"]
    lattice = options["main"]["lattice"]
    conc = options["main"]["concentration"]

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
    errfiles = []
    
    nocalcs = len(temp)*len(press)*len(lattice)*len(conc)
    print("Total number of %d calculations registered" % nocalcs)

    #main looping starts
    for t in temp:
        for p in press:
            for count, l in enumerate(lattice):
                for c in conc:
                    identifier = "-".join([str(l), str(int(temp)), 
                                        str(int(press)), "%.02f"%c])
                    scriptpath = os.path.join(os.getcwd(), ".".join([identifier, "sub"]))
                    errfile = os.path.join(os.getcwd(), ".".join([identifier, "err"]))
                    errfiles.append(errfile)
                    #for lattice just provide the number of position
                    scheduler.maincommand = "tint_kernel -i %s -t %f -p %f -l %d -c %f"%(inputfile, t, p, count, c)
                    scheduler.write_script(scriptpath)
                    _ = scheduler.submit()


    if monitor:
        raise NotImplementedError("feature not implemented")

def integrate():
    #WARNING: This method is not updated
    #grab the values
    sfe = []
    for rep in sreports:
        with open(rep, 'r') as fout:
            data = yaml.load(fout, Loader=yaml.FullLoader)
        sfe.append(data["fe"])

    lfe = []
    for rep in lreports:
        with open(rep, 'r') as fout:
            data = yaml.load(fout, Loader=yaml.FullLoader)
        lfe.append(data["fe"])

    #now we have to do linear fits
    #WARNING: check quality of fit
    sfit = np.polyfit(temparray, sfe, 1)
    lfit = np.polyfit(temparray, lfe, 1)

    ntemp = np.arange(options["main"]["tm"][0], options["main"]["tm"][1]+1, 1)
    diff = np.polyval(lfit, ntemp) - np.polyval(sfit, ntemp)
    minarg = np.argsort(np.abs(diff))[0]

    res = np.column_stack((temparray, sfe, lfe))
    resfile = os.path.join(os.getcwd(), "results.dat")
    np.savetxt(resfile, res, header="temp solid liquid")

    if not (sfe[0]-lfe[0])*(sfe[-1]-lfe[-1]) < 0:
        raise RuntimeError("Melting temp not within range, or calculations not converged")
        
    print("Calculated Tm = %f with Dg = %f"%(ntemp[minarg], diff[minarg]))
    print("Results saved in results.dat")

def main():
    arg = ap.ArgumentParser()
    
    #argument name of input file
    arg.add_argument("-i", "--input", required=True, type=str,
    help="name of the input file")
    args = vars(arg.parse_args())
    
    spawn_jobs(args["input"])    