import os
import numpy as np
import time
import yaml

from pytint.input import read_yamlfile
import pytint.queue as pq
import argparse as ap

def spawn_jobs(inputfile):
    options = read_yamlfile(inputfile)

    #check the reqd temps
    if not len(options["main"]["tm"]) == 2:
        raise valueError("Length of input temperature should be 2")

    if not options["main"]["tsims"] > 1:
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
    reports = []
    errfiles = []
    sreports = []
    lreports = []
    
    for temp in temparray:
        for conc in [0.0]:
            #spawn jobs
            #clear jobs if they exist
            identistring = ".".join(["solid", str(temp), "%.02f"%conc])
            reportfile = os.path.join(os.getcwd(), ".".join([identistring, "yaml"]))
            if os.path.exists(reportfile):
                os.remove(reportfile)

            #now make a scriptfile
            scriptpath = os.path.join(os.getcwd(), ".".join([identistring, "sub"]))
            errfile = os.path.join(os.getcwd(), ".".join([identistring, "sub", "err"]))
            errfiles.append(errfile)
            scheduler.maincommand = "tint_kernel -i %s -t %f -c %f -s yes"%(inputfile, temp, conc)
            scheduler.write_script(scriptpath)
            _ = scheduler.submit()
            reports.append(reportfile)
            sreports.append(reportfile)

            identistring = ".".join(["liquid", str(temp), "%.02f"%conc])
            reportfile = os.path.join(os.getcwd(), ".".join([identistring, "yaml"]))
            if os.path.exists(reportfile):
                os.remove(reportfile)

            #now make a scriptfile
            scriptpath = os.path.join(os.getcwd(), ".".join([identistring, "sub"]))
            errfile = os.path.join(os.getcwd(), ".".join([identistring, "sub", "err"]))
            errfiles.append(errfile)
            scheduler.maincommand = "tint_kernel -i %s -t %f -c %f -s no"%(inputfile, temp, conc)
            scheduler.write_script(scriptpath)
            _ = scheduler.submit()
            reports.append(reportfile)
            lreports.append(reportfile)

    #array of jobs are created
    #now monitor jobs regularly
    
    errored = []
    messages = []

    while(True):
        done = 0
        for count, report in enumerate(reports):
            #print(report)
            if os.path.exists(report):
                done += 1
        if (done == len(reports)):
            break
        #print(done, len(reports))
        time.sleep(options["main"]["updatetime"])
        #check if some errors exist
        errfile = errfiles[count]
        if os.path.exists(errfile):
            #check if the file is empty
            if not (os.stat(errfile).st_size == 0):
                file = open(errfile, mode='r')
                contents = file.read()
                errored.append(count)
                messages.append(contents)
        if len(errored) > 0:
            for c, err in enumerate(errored):
                print(err, messages[c])
            raise RuntimeError("Jobs failed")


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