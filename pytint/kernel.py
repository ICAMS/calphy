import os
import numpy as np
import time
import yaml

from pytint.input import read_yamlfile
import pytint.queue as pq
import argparse as ap
import pytint.lattice as pl

def spawn_jobs(inputfile, monitor=False):
    
    options = read_yamlfile(inputfile)
    #gather job arrays
    temp = options["main"]["temperature"]
    press = options["main"]["pressure"]
    lattice = options["main"]["lattice"]
    conc = options["main"]["concentration"]
    element = options["main"]["element"]

    #Check lattice values
    lattice = [x.upper() for x in lattice]
    unidentified = [l for l in lattice if l not in ["BCC", "FCC", "HCP", "DIA", "SC"]]
    if len(unidentified) > 0:
        raise ValueError("Unknown lattice found. Allowed options are BCC, FCC, HCP, DIA or SC.")

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

    lattice_constants, atoms_per_cell, lammps_lattice = pl.get_lattice(element, lattice)

    #main looping starts
    for count, l in enumerate(lattice):
        for t in temp:
            for p in press:
                for c in conc:
                    ts = int(t)
                    ps = "%.02f"%p
                    cs = "%.02f"%c
                    identistring = "-".join([l, str(ts), ps, cs])
                    scriptpath = os.path.join(os.getcwd(), ".".join([identistring, "sub"]))
                    errfile = os.path.join(os.getcwd(), ".".join([identistring, "err"]))
                    errfiles.append(errfile)

                    #get the other info which is required
                    apc = atoms_per_cell[count]
                    a = lattice_constants[count]
                    ml = lammps_lattice[count]

                    #for lattice just provide the number of position
                    scheduler.maincommand = "tint_kernel -i %s -t %f -p %f -l %s -apc %d -a %f -c %f -m %s"%(inputfile, 
                        t, p, l, apc, a, c, ml)
                    scheduler.write_script(scriptpath)
                    _ = scheduler.submit()


    if monitor:
        raise NotImplementedError("feature not implemented")

def integrate(inputfile):
    """
    Integrate results
    """
    options = read_yamlfile(inputfile)
    #gather job arrays
    temp = options["main"]["temperature"]
    press = options["main"]["pressure"]
    lattice = options["main"]["lattice"]
    conc = options["main"]["concentration"]

    #here we have to loop first
    total = 0
    success = 0

    for count, l in enumerate(lattice):
        out_t = []
        out_p = []
        out_c = []
        out_fe = []
        out_ferr = []
        for t in temp:
            for p in press:
                for c in conc:
                    ts = int(t)
                    ps = "%.02f"%p
                    cs = "%.02f"%c
                    identistring = "-".join([l, str(ts), ps, cs])
                    repfile = os.path.join(os.getcwd(), identistring, "report.yaml")
                    total += 1
                    if not os.path.exists(repfile):
                        print("%s not found, skipping.."%identistring)

                    with open(repfile, 'r') as fout:
                        data = yaml.load(fout, Loader=yaml.FullLoader)

                    out_t.append(t)
                    out_p.append(p)
                    out_c.append(c)
                    out_fe.append(data["fe"])
                    out_ferr.append(data["fe_err"])
                    success += 1
        x = np.column_stack((out_t, out_p, out_c, out_fe, out_ferr))
        outfile = os.path.join(os.getcwd(), "fe_%s.dat"%l)
        np.savetxt(outfile, x, fmt=('%.02f', '%.02f', 
                        '%.02f', '%.05f', '%.05f'),
                        header="temperature pressure concentration fe fe_err")
    print("%d/%d results saved."%(success, total))

def main():
    arg = ap.ArgumentParser()
    
    #argument name of input file
    arg.add_argument("-i", "--input", required=True, type=str,
    help="name of the input file")
    
    arg.add_argument("-m", "--mode", required=False, choices=["run", "integrate"],
    default="run", help="name of the input file")

    args = vars(arg.parse_args())
    
    if args["mode"] == "run":
        spawn_jobs(args["input"])
    else:
        integrate(args["input"]) 
