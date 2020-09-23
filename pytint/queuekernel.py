import os
import numpy as np
import shutil
import argparse as ap
import subprocess
import yaml

import pytint.average_scripts as pavg
import pytint.integrate_scripts as pint
from pytint.input import read_yamlfile
from pytint.integrators import *
import pyscal.core as pc
import pyscal.traj_process as ptp

def main():
    arg = ap.ArgumentParser()
    
    #argument name of input file
    arg.add_argument("-i", "--input", required=True, type=str,
    help="name of the input file")

    arg.add_argument("-t", "--temperature", required=True, type=float,
    help="temperature for the simulation")

    arg.add_argument("-c", "--concentration", required=False, type=float,
    default=0.00, help="concentration for the simulation")

    arg.add_argument("-s", "--solid", required=True, type=str,
    choices=["yes", "no"], help="is the system in solid state?")

    args = vars(arg.parse_args())
    options = read_yamlfile(args["input"])
    
    if args["solid"] == "yes":
        args["solid"] = True
        basename = os.path.join(os.getcwd(), "solid")
        identistring = ".".join(["solid", str(int(args["temperature"])), "%.02f"%args["concentration"]])
    else:
        args["solid"] = False
        basename = os.path.join(os.getcwd(), "liquid")
        identistring = ".".join(["liquid", str(int(args["temperature"])), "%.02f"%args["concentration"]])

    basedir = os.getcwd()
    simfolder = ".".join([basename, str(int(args["temperature"])), "%.02f"%args["concentration"]])
    

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

    #so far so good
    #at this point, we have to run md with lammps - using a subprocess call
    #switch folder
    os.chdir(simfolder)

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


    #at this point, run is over
    #initial stage of interaction to be done
    avgfile = os.path.join(simfolder, "avg.dat")
    vol = np.loadtxt(avgfile, usecols=(2,), unpack=True)
    avgvol = np.mean(vol[-100:])
    ncells = options["md"]["nx"]*options["md"]["ny"]*options["md"]["nz"]
    natoms = ncells*options["md"]["atoms_per_cell"]

    if args["solid"]:
        avglat = (avgvol/ncells)**(1/3)
        msdfile = os.path.join(simfolder, "msd.dat")
        msd = np.loadtxt(msdfile, usecols=(1,), unpack=True)
        k = 3*kb*args["temperature"]/np.mean(msd[-100:])
        #WARNING: Add melt checks with pyscal
    else:
        #WARNING: This can be automated - learn atoms per cell
        #automatically
        rho = natoms/avgvol
        #WARNING: hard coded ufm parameter
        eps = args["temperature"]*50.0*kb


    #initial stage is over
    #we can move on to submitting integrate scripts
    scriptpath = os.path.join(simfolder, "mdint.in")

    if args["solid"]:
        pint.write_script_solid(scriptpath, args["temperature"], 
            k, avglat, options)
    else:
        #we need to grab a conf for the liquid
        #possibly implement for solid too - we will have to see
        trajfile = os.path.join(simfolder, "traj.dat")
        files = ptp.split_trajectory(trajfile)
        conf = os.path.join(simfolder, "conf.dump")
        #now rewrite conf
        sys = pc.System()
        sys.read_inputfile(files[-1], customkeys=["vx", "vy", "vz", "mass"])
        sys.to_file(conf, customkeys=["vx", "vy", "vz", "mass"])
        #remove unnecessary files
        os.remove(trajfile)
        for file in files:
            os.remove(file)

        pint.write_script_liquid(scriptpath, args["temperature"], 
            eps, conf, options)
    
    #great! now we submit calculation again
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
    
    # now we can do final integration
    # we need string changes here
    if args["solid"]:
        f1 = get_einstein_crystal_fe(args["temperature"], 
            natoms, options["md"]["mass"], 
            avglat, k, options["md"]["atoms_per_cell"])
        w, q, qerr = find_w(simfolder, nsims=options["main"]["nsims"], 
            full=True, 
            temp=args["temperature"])
        fe = f1 + w
    else:
        w, q, qerr = find_w(simfolder, nsims=options["main"]["nsims"], 
            full=True)  
        #WARNING: hardcoded UFM parameters           
        f1 = get_uhlenbeck_ford_fe(args["temperature"], 
            rho, 50, 1.5)
        f2 = get_ideal_gas_fe(args["temperature"], rho, 
            natoms, options["md"]["mass"], xa=(1-args["concentration"]), 
            xb=args["concentration"])
        fe = f2 + f1 - w

    #thats a wrap
    #make a report
    report = {}
    report["temperature"] = args["temperature"]
    report["concentration"] = args["concentration"]
    report["solid"] = args['solid']

    if args["solid"]:
        report["avglat"] = float(avglat)
        report["k"] = float(k)
    else:
        report["rho"] = float(rho)

    report["fe"] = float(fe)
    report["fe_err"] = float(qerr)

    reportfile = os.path.join(simfolder, "report.yaml")
    with open(reportfile, 'w') as f:
        yaml.dump(report, f)

    #write a copy in home folder too
    reportfile = os.path.join(basedir, ".".join([identistring, "yaml"]))
    with open(reportfile, 'w') as f:
        yaml.dump(report, f)
