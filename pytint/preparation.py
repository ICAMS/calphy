from pylammpsmpi import LammpsLibrary
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
import pyscal.core as pc
import os
from scipy.constants import physical_constants
import pandas as pd
import shutil

kb = physical_constants["Boltzmann constant in eV/K"]

def find_thermo_solid(temp, **kwargs):
    """
    Find the average lattice constant and spring constant
    for the solid structure at the required temperature.
    """
    timestep = kwargs.get('timestep', 0.004)
    alat = 
    nx = 
    ny =
    nz = 
    short_runs =
    short_freq =
    structure =
    
    timestep = 0.004
    alat = 4.0
    nx = 10
    ny = 10
    nz = 10
    short_runs = 10000
    short_freq = 100
    short_intv = short_runs//short_freq 
    lmp = LammpsLibrary(cores=20, mode='local')
    lmp.echo("log")
    lmp.units("metal")
    lmp.boundary("p p p")
    lmp.atom_style("atomic")
    lmp.timestep(timestep)
    lmp.lattice("fcc", alat)
    lmp.region("box", "block", 0, nx, 0, ny, 0, nz)
    lmp.create_box(2, "box")
    lmp.command("create_atoms 1 box")
    lmp.variable("fract  atom random(0,1,9876)<=%s"%str(conc))
    lmp.group("part variable fract")
    lmp.set("group part type 2")
    lmp.pair_style("sw")
    lmp.pair_coeff("* * AfAd.sw Af Ad")
    lmp.neighbor("1.0 bin")
    lmp.neigh_modify("every 1  delay 1  check yes")
    lmp.mass("*  28.08")
    lmp.velocity("all create", temp, np.random.randint(0, 10000))
    lmp.fix("1 all npt temp", temp, temp, 0.1, "iso", 0.0, 0.0, 1.0) 
    temp1 = []
    vol1 = []
    for i in tqdm(range(short_intv)):
        lmp.run(short_freq)
        temp1.append(lmp.temp)
        vol1.append(lmp.vol)
    avgvol = np.mean(vol1[:30])
    avglat = (avgvol/(nx*ny*nz))**(1/3)
    lmp.dump("1 all custom 100 traj3.dat id type mass x y z vx vy vz")
    lmp.run("0")
    lmp.undump(1)
    lmp.compute("1 all msd com yes")
    
    temp2 = []
    msd2 = []
    for i in tqdm(range(short_intv)):
        lmp.run(short_freq)
        temp2.append(lmp.temp)
        mm = lmp.extract_compute("1", 0, 1, length=4)
        msd2.append(mm[3])  
    k = 3*kb[0]*temp/np.mean(msd2)
    lmp.close()
    tfile = ".".join(["conf", str(temp), str(conc), "dump"])
    os.system("mv traj3.dat %s"%tfile)
    return avglat, k