"""
Contains methods for solid
"""

import numpy as np
import yaml
from pytint.integrators import *
from pylammpsmpi import LammpsLibrary
import pyscal.core as pc
from scipy.integrate import cumtrapz

class Solid:
    """
    Solid class
    """
    def __init__(self, t=None, p=None, l=None, apc=None,
                    alat=None, c=None, options=None, simfolder=None):
        """
        Set up class
        """
        self.t = t
        self.p = p
        self.l = l
        self.apc = apc
        self.alat = alat
        self.c = c
        self.options = options
        self.simfolder = simfolder

    def run_averaging(self):
        """
        Write averagin script for solid
        """

        cores = self.options["queue"]["cores"]
        lmp = LammpsLibrary(mode="local", cores=cores)

        lmp.command("units            metal")
        lmp.command("boundary         p p p")
        lmp.command("atom_style       atomic")
        lmp.command("timestep         %f"%self.options["md"]["timestep"])

        lmp.command("lattice          %s %f"%(self.l, self.alat))
        lmp.command("region           box block 0 %d 0 %d 0 %d"%(self.options["md"]["nx"], self.options["md"]["ny"], self.options["md"]["nz"]))
        lmp.command("create_box       1 box")
        lmp.command("create_atoms     1 box")

        lmp.command("pair_style       %s"%self.options["md"]["pair_style"])
        lmp.command("pair_coeff       %s"%self.options["md"]["pair_coeff"])

        lmp.command("mass             * %f"%self.options["md"]["mass"])

        lmp.command("velocity         all create %f %d"%(self.t, np.random.randint(0, 10000)))
        lmp.command("fix              1 all npt temp %f %f %f iso %f %f %f"%(self.t, self.t, self.options["md"]["tdamp"], 
                                            0, self.p, self.options["md"]["pdamp"]))
        lmp.command("thermo_style     custom step pe press vol etotal temp")
        lmp.command("thermo           10")
        lmp.command("run              %d"%int(self.options["md"]["nsmall"])) 
        lmp.command("unfix            1")

        lmp.command("velocity         all create %f %d"%(self.t, np.random.randint(0, 10000)))
        lmp.command("fix              1 all npt temp %f %f %f iso %f %f %f"%(self.t, self.t, self.options["md"]["tdamp"], 
                                            self.p, self.p, self.options["md"]["pdamp"]))
        lmp.command("run              %d"%int(self.options["md"]["nsmall"])) 

        lmp.command("fix              2 all print 10 \"$(step) $(press) $(vol) $(temp)\" file avg.dat")
        lmp.command("run              %d"%int(self.options["md"]["nlarge"]))

        #now run for msd
        lmp.command("unfix            1")
        lmp.command("unfix            2")

        lmp.command("fix              3 all nvt temp %f %f %f"%(self.t, self.t, self.options["md"]["tdamp"]))
        lmp.command("compute          1 all msd com yes")

        lmp.command("variable         msd equal c_1[4]")

        lmp.command("fix              4 all print 10 \"$(step) ${msd}\" file msd.dat")
        lmp.command("run              %d"%int(self.options["md"]["nlarge"]))
        lmp.close()

    def gather_average_data(self):
        """
        Gather average daya
        """
        avgfile = os.path.join(self.simfolder, "avg.dat")
        vol = np.loadtxt(avgfile, usecols=(2,), unpack=True)
        avgvol = np.mean(vol[-100:])
        ncells = self.options["md"]["nx"]*self.options["md"]["ny"]*self.options["md"]["nz"]
        self.natoms = ncells*self.apc

        self.avglat = (avgvol/ncells)**(1/3)
        msdfile = os.path.join(self.simfolder, "msd.dat")
        msd = np.loadtxt(msdfile, usecols=(1,), unpack=True)
        self.k = 3*kb*self.t/np.mean(msd[-100:])

    def run_integration(self, iteration=1):
        """
        Write TI integrate script
        """
        cores = self.options["queue"]["cores"]
        lmp = LammpsLibrary(mode="local", cores=cores)

        lmp.command("echo              log")

        lmp.command("variable          T0 equal 0.7*%f"%self.t)  # Initial temperature.
        lmp.command("variable          te equal %d"%self.options["md"]["te"])   # Equilibration time (steps).
        lmp.command("variable          ts equal %d"%self.options["md"]["ts"])  # Switching time (steps).
        lmp.command("variable          k equal %f"%self.k)
        lmp.command("variable          rand equal %d"%np.random.randint(0, 1000))


        #-------------------------- Atomic Setup --------------------------------------#  
        lmp.command("units            metal")
        lmp.command("boundary         p p p")
        lmp.command("atom_style       atomic")

        lmp.command("lattice          %s %f"%(self.l, self.alat))
        lmp.command("region           box block 0 %d 0 %d 0 %d"%(self.options["md"]["nx"], self.options["md"]["ny"], self.options["md"]["nz"]))
        lmp.command("create_box       1 box")
        lmp.command("create_atoms     1 box")


        lmp.command("neigh_modify    every 1 delay 0 check yes once no")

        lmp.command("pair_style       %s"%self.options["md"]["pair_style"])
        lmp.command("pair_coeff       %s"%self.options["md"]["pair_coeff"])
        lmp.command("mass             * %f"%self.options["md"]["mass"])

        #---------------------- Thermostat & Barostat ---------------------------------#
        lmp.command("fix               f1 all nve")
        lmp.command("fix               f2 all ti/spring 10.0 1000 1000 function 2")
        lmp.command("fix               f3 all langevin %f %f %f %d zero yes"%(self.t, self.t, self.options["md"]["tdamp"], 
                                        np.random.randint(0, 10000)))

        #------------------ Computes, variables & modifications -----------------------#
        lmp.command("compute           Tcm all temp/com")
        lmp.command("fix_modify        f3 temp Tcm")

        lmp.command("variable          step    equal step")
        lmp.command("variable          dU      equal pe/atoms-f_f2/atoms")
        lmp.command("variable          lambda  equal f_f2[1]")

        lmp.command("variable          te_run  equal %d-1"%self.options["md"]["te"]) # Print correctly on fix print.
        lmp.command("variable          ts_run  equal %d+1"%self.options["md"]["ts"]) # Print correctly on fix print.

        #------------------------- Thermo stuff ---------------------------------------#
        lmp.command("thermo_style      custom step pe c_Tcm")
        lmp.command("timestep          0.001")
        lmp.command("thermo            10000")

        #------------------------- Running the Simulation -----------------------------#
        lmp.command("velocity          all create ${T0} ${rand} mom yes rot yes dist gaussian")

        lmp.command("fix               f2 all ti/spring %f %d %d function 2"%(self.k, self.options["md"]["ts"], self.options["md"]["te"]))

        # Forward. 
        lmp.command("run               ${te_run}")
        lmp.command("fix               f4 all print 1 \"${dU} ${lambda}\" screen no file forward_%d.dat "%iteration)
        lmp.command("run               ${ts_run}")
        lmp.command("unfix             f4")

        # Backward. 
        lmp.command("run               ${te_run}")
        lmp.command("fix               f4 all print 1 \"${dU} ${lambda}\" screen no file backward_%d.dat"%iteration)
        lmp.command("run               ${ts_run}")
        lmp.command("unfix             f4")

        lmp.close()


    def thermodynamic_integration(self):
        f1 = get_einstein_crystal_fe(self.t, 
            self.natoms, self.options["md"]["mass"], 
            self.avglat, self.k, self.apc)
        w, q, qerr = find_w(self.simfolder, nsims=self.options["main"]["nsims"], 
            full=True)
        fe = f1 + w
        self.fe = fe
        self.ferr = qerr


    def submit_report(self):
        report = {}
        report["temperature"] = int(self.t)
        report["pressure"] = float(self.p)
        report["lattice"] = str(self.l)
        report["concentration"] = float(self.c)
        report["avglat"] = float(self.avglat)
        report["k"] = float(self.k)
        report["fe"] = float(self.fe)
        report["fe_err"] = float(self.ferr)

        reportfile = os.path.join(self.simfolder, "report.yaml")
        with open(reportfile, 'w') as f:
            yaml.dump(report, f)


    def reversible_scaling(self, iteration=1):
        """
        Write TI integrate script
        """
        #rev scale needs tstart and tend; here self.t will be start
        #tend will be the final temp
        
        t0 = self.t
        tf = self.options["main"]["temperature"][-1]
        li = 1
        lf = t0/tf

        cores = self.options["queue"]["cores"]
        lmp = LammpsLibrary(mode="local", cores=cores)

        lmp.command("echo              log")

        lmp.command("variable          T0 equal %f"%t0)  # Initial temperature.
        lmp.command("variable          te equal %d"%self.options["md"]["te"])   # Equilibration time (steps).
        lmp.command("variable          ts equal %d"%self.options["md"]["ts"])  # Switching time (steps).
        lmp.command("variable          li equal %f"%li)
        lmp.command("variable          lf equal %f"%lf)
        lmp.command("variable          rand equal %d"%np.random.randint(0, 1000))


    #-------------------------- Atomic Setup --------------------------------------#  
        lmp.command("units            metal")
        lmp.command("boundary         p p p")
        lmp.command("atom_style       atomic")

        lmp.command("lattice          %s %f"%(self.l, self.alat))
        lmp.command("region           box block 0 %d 0 %d 0 %d"%(self.options["md"]["nx"], self.options["md"]["ny"], self.options["md"]["nz"]))
        lmp.command("create_box       1 box")
        lmp.command("create_atoms     1 box")


        lmp.command("neigh_modify    every 1 delay 0 check yes once no")

        lmp.command("pair_style       %s"%self.options["md"]["pair_style"])
        lmp.command("pair_coeff       %s"%self.options["md"]["pair_coeff"])
        lmp.command("mass             * %f"%self.options["md"]["mass"])

    #---------------------- Thermostat & Barostat ---------------------------------#
        lmp.command("fix               f1 all nph aniso %f %f %f"%(self.p, self.p, self.options["md"]["pdamp"]))
        lmp.command("fix               f2 all langevin ${T0} ${T0} %f %d zero yes"%(self.options["md"]["tdamp"], np.random.randint(0, 10000)))
        lmp.command("run               ${te}")
        lmp.command("unfix             f1")
        lmp.command("unfix             f2")

        lmp.command("variable         xcm equal xcm(all,x)")
        lmp.command("variable         ycm equal xcm(all,y)")
        lmp.command("variable         zcm equal xcm(all,z)")
        
        lmp.command("fix              f1 all nph iso %f %f %f fixedpoint ${xcm} ${ycm} ${zcm}"%(self.p, self.p, self.options["md"]["pdamp"]))
        lmp.command("fix              f2 all langevin ${T0} ${T0} %f %d zero yes"%(self.options["md"]["tdamp"], np.random.randint(0, 10000)))
        
    #------------------ Computes, variables & modifications -----------------------#
        lmp.command("compute           tcm all temp/com")
        lmp.command("fix_modify        f1 temp tcm")
        lmp.command("fix_modify        f2 temp tcm")

        lmp.command("variable          step    equal step")
        lmp.command("variable          dU      equal c_thermo_pe/atoms")
        lmp.command("variable          te_run  equal ${te}-1")
        lmp.command("variable          ts_run  equal ${ts}+1")
        lmp.command("thermo_style      custom step pe c_tcm")
        lmp.command("timestep          %f"%self.options["md"]["timestep"])
        lmp.command("thermo            10000")
        

        lmp.command("velocity          all create ${T0} ${rand} mom yes rot yes dist gaussian")   
        lmp.command("variable          i loop %d"%self.options["main"]["nsims"])

        lmp.command("run               ${te}")
        lmp.command("variable          lambda equal ramp(${li},${lf})")

        #we need to similar to liquid here

        lmp.command("fix               f3 all adapt 1 pair %s scale * * v_lambda"%self.options["md"]["pair_style"])
        lmp.command("fix               f4 all print 1 \"${dU} ${lambda}\" screen no file forward_%d.dat"%iteration)
        lmp.command("run               ${ts}")
        lmp.command("unfix             f3")
        lmp.command("unfix             f4")

        #check for melting
        lmp.command("dump              2 all custom 1 traj.dat id type mass x y z vx vy vz")
        lmp.command("run               0")
        lmp.command("undump            2")
        
        sys = pc.System()
        sys.read_inputfile("traj.dat")
        sys.find_neighbors(method="cutoff", cutoff=0)
        solids = sys.find_solids()
        if (solids/lmp.natoms < 0.7):
            lmp.close()
            raise RuntimeError("System melted, increase size or reduce scaling!")

        lmp.command("run               ${te}")
        lmp.command("variable          lambda equal ramp(${lf},${li})")
        lmp.command("fix               f3 all adapt 1 pair %s scale * * v_lambda"%self.options["md"]["pair_style"])
        lmp.command("fix               f4 all print 1 \"${dU} ${lambda}\" screen no file backward_%d.dat"%iteration)
        lmp.command("run               ${ts}")
        lmp.command("unfix             f3")
        lmp.command("unfix             f4")
        
        lmp.close()


    def integrate_reversible_scaling(self, f0, scale_energy=True):
        """
        Carry out the reversible scaling operation
        """
        ws = []
        for i in range(1, self.options["main"]["nsims"]+1):
            fdx, flambda = np.loadtxt(os.path.join(self.simfolder, "forward_%d.dat"%i), unpack=True, comments="#")
            bdx, blambda = np.loadtxt(os.path.join(self.simfolder, "backward_%d.dat"%i), unpack=True, comments="#")
            
            if scale_energy:
                fdx /= flambda
                bdx /= blambda
            wf = cumtrapz(fdx, flambda,initial=0)
            wb = cumtrapz(bdx[::-1], blambda[::-1],initial=0)
            w = (wf + wb) / (2*flambda)
            ws.append(w)
     
        wmean = np.mean(ws, axis=0)
        werr = np.std(ws, axis=0)
        temp = self.t/flambda
        f = f0/flambda + 1.5*kb*temp*np.log(flambda) + wmean
        outfile = os.path.join(self.simfolder, "reversible_scaling.dat")
        np.savetxt(outfile, np.column_stack((temp, f, werr)))