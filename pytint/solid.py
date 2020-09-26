"""
Contains methods for solid
"""

import numpy as np
import yaml
from pytint.integrators import *


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

    def write_average_script(self, mdscriptfile):
        """
        Write averagin script for solid
        """
        with open(mdscriptfile, 'w') as fout:
            fout.write("units            metal\n")
            fout.write("boundary         p p p\n")
            fout.write("atom_style       atomic\n")
            fout.write("timestep         %f\n"%self.options["md"]["timestep"])

            fout.write("lattice          %s %f\n"%(self.l, self.alat))
            fout.write("region           box block 0 %d 0 %d 0 %d\n"%(self.options["md"]["nx"], self.options["md"]["ny"], self.options["md"]["nz"]))
            fout.write("create_box       1 box\n")
            fout.write("create_atoms     1 box\n")

            fout.write("pair_style       %s\n"%self.options["md"]["pair_style"])
            fout.write("pair_coeff       %s\n"%self.options["md"]["pair_coeff"])

            fout.write("mass             * %f\n"%self.options["md"]["mass"])

            fout.write("velocity         all create %f %d\n"%(1, np.random.randint(0, 10000)))
            fout.write("fix              1 all npt temp %f %f %f iso %f %f %f\n"%(1, self.t, self.options["md"]["tdamp"], 
                                                0, self.p, self.options["md"]["pdamp"]))
            fout.write("thermo_style     custom step pe press vol etotal temp\n")
            fout.write("thermo           10\n")
            fout.write("run              %d\n"%int(self.options["md"]["nsmall"])) 
            fout.write("unfix            1\n")

            fout.write("velocity         all create %f %d\n"%(self.t, np.random.randint(0, 10000)))
            fout.write("fix              1 all npt temp %f %f %f iso %f %f %f\n"%(self.t, self.t, self.options["md"]["tdamp"], 
                                                self.p, self.p, self.options["md"]["pdamp"]))
            fout.write("run              %d\n"%int(self.options["md"]["nsmall"])) 

            fout.write("fix              2 all print 10 \"$(step) $(press) $(vol) $(temp)\" file avg.dat\n")
            fout.write("run              %d\n"%int(self.options["md"]["nlarge"]))

            #now run for msd
            fout.write("unfix            1\n")
            fout.write("unfix            2\n")

            fout.write("fix              3 all nvt temp %f %f %f\n"%(self.t, self.t, self.options["md"]["tdamp"]))
            fout.write("compute          1 all msd com yes\n")

            fout.write("variable         msd equal c_1[4]\n")

            fout.write("fix              4 all print 10 \"$(step) ${msd}\" file msd.dat\n")
            fout.write("run              %d\n"%int(self.options["md"]["nlarge"]))


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

    def write_integrate_script(self, mdscriptfile):
        """
        Write TI integrate script
        """
        with open(mdscriptfile, 'w') as fout:
            fout.write("echo              log\n")

            fout.write("variable          T0 equal 0.7*%f\n"%self.t)  # Initial temperature.
            fout.write("variable          te equal %d\n"%self.options["md"]["te"])   # Equilibration time (steps).
            fout.write("variable          ts equal %d\n"%self.options["md"]["ts"])  # Switching time (steps).
            fout.write("variable          k equal %f\n"%self.k)
            fout.write("variable          rand equal %d\n"%np.random.randint(0, 1000))


        #-------------------------- Atomic Setup --------------------------------------#  
            fout.write("units            metal\n")
            fout.write("boundary         p p p\n")
            fout.write("atom_style       atomic\n")

            fout.write("lattice          %s %f\n"%(self.l, self.alat))
            fout.write("region           box block 0 %d 0 %d 0 %d\n"%(self.options["md"]["nx"], self.options["md"]["ny"], self.options["md"]["nz"]))
            fout.write("create_box       1 box\n")
            fout.write("create_atoms     1 box\n")


            fout.write("neigh_modify    every 1 delay 0 check yes once no\n")

            fout.write("pair_style       %s\n"%self.options["md"]["pair_style"])
            fout.write("pair_coeff       %s\n"%self.options["md"]["pair_coeff"])
            fout.write("mass             * %f\n"%self.options["md"]["mass"])

        #---------------------- Thermostat & Barostat ---------------------------------#
            fout.write("fix               f1 all nve\n")
            fout.write("fix               f2 all ti/spring 10.0 1000 1000 function 2\n")
            fout.write("fix               f3 all langevin %f %f %f %d zero yes\n"%(self.t, self.t, self.options["md"]["tdamp"], 
                                            np.random.randint(0, 10000)))

        #------------------ Computes, variables & modifications -----------------------#
            fout.write("compute           Tcm all temp/com\n")
            fout.write("fix_modify        f3 temp Tcm\n")

            fout.write("variable          step    equal step\n")
            fout.write("variable          dU      equal pe/atoms-f_f2/atoms\n")
            fout.write("variable          lambda  equal f_f2[1]\n")

            fout.write("variable          te_run  equal %d-1\n"%self.options["md"]["te"]) # Print correctly on fix print.
            fout.write("variable          ts_run  equal %d+1\n"%self.options["md"]["ts"]) # Print correctly on fix print.

        #------------------------- Thermo stuff ---------------------------------------#
            fout.write("thermo_style      custom step pe c_Tcm\n")
            fout.write("timestep          0.001\n")
            fout.write("thermo            10000\n")

        #------------------------- Running the Simulation -----------------------------#
            fout.write("velocity          all create ${T0} ${rand} mom yes rot yes dist gaussian\n")


            fout.write("variable          i loop %d\n"%self.options["main"]["nsims"])
            fout.write("label repetitions\n")
            fout.write("fix               f2 all ti/spring %f %d %d function 2\n"%(self.k, self.options["md"]["ts"], self.options["md"]["te"]))

            # Forward. 
            fout.write("run               ${te_run}\n")
            fout.write("fix               f4 all print 1 \"${dU} ${lambda}\" screen no file forward_$i.dat \n")
            fout.write("run               ${ts_run}\n")
            fout.write("unfix             f4\n")

            # Backward. 
            fout.write("run               ${te_run}\n")
            fout.write("fix               f4 all print 1 \"${dU} ${lambda}\" screen no file backward_$i.dat\n")
            fout.write("run               ${ts_run}\n")
            fout.write("unfix             f4\n")

            fout.write("next i\n")
            fout.write("jump SELF repetitions\n")


    def thermodynamic_integration(self):
        f1 = get_einstein_crystal_fe(self.t, 
            self.natoms, self.options["md"]["mass"], 
            self.avglat, self.k, self.apc)
        w, q, qerr = find_w(self.simfolder, nsims=self.options["main"]["nsims"], 
            full=True, 
            temp=self.t)
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
