"""
Contains methods for solid
"""

import numpy as np
import yaml
from pytint.integrators import *
from pylammpsmpi import LammpsLibrary
import pyscal.core as pc
import logging

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

        logfile = os.path.join(self.simfolder, "tint.log")
        self.prepare_log(logfile)

    def prepare_log(self, file):
        logger = logging.getLogger(__name__)
        handler = logging.FileHandler(file)
        formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.DEBUG)
        logger.propagate = False
        self.logger = logger

    def run_averaging(self):
        """
        Write averagin script for solid
        """

        cores = self.options["queue"]["cores"]
        lmp = LammpsLibrary(mode="local", cores=cores, working_directory=self.simfolder)

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

        #add some computes
        lmp.command("variable         mvol equal vol")
        lmp.command("variable         mpress equal press")

        if self.p == 0:
            #This routine should be followed for zero pressure
            lmp.command("velocity         all create %f %d"%(self.t, np.random.randint(0, 10000)))
            lmp.command("fix              1 all npt temp %f %f %f iso %f %f %f"%(self.t, self.t, self.options["md"]["tdamp"], 
                                                self.p, self.p, self.options["md"]["pdamp"]))
            lmp.command("thermo_style     custom step pe press vol etotal temp")
            lmp.command("thermo           10")
            lmp.command("run              %d"%int(self.options["md"]["nsmall"])) 
            lmp.command("unfix            1")

            lmp.command("velocity         all create %f %d"%(self.t, np.random.randint(0, 10000)))
            lmp.command("fix              1 all npt temp %f %f %f iso %f %f %f"%(self.t, self.t, self.options["md"]["tdamp"], 
                                                self.p, self.p, self.options["md"]["pdamp"]))
            lmp.command("run              %d"%int(self.options["md"]["nsmall"])) 
        else:
            #Now this routine is for non-zero pressure
            #one has to equilibriate at a low temperature, but high pressure and then increase temp gradually
            #start at 0.25 temp, and increase to 0.50, while keeping high pressure
            lmp.command("velocity         all create %f %d"%(0.25*self.t, np.random.randint(0, 10000)))
            lmp.command("fix              1 all npt temp %f %f %f iso %f %f %f"%(0.25*self.t, 0.5*self.t, self.options["md"]["tdamp"], 
                                                self.p, self.p, self.options["md"]["pdamp"]))
            lmp.command("thermo_style     custom step pe press vol etotal temp")
            lmp.command("thermo           10")
            lmp.command("run              %d"%int(self.options["md"]["nsmall"])) 
            lmp.command("unfix            1")

            #now heat again
            lmp.command("fix              1 all npt temp %f %f %f iso %f %f %f"%(0.5*self.t, self.t, self.options["md"]["tdamp"], 
                                                self.p, self.p, self.options["md"]["pdamp"]))
            lmp.command("run              %d"%int(self.options["md"]["nsmall"])) 
            lmp.command("unfix            1")

            #now run normal cycle
            lmp.command("fix              1 all npt temp %f %f %f iso %f %f %f"%(self.t, self.t, self.options["md"]["tdamp"], 
                                                self.p, self.p, self.options["md"]["pdamp"]))
            lmp.command("run              %d"%int(self.options["md"]["nsmall"])) 


        #this is when the averaging routine starts
        lmp.command("fix              2 all ave/time %d %d %d v_mvol v_mpress file avg.dat"%(int(self.options["md"]["nevery"]),
            int(self.options["md"]["nrepeat"]), int(self.options["md"]["nevery"]*self.options["md"]["nrepeat"])))
        
        laststd = 0.00
        for i in range(int(self.options["md"]["ncycles"])):
            lmp.command("run              %d"%int(self.options["md"]["nsmall"]))
            ncount = int(self.options["md"]["nsmall"])//int(self.options["md"]["nevery"]*self.options["md"]["nrepeat"])
            #now we can check if it converted
            file = os.path.join(self.simfolder, "avg.dat")
            quant, ipress = np.loadtxt(file, usecols=(1,2), unpack=True)
            lx = (quant/(self.options["md"]["nx"]*self.options["md"]["ny"]*self.options["md"]["nz"]))**(1/3)
            lx = lx[-ncount+1:]
            mean = np.mean(lx)
            std = np.std(lx)
            self.logger.info("At count %d mean lattice constant is %f std is %f"%(i+1, mean, std))
            if (np.abs(laststd - std) < self.options["conv"]["alat_tol"]):
                self.avglat = np.round(mean, decimals=3)
                self.logger.info("finalized lattice constant %f pressure %f"%(self.avglat, np.mean(ipress)))
                break
            laststd = std

        #now run for msd
        lmp.command("unfix            1")
        lmp.command("unfix            2")

        #check for melting
        #check for melting
        lmp.command("dump              2 all custom 1 traj.dat id type mass x y z vx vy vz")
        lmp.command("run               0")
        lmp.command("undump            2")
        
        sys = pc.System()
        sys.read_inputfile("traj.dat")
        sys.find_neighbors(method="cutoff", cutoff=0)
        solids = sys.find_solids()
        if (solids/lmp.natoms < self.options["conv"]["solid_frac"]):
            lmp.close()
            raise RuntimeError("System melted, increase size or reduce temp!")


        lmp.command("fix              3 all nvt temp %f %f %f"%(self.t, self.t, self.options["md"]["tdamp"]))
        lmp.command("compute          1 all msd com yes")

        lmp.command("variable         msd equal c_1[4]")

        #we need a similar averaging routine here
        lmp.command("fix              4 all ave/time %d %d %d v_msd file msd.dat"%(int(self.options["md"]["nevery"]),
            int(self.options["md"]["nrepeat"]), int(self.options["md"]["nevery"]*self.options["md"]["nrepeat"])))
        laststd = 0.00
        for i in range(self.options["md"]["ncycles"]):
            lmp.command("run              %d"%int(self.options["md"]["nsmall"]))
            ncount = int(self.options["md"]["nsmall"])//int(self.options["md"]["nevery"]*self.options["md"]["nrepeat"])
            #now we can check if it converted
            file = os.path.join(self.simfolder, "msd.dat")
            quant = np.loadtxt(file, usecols=(1,), unpack=True)[-ncount+1:]
            quant = 3*kb*self.t/quant
            #self.logger.info(quant)
            mean = np.mean(quant)
            std = np.std(quant)
            self.logger.info("At count %d mean k is %f std is %f"%(i+1, mean, std))
            if (np.abs(laststd - std) < self.options["conv"]["k_tol"]):
                self.k = np.round(mean, decimals=2)
                self.logger.info("finalized sprint constant %f"%(self.k))
                break
            laststd = std

        #check for melting
        lmp.command("dump              2 all custom 1 traj.dat id type mass x y z vx vy vz")
        lmp.command("run               0")
        lmp.command("undump            2")
        
        sys = pc.System()
        sys.read_inputfile("traj.dat")
        sys.find_neighbors(method="cutoff", cutoff=0)
        solids = sys.find_solids()
        if (solids/lmp.natoms < self.options["conv"]["solid_frac"]):
            lmp.close()
            raise RuntimeError("System melted, increase size or reduce temp!")

        lmp.close()

        #now some housekeeping
        ncells = self.options["md"]["nx"]*self.options["md"]["ny"]*self.options["md"]["nz"]
        self.natoms = ncells*self.apc


    def run_integration(self, iteration=1):
        """
        Write TI integrate script
        """
        cores = self.options["queue"]["cores"]
        lmp = LammpsLibrary(mode="local", cores=cores, working_directory=self.simfolder)

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

        lmp.command("lattice          %s %f"%(self.l, self.avglat))
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
        lmp.command("variable          dU1      equal pe/atoms")
        lmp.command("variable          dU2      equal f_f2/atoms")
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
        lmp.command("fix               f4 all print 1 \"${dU1} ${dU2} ${lambda}\" screen no file forward_%d.dat "%iteration)
        lmp.command("run               ${ts_run}")
        lmp.command("unfix             f4")

        # Backward. 
        lmp.command("run               ${te_run}")
        lmp.command("fix               f4 all print 1 \"${dU1} ${dU2} ${lambda}\" screen no file backward_%d.dat"%iteration)
        lmp.command("run               ${ts_run}")
        lmp.command("unfix             f4")

        lmp.command("unfix             f2")

        lmp.close()


    def thermodynamic_integration(self):
        f1 = get_einstein_crystal_fe(self.t, 
            self.natoms, self.options["md"]["mass"], 
            self.avglat, self.k, self.apc)
        w, q, qerr = find_w(self.simfolder, nsims=self.options["main"]["nsims"], 
            full=True, solid=True)
        
        self.fref = f1
        self.w = w
        self.ferr = qerr

        if self.p != 0:
            #add pressure contribution
            p = self.p/(10000*160.21766208)
            v = (self.avglat**3)/self.apc
            self.pv = p*v
        else:
            self.pv = 0 

        self.fe = self.fref + self.w + self.pv


    def submit_report(self):

        report = {}

        #input quantities
        report["input"] = {}
        report["input"]["temperature"] = int(self.t)
        report["input"]["pressure"] = float(self.p)
        report["input"]["lattice"] = str(self.l)
        report["input"]["concentration"] = float(self.c)

        #average quantities
        report["average"] = {}
        report["average"]["lattice_constant"] = float(self.avglat)
        report["average"]["spring_constant"] = float(self.k)
        
        #results
        report["results"] = {}
        report["results"]["free_energy"] = float(self.fe)
        report["results"]["error"] = float(self.ferr)
        report["results"]["reference_system"] = float(self.fref)
        report["results"]["work"] = float(self.w)
        report["results"]["pv"] = float(self.pv)

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
        lmp = LammpsLibrary(mode="local", cores=cores, working_directory=self.simfolder)

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

        lmp.command("lattice          %s %f"%(self.l, self.avglat))
        lmp.command("region           box block 0 %d 0 %d 0 %d"%(self.options["md"]["nx"], self.options["md"]["ny"], self.options["md"]["nz"]))
        lmp.command("create_box       1 box")
        lmp.command("create_atoms     1 box")


        lmp.command("neigh_modify    every 1 delay 0 check yes once no")

        lmp.command("pair_style       %s"%self.options["md"]["pair_style"])
        lmp.command("pair_coeff       %s"%self.options["md"]["pair_coeff"])
        lmp.command("mass             * %f"%self.options["md"]["mass"])

        lmp.command("velocity          all create ${T0} ${rand} mom yes rot yes dist gaussian")

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


    def integrate_reversible_scaling(self, scale_energy=False):
        """
        Carry out the reversible scaling operation
        """
        integrate_rs(self.simfolder, self.fe, self.t,
            nsims=self.options["main"]["nsims"], scale_energy=scale_energy)