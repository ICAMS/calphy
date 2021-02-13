"""
Contains methods for Solid thermodynamic integration and
reversible scaling.

"""

import numpy as np
import yaml

from pytint.integrators import *
import pyscal.traj_process as ptp
import pytint.lattice as pl
import pytint.helpers as ph

class Solid:
    """
    Solid method class

    Parameters
    ----------
    t : float
        Temperature for the calculation
        Unit: K

    p : float
        pressure for the calculation
        Unit: bar

    l : string
        lattice to be used for the calculation

    apc : int
        number of atoms in a single unit cell
    
    alat : float
        lattice constant

    options : dict
        dict of input options

    simfolder : string
        base folder for running calculations

    Attributes
    ----------
    t : float
        temperature

    p : float
        pressure
    
    l : string
        lattice

    apc : int
        number of atoms in unit cell

    alat : float
        lattice constant

    c : float
        concentration

    options : dict
        dict containing options

    simfolder : string
        main simulation directory

    """
    def __init__(self, options=None, kernel=None, simfolder=None):
        self.options = options
        self.simfolder = simfolder
        self.kernel = kernel

        self.calc = options["calculations"][kernel]
        self.nsims = self.calc["nsims"]

        self.t = self.calc["temperature"]
        self.tend = self.calc["temperature_stop"]
        self.thigh = self.calc["thigh"] 
        self.p = self.calc["pressure"]

        self.l = None
        self.alat = None
        self.apc = None
        self.vol = None
        self.concentration = None
        self.prepare_lattice()

        logfile = os.path.join(self.simfolder, "tint.log")
        self.logger = ph.prepare_log(logfile)

        #other properties
        self.cores = self.options["queue"]["cores"]
        self.ncells = np.prod(self.calc["repeat"])
        self.natoms = self.ncells*self.apc        
        
        #properties that will be calculated later
        self.volatom = None
        self.k = None
        self.ferr = None
        self.fref = None
        self.fideal = None
        self.w = None
        self.pv = None
        self.fe = None

        #box dimensions that need to be stored
        self.lx = None
        self.ly = None
        self.lz = None

    def prepare_lattice(self):
        #process lattice
        l, alat, apc, conc = pl.prepare_lattice(self.calc)
        self.l = l
        self.alat = alat
        self.apc = apc
        self.concentration = conc

    def run_averaging(self):
        """
        Run averaging routine

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        lmp = ph.create_object(self.cores, self.simfolder, self.options["md"]["timestep"])

        #set up structure
        lmp = ph.create_structure(lmp, self.calc)

        #set up potential
        lmp = ph.set_potential(lmp, self.options)

        #add some computes
        lmp.command("variable         mvol equal vol")
        lmp.command("variable         mlx equal lx")
        lmp.command("variable         mly equal ly")
        lmp.command("variable         mlz equal lz")
        lmp.command("variable         mpress equal press")

        if self.p == 0:
            #This routine should be followed for zero pressure
            lmp.command("velocity         all create %f %d"%(self.t, np.random.randint(0, 10000)))
            lmp.command("fix              1 all npt temp %f %f %f iso %f %f %f"%(self.t, self.t, self.options["md"]["tdamp"], 
                                                self.p, self.p, self.options["md"]["pdamp"]))
            lmp.command("thermo_style     custom step pe press vol etotal temp lx ly lz")
            lmp.command("thermo           10")
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
        lmp.command("fix              2 all ave/time %d %d %d v_mlx v_mly v_mlz v_mpress file avg.dat"%(int(self.options["md"]["nevery"]),
            int(self.options["md"]["nrepeat"]), int(self.options["md"]["nevery"]*self.options["md"]["nrepeat"])))
        
        laststd = 0.00
        for i in range(int(self.options["md"]["ncycles"])):
            lmp.command("run              %d"%int(self.options["md"]["nsmall"]))
            ncount = int(self.options["md"]["nsmall"])//int(self.options["md"]["nevery"]*self.options["md"]["nrepeat"])
            #now we can check if it converted
            file = os.path.join(self.simfolder, "avg.dat")
            lx, ly, lz, ipress = np.loadtxt(file, usecols=(1, 2, 3, 4), unpack=True)
            
            #lxpc = ((lx*ly*lz)/self.ncells)**(1/3)
            #lxpc = ipress[-ncount+1:]
            lxpc = ipress
            mean = np.mean(lxpc)
            std = np.std(lxpc)
            volatom = np.mean((lx*ly*lz)/self.natoms)
            self.logger.info("At count %d mean pressure is %f with %f vol/atom"%(i+1, mean, volatom))
            
            #if (np.abs(laststd - std) < self.options["conv"]["alat_tol"]):
            if (np.abs(mean - self.p)) < self.options["conv"]["p_tol"]:

                #process other means
                self.lx = np.round(np.mean(lx[-ncount+1:]), decimals=3)
                self.ly = np.round(np.mean(ly[-ncount+1:]), decimals=3)
                self.lz = np.round(np.mean(lz[-ncount+1:]), decimals=3)
                self.volatom = volatom
                self.vol = self.lx*self.ly*self.lz
                self.logger.info("finalized vol/atom %f at pressure %f"%(self.volatom, mean))
                self.logger.info("Avg box dimensions x: %f, y: %f, z:%f"%(self.lx, self.ly, self.lz))
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
        
        solids = ph.find_solid_fraction("traj.dat")
        if (solids/lmp.natoms < self.options["conv"]["solid_frac"]):
            lmp.close()
            raise RuntimeError("System melted, increase size or reduce temp!")


        lmp.command("fix              3 all nvt temp %f %f %f"%(self.t, self.t, self.options["md"]["tdamp"]))
        
        lmp = ph.compute_msd(lmp, self.options)
        #we need a similar averaging routine here
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
                #now reevaluate spring constants
                k = []
                for i in range(self.options["nelements"]):
                    quant = np.loadtxt(file, usecols=(i+1, ), unpack=True)[-ncount+1:]
                    quant = 3*kb*self.t/quant
                    k.append(np.round(np.mean(quant), decimals=2))

                self.k = k
                self.logger.info("finalized sprint constants")
                self.logger.info(self.k)
                break
            laststd = std

        #check for melting
        lmp.command("dump              2 all custom 1 traj.dat id type mass x y z vx vy vz")
        lmp.command("run               0")
        lmp.command("undump            2")
        
        solids = ph.find_solid_fraction("traj.dat")
        if (solids/lmp.natoms < self.options["conv"]["solid_frac"]):
            lmp.close()
            raise RuntimeError("System melted, increase size or reduce temp!")

        lmp.close()
        self.process_traj()


    def process_traj(self):
        """
        Process the out trajectory after averaging cycle and 
        extract a configuration to run integration

        Parameters
        ----------
        None

        Returns
        -------
        None
        
        """
        trajfile = os.path.join(self.simfolder, "traj.dat")
        files = ptp.split_trajectory(trajfile)
        conf = os.path.join(self.simfolder, "conf.dump")

        ph.reset_timestep(files[-1], conf)

        os.remove(trajfile)
        for file in files:
            os.remove(file)


    def run_integration(self, iteration=1):
        """
        Run integration routine

        Parameters
        ----------
        iteration : int, optional
            iteration to run, default 1

        Returns
        -------
        None
        """
        lmp = ph.create_object(self.cores, self.simfolder, self.options["md"]["timestep"])


        lmp.command("variable          T0 equal 0.7*%f"%self.t)  # Initial temperature.
        lmp.command("variable          te equal %d"%self.options["md"]["te"])   # Equilibration time (steps).
        lmp.command("variable          ts equal %d"%self.options["md"]["ts"])  # Switching time (steps).
        #lmp.command("variable          k equal %f"%self.k)
        lmp.command("variable          rand equal %d"%np.random.randint(0, 1000))


        conf = os.path.join(self.simfolder, "conf.dump")
        lmp = ph.read_dump(lmp, conf, species=self.options["nelements"])

        #set up potential
        lmp = ph.set_potential(lmp, self.options)

        #remap the box to get the correct pressure
        lmp = ph.remap_box(lmp, self.lx, self.ly, self.lz)

        #create groups
        for i in range(self.options["nelements"]):
            lmp.command("group  g%d type %d"%(i+1, i+1))

        for i in range(self.options["nelements"]):
            lmp.command("variable   count%d equal count(g%d)"%(i+1, i+1))

        lmp.command("run               0")
        #---------------------- Thermostat & Barostat ---------------------------------#
        lmp.command("fix               f1 all nve")
        
        #nelements 
        for i in range(self.options["nelements"]):
            lmp.command("fix               ff%d g%d ti/spring 10.0 1000 1000 function 2"%(i+1, i+1))
        
        lmp.command("fix               f3 all langevin %f %f %f %d zero yes"%(self.t, self.t, self.options["md"]["tdamp"], 
                                        np.random.randint(0, 10000)))

        #------------------ Computes, variables & modifications -----------------------#
        lmp.command("compute           Tcm all temp/com")
        lmp.command("fix_modify        f3 temp Tcm")

        lmp.command("variable          step    equal step")
        
        #nelements
        lmp.command("variable          dU1      equal pe/atoms")
        for i in range(self.options["nelements"]):
            lmp.command("variable          dU%d      equal f_ff%d/v_count%d"%(i+2, i+1, i+1))
        
        lmp.command("variable          lambda  equal f_ff1[1]")

        lmp.command("variable          te_run  equal %d-1"%self.options["md"]["te"]) # Print correctly on fix print.
        lmp.command("variable          ts_run  equal %d+1"%self.options["md"]["ts"]) # Print correctly on fix print.

        #------------------------- Thermo stuff ---------------------------------------#
        lmp.command("thermo_style      custom step pe c_Tcm")
        lmp.command("timestep          0.001")
        lmp.command("thermo            10000")

        #------------------------- Running the Simulation -----------------------------#
        lmp.command("velocity          all create ${T0} ${rand} mom yes rot yes dist gaussian")

        for i in range(self.options["nelements"]):
            lmp.command("fix               ff%d g%d ti/spring %f %d %d function 2"%(i+1, i+1, self.k[i], 
                self.options["md"]["ts"], self.options["md"]["te"]))

        # Forward. 
        lmp.command("run               ${te_run}")
        str1 = "fix f4 all print 1 \"${dU1} "
        str2 = []
        for i in range(self.options["nelements"]):
            str2.append("${dU%d}"%(i+2))
        str2.append("${lambda}\"")
        str2 = " ".join(str2)
        str3 = " screen no file forward_%d.dat"%iteration
        command = str1 + str2 + str3
        lmp.command(command)
        lmp.command("run               ${ts_run}")
        lmp.command("unfix             f4")

        # Backward. 
        lmp.command("run               ${te_run}")
        str1 = "fix f4 all print 1 \"${dU1} "
        str2 = []
        for i in range(self.options["nelements"]):
            str2.append("${dU%d}"%(i+2))
        str2.append("${lambda}\"")
        str2 = " ".join(str2)
        str3 = " screen no file backward_%d.dat"%iteration
        command = str1 + str2 + str3
        lmp.command(command)
        lmp.command("run               ${ts_run}")
        lmp.command("unfix             f4")

        #lmp.command("unfix             f2")

        lmp.close()


    def thermodynamic_integration(self):
        """
        Calculate free energy

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        f1 = get_einstein_crystal_fe(self.t, 
            self.natoms, self.options["mass"], 
            self.vol, self.k, self.concentration)
        w, q, qerr = find_w(self.simfolder, nelements=self.options["nelements"], concentration=self.concentration, nsims=self.nsims, 
            full=True, solid=True)
        
        self.fref = f1
        self.w = w
        self.ferr = qerr

        if self.p != 0:
            #add pressure contribution
            p = self.p/(10000*160.21766208)
            v = self.vol/self.natoms
            self.pv = p*v
        else:
            self.pv = 0 

        self.fe = self.fref + self.w + self.pv


    def submit_report(self):
        """
        Submit final report containing results

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        report = {}

        #input quantities
        report["input"] = {}
        report["input"]["temperature"] = int(self.t)
        report["input"]["pressure"] = float(self.p)
        report["input"]["lattice"] = str(self.l)
        report["input"]["element"] = " ".join(np.array(self.options["element"]).astype(str))
        report["input"]["concentration"] = " ".join(np.array(self.concentration).astype(str))

        #average quantities
        report["average"] = {}
        report["average"]["vol/atom"] = float(self.volatom)
        report["average"]["spring_constant"] = " ".join(np.array(self.k).astype(str))
        
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
        Perform reversible scaling calculation in NPT

        Parameters
        ----------
        iteration : int, optional
            iteration of the calculation. Default 1

        Returns
        -------
        None
        """
        #rev scale needs tstart and tend; here self.t will be start
        #tend will be the final temp
        # solid cannot go directly to nph/langevin
        # pressure needs to be scaled up from initial structure- then temperature
        
        t0 = self.t
        tf = self.tend
        li = 1
        lf = t0/tf
        pi = self.p
        pf = lf*pi

        #create lammps object
        lmp = ph.create_object(self.cores, self.simfolder, self.options["md"]["timestep"])

        lmp.command("echo              log")
        lmp.command("variable          T0 equal %f"%t0)  # Initial temperature.
        lmp.command("variable          te equal %d"%self.options["md"]["te"])   # Equilibration time (steps).
        lmp.command("variable          ts equal %d"%self.options["md"]["ts"])  # Switching time (steps).
        lmp.command("variable          li equal %f"%li)
        lmp.command("variable          lf equal %f"%lf)
        lmp.command("variable          rand equal %d"%np.random.randint(0, 1000))

        conf = os.path.join(self.simfolder, "conf.dump")
        lmp = ph.read_dump(lmp, conf, species=self.options["nelements"])

        #set up potential
        lmp = ph.set_potential(lmp, self.options)

        #remap the box to get the correct pressure
        lmp = ph.remap_box(lmp, self.lx, self.ly, self.lz)

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
        lmp.command("thermo_style      custom step pe c_tcm press vol")
        lmp.command("timestep          %f"%self.options["md"]["timestep"])
        lmp.command("thermo            10000")
        

        lmp.command("velocity          all create ${T0} ${rand} mom yes rot yes dist gaussian")   
        lmp.command("run               ${te}")
        lmp.command("unfix             f1")


        lmp.command("variable          lambda equal ramp(${li},${lf})")

        #we need to similar to liquid here
        lmp.command("fix              f1 all nph iso %f %f %f fixedpoint ${xcm} ${ycm} ${zcm}"%(pi, 
            pf, self.options["md"]["pdamp"]))
        lmp.command("fix_modify        f1 temp tcm")
        lmp.command("fix               f3 all adapt 1 pair %s scale * * v_lambda"%self.options["md"]["pair_style"])
        lmp.command("fix               f4 all print 1 \"${dU} $(press) $(vol) ${lambda}\" screen no file forward_%d.dat"%iteration)
        lmp.command("run               ${ts}")
        lmp.command("unfix             f3")
        lmp.command("unfix             f4")
        lmp.command("unfix             f1")

        #check for melting
        lmp.command("dump              2 all custom 1 traj.dat id type mass x y z vx vy vz")
        lmp.command("run               0")
        lmp.command("undump            2")
        
        solids = ph.find_solid_fraction("traj.dat")
        if (solids/lmp.natoms < 0.7):
            lmp.close()
            raise RuntimeError("System melted, increase size or reduce scaling!")
        
        lmp.command("fix              f1 all nph iso %f %f %f fixedpoint ${xcm} ${ycm} ${zcm}"%(pf, 
            pf, self.options["md"]["pdamp"]))
        lmp.command("fix_modify        f1 temp tcm")
        lmp.command("run               ${te}")
        lmp.command("unfix             f1")

        lmp.command("variable          lambda equal ramp(${lf},${li})")
        
        lmp.command("fix              f1 all nph iso %f %f %f fixedpoint ${xcm} ${ycm} ${zcm}"%(pf, 
            pi, self.options["md"]["pdamp"]))
        lmp.command("fix_modify        f1 temp tcm")
        lmp.command("fix               f3 all adapt 1 pair %s scale * * v_lambda"%self.options["md"]["pair_style"])
        lmp.command("fix               f4 all print 1 \"${dU} $(press) $(vol) ${lambda}\" screen no file backward_%d.dat"%iteration)
        lmp.command("run               ${ts}")
        lmp.command("unfix             f3")
        lmp.command("unfix             f4")
        lmp.command("unfix             f1")
        
        lmp.close()


    def integrate_reversible_scaling(self, scale_energy=False, return_values=False):
        """
        Carry out the reversible scaling operation
        """
        integrate_rs(self.simfolder, self.fe, self.t, self.natoms, p=self.p,
            nsims=self.nsims, scale_energy=scale_energy, return_values=return_values)
