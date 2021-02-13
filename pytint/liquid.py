"""
Contains methods for liquid
"""

import numpy as np
import yaml

from pytint.integrators import *
import pyscal.traj_process as ptp
import pytint.lattice as pl
import pytint.helpers as ph



class Liquid:
    """
    Liquid class

    Parameters
    ----------
    t : float
        simulation temperature

    p : float
        pressure

    l : int
        selected lattice indicator

    apc : int
        atoms per cell

    alat : float
        lattice constant

    c : float
        concentration

    options: options class
        input options

    simfolder: string
        simulation folder

    thigh : float
        temperature to melt the structure 
    """
    def __init__(self, options=None, kernel=None, simfolder=None):
        """
        Set up class
        """
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
        self.prepare_lattice()

        logfile = os.path.join(self.simfolder, "tint.log")
        self.logger = ph.prepare_log(logfile)

        #other properties
        self.cores = self.options["queue"]["cores"]
        self.ncells = np.prod(self.calc["repeat"])
        
        if self.l == "file":
            self.natoms = ph.check_data_file(self.calc["lattice"])
            #reset apc
            self.apc = self.natoms
        else:
            self.natoms = self.ncells*self.apc
        
        #the UFM system properties
        self.eps = self.t*50.0*kb

        #properties that will be calculated later
        self.volatom = None
        self.rho = None
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
        l, alat, apc = pl.prepare_lattice(self.calc)
        self.l = l
        self.alat = alat
        self.apc = apc


    def run_averaging(self):
        """
        Run averaging cycle

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        Run the averaging cycle to find the equilibrium number
        density at the given temperature.
        """
        
        #create lammps object
        lmp = ph.create_object(self.cores, self.simfolder, self.options["md"]["timestep"])

        #set up structure
        lmp = ph.create_structure(lmp, self.calc)

        #set up potential
        lmp = ph.set_potential(lmp, self.options)

        #Melt regime for the liquid
        lmp.velocity("all create", self.thigh, np.random.randint(0, 10000))
        
        #add some computes
        lmp.command("variable         mvol equal vol")
        lmp.command("variable         mlx equal lx")
        lmp.command("variable         mly equal ly")
        lmp.command("variable         mlz equal lz")
        lmp.command("variable         mpress equal press") 

        #melting cycle
        # try with different multiples of thmult until the structure melts.
        melted = False
        for thmult in np.arange(1.0, 2.0, 0.1):
            
            trajfile = os.path.join(self.simfolder, "traj.melt")
            if os.path.exists(trajfile):
                os.remove(trajfile)

            self.logger.info("Starting melting cycle with thigh temp %f, factor %f"%(self.thigh, thmult))
            lmp.fix("1 all npt temp", self.thigh*thmult, self.thigh*thmult, self.options["md"]["tdamp"], 
                                          "iso", self.p, self.p, self.options["md"]["pdamp"])
            lmp.run(int(self.options["md"]["nsmall"]))
            lmp.unfix(1)
            lmp.dump("2 all custom", 1, trajfile,"id type mass x y z vx vy vz")
            lmp.run(0)
            lmp.undump(2)
            
            #we have to check if the structure melted
            solids = ph.find_solid_fraction("traj.melt")
            self.logger.info("fraction of solids found: %f", solids/self.natoms)
            if (solids/self.natoms < self.options["conv"]["liquid_frac"]):
                melted = True
                break
        
        #if not melted throw error
        if not melted:
            raise ValueError("Liquid system did not melt, maybe try a higher thigh temperature.")

        #now assign correct temperature and equilibrate
        lmp.velocity("all create", self.t, np.random.randint(0, 10000))
        lmp.fix("1 all npt temp", self.t, self.t, self.options["md"]["tdamp"], 
                                      "iso", self.p, self.p, self.options["md"]["pdamp"])
        lmp.run(int(self.options["md"]["nsmall"])) 
        
        #start recording average values
        lmp.command("fix              2 all ave/time %d %d %d v_mlx v_mly v_mlz v_mpress file avg.dat"%(int(self.options["md"]["nevery"]),
            int(self.options["md"]["nrepeat"]), int(self.options["md"]["nevery"]*self.options["md"]["nrepeat"])))

        #monitor the average values until calculation is converged
        laststd = 0.00
        for i in range(self.options["md"]["ncycles"]):
            lmp.command("run              %d"%int(self.options["md"]["nsmall"]))
            ncount = int(self.options["md"]["nsmall"])//int(self.options["md"]["nevery"]*self.options["md"]["nrepeat"])
            #now we can check if it converted
            file = os.path.join(self.simfolder, "avg.dat")
            lx, ly, lz, ipress = np.loadtxt(file, usecols=(1,2,3,4), unpack=True)

            #lxpc = ((lx*ly*lz)/self.ncells)**(1/3)
            #lxpc = lxpc[-ncount+1:]
            lxpc = ipress
            mean = np.mean(lxpc)
            std = np.std(lxpc)
            volatom = np.mean((lx*ly*lz)/natoms)            
            self.logger.info("At count %d mean pressure is %f with vol/atom %f"%(i+1, mean, volatom))

            #if (np.abs(laststd - std) < self.options["conv"]["alat_tol"]):
            if (np.abs(mean - self.p)) < self.options["conv"]["p_tol"]:
                #self.avglat = np.round(mean, decimals=3)

                #process other means
                self.lx = np.round(np.mean(lx[-ncount+1:]), decimals=3)
                self.ly = np.round(np.mean(ly[-ncount+1:]), decimals=3)
                self.lz = np.round(np.mean(lz[-ncount+1:]), decimals=3)
                self.volatom = volatom
                self.vol = self.lx*self.ly*self.lz
                self.rho = self.natoms/(self.lx*self.ly*self.lz)

                self.logger.info("finalized vol/atom %f at pressure %f"%(self.volatom, mean))
                self.logger.info("Avg box dimensions x: %f, y: %f, z:%f"%(self.lx, self.ly, self.lz))
                break
            
            laststd = std

        lmp.command("dump              2 all custom 1 traj.dat id type mass x y z vx vy vz")
        lmp.command("run               0")
        lmp.command("undump            2")

        #finish run
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
        Write TI integrate script
        """

        #create lammps object
        lmp = ph.create_object(self.cores, self.simfolder, self.options["md"]["timestep"])

        
        lmp.command("variable        rnd      equal   round(random(0,999999,%d))"%np.random.randint(0, 10000))
        lmp.command("variable        dt       equal   %f"%self.options["md"]["timestep"])             # Timestep (ps).
        # Adiabatic switching parameters.
        lmp.command("variable        li       equal   1.0")               # Initial lambda.
        lmp.command("variable        lf       equal   0.0")               # Final lambda.
        #------------------------------------------------------------------------------------------------------#

        conf = os.path.join(self.simfolder, "conf.dump")
        lmp = ph.read_dump(lmp, conf, species=self.options["nelements"])

        # Define MEAM and UF potentials parameters.
        lmp = ph.set_hybrid_potential(lmp, self.options, self.eps)

        #remap the box to get the correct pressure
        lmp = ph.remap_box(lmp, self.lx, self.ly, self.lz)

        ################################     Fixes, computes and constraints     ###############################
        # Integrator & thermostat.
        lmp.command("fix             f1 all nve")                              
        lmp.command("fix             f2 all langevin %f %f %f ${rnd}"%(self.t, self.t, self.options["md"]["tdamp"]))
        lmp.command("variable        rnd equal round(random(0,999999,0))")

        # Compute the potential energy of each pair style.
        lmp.command("compute         c1 all pair %s"%self.options["md"]["pair_style"])
        lmp.command("compute         c2 all pair ufm")
        #------------------------------------------------------------------------------------------------------#


        ##########################################     Output setup     ########################################
        # Output variables.
        lmp.command("variable        step equal step")
        lmp.command("variable        dU1 equal c_c1/atoms")             # Driving-force obtained from NEHI procedure.
        lmp.command("variable        dU2 equal c_c2/atoms")

        # Thermo output.
        lmp.command("thermo_style    custom step v_dU1 v_dU2")
        lmp.command("thermo          1000")
        #------------------------------------------------------------------------------------------------------#


        ##########################################     Run simulation     ######################################
        # Turn UF potential off (completely) to equilibrate the Sw potential.
        lmp.command("variable        zero equal 0")
        lmp.command("fix             f0 all adapt 0 pair ufm scale * * v_zero")
        lmp.command("run             0")
        lmp.command("unfix           f0")

        # Equilibrate the fluid interacting by Sw potential and switch to UF potential (Forward realization).
        lmp.command("run             %d"%self.options["md"]["te"])

        lmp.command("print           \"${dU1} ${dU2} ${li}\" file forward_%d.dat"%iteration)
        lmp.command("variable        lambda_sw equal ramp(${li},${lf})")                 # Linear lambda protocol from 1 to 0.
        lmp.command("fix             f3 all adapt 1 pair %s scale * * v_lambda_sw"%self.options["md"]["pair_style"])
        lmp.command("variable        lambda_ufm equal ramp(${lf},${li})")                  # Linear lambda protocol from 0 to 1.
        lmp.command("fix             f4 all adapt 1 pair ufm scale * * v_lambda_ufm")
        lmp.command("fix             f5 all print 1 \"${dU1} ${dU2} ${lambda_sw}\" screen no append forward_%d.dat"%iteration)
        lmp.command("run             %d"%self.options["md"]["ts"])

        lmp.command("unfix           f3")
        lmp.command("unfix           f4")
        lmp.command("unfix           f5")

        # Equilibrate the fluid interacting by UF potential and switch to sw potential (Backward realization).
        lmp.command("run             %d"%self.options["md"]["te"])

        lmp.command("print           \"${dU1} ${dU2} ${lf}\" file backward_%d.dat"%iteration)
        lmp.command("variable        lambda_sw equal ramp(${lf},${li})")                 # Linear lambda protocol from 0 to 1.
        lmp.command("fix             f3 all adapt 1 pair %s scale * * v_lambda_sw"%self.options["md"]["pair_style"])
        lmp.command("variable        lambda_ufm equal ramp(${li},${lf})")                  # Linear lambda protocol from 1 to 0.
        lmp.command("fix             f4 all adapt 1 pair ufm scale * * v_lambda_ufm")
        lmp.command("fix             f5 all print 1 \"${dU1} ${dU2} ${lambda_sw}\" screen no append backward_%d.dat"%iteration)
        lmp.command("run             %d"%self.options["md"]["ts"])

        lmp.command("unfix           f3")
        lmp.command("unfix           f4")
        lmp.command("unfix           f5")
        #------------------------------------------------------------------------------------------------------#
        lmp.close()
    
    def thermodynamic_integration(self):
        """
        Perform thermodynamic integration
        """
        w, q, qerr = find_w(self.simfolder, nsims=self.nsims, 
            full=True, solid=False)  
        #WARNING: hardcoded UFM parameters           
        f1 = get_uhlenbeck_ford_fe(self.t, 
            self.rho, 50, 1.5)

        #we need to find concentration too to find ideal gas fe for multi species
        f2 = get_ideal_gas_fe(self.t, self.rho, 
            self.natoms, self.options["mass"], self.calc["concentration"])
        
        self.ferr = qerr
        self.fref = f1
        self.fideal = f2
        self.w = w

        if self.p != 0:
            #add pressure contribution
            p = self.p/(10000*160.21766208)
            v = self.vol/self.natoms
            self.pv = p*v
        else:
            self.pv = 0

        self.fe = self.fideal + self.fref - self.w + self.pv


    def submit_report(self):

        report = {}

        #input quantities
        report["input"] = {}
        report["input"]["temperature"] = int(self.t)
        report["input"]["pressure"] = float(self.p)
        report["input"]["lattice"] = str(self.l)
        report["input"]["element"] = " ".join(np.array(self.options["element"]).astype(str))
        report["input"]["concentration"] = " ".join(np.array(self.calc["concentration"]).astype(str))

        #average quantities
        report["average"] = {}
        report["average"]["vol/atom"] = float(self.volatom)
        report["average"]["density"] = float(self.rho)
        
        #results
        report["results"] = {}
        report["results"]["free_energy"] = float(self.fe)
        report["results"]["error"] = float(self.ferr)
        report["results"]["reference_system"] = float(self.fref) + float(self.fideal)
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
        
        #first we need to do an averaging scheme
        #but only if iteration is 1
        
        #Now do reversible scaling
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
        lmp.command("fix               f1 all nph iso %f %f %f"%(self.p, self.p, self.options["md"]["pdamp"]))
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
