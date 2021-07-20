"""
calphy: a Python library and command line interface for automated free
energy calculations.

Copyright 2021  (c) Sarath Menon^1, Yury Lysogorskiy^1, Ralf Drautz^1
^1: Ruhr-University Bochum, Bochum, Germany

More information about the program can be found in:
Menon, Sarath, Yury Lysogorskiy, Jutta Rogal, and Ralf Drautz. 
“Automated Free Energy Calculation from Atomistic Simulations.” 
ArXiv:2107.08980 [Cond-Mat], July 19, 2021. 
http://arxiv.org/abs/2107.08980.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

See the LICENSE file.

For more information contact:
sarath.menon@ruhr-uni-bochum.de
"""

import numpy as np
import yaml

from calphy.integrators import *
import pyscal.traj_process as ptp
import calphy.lattice as pl
import calphy.helpers as ph

class Alchemy:
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
        
        if self.calc["iso"]:
            self.iso = "iso"
        else:
            self.iso = "aniso"


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

        #now backup pair styles: to use in integration mode
        self.pair_style = self.options["md"]["pair_style"]
        self.pair_coeff = self.options["md"]["pair_coeff"]

        #now manually tune pair styles
        self.options["md"]["pair_style"] = self.options["md"]["pair_style"][0]
        self.options["md"]["pair_coeff"] = self.options["md"]["pair_coeff"][0]


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
            lmp.command("fix              1 all npt temp %f %f %f %s %f %f %f"%(self.t, self.t, self.options["md"]["tdamp"], 
                                                self.iso, self.p, self.p, self.options["md"]["pdamp"]))
            lmp.command("thermo_style     custom step pe press vol etotal temp lx ly lz")
            lmp.command("thermo           10")
            lmp.command("run              %d"%int(self.options["md"]["nsmall"])) 

        else:
            #Now this routine is for non-zero pressure
            #one has to equilibriate at a low temperature, but high pressure and then increase temp gradually
            #start at 0.25 temp, and increase to 0.50, while keeping high pressure
            lmp.command("velocity         all create %f %d"%(0.25*self.t, np.random.randint(0, 10000)))
            lmp.command("fix              1 all npt temp %f %f %f %s %f %f %f"%(0.25*self.t, 0.5*self.t, self.options["md"]["tdamp"], 
                                                self.iso, self.p, self.p, self.options["md"]["pdamp"]))
            lmp.command("thermo_style     custom step pe press vol etotal temp")
            lmp.command("thermo           10")
            lmp.command("run              %d"%int(self.options["md"]["nsmall"])) 
            lmp.command("unfix            1")

            #now heat again
            lmp.command("fix              1 all npt temp %f %f %f %s %f %f %f"%(0.5*self.t, self.t, self.options["md"]["tdamp"], 
                                                self.iso, self.p, self.p,  self.options["md"]["pdamp"]))
            lmp.command("run              %d"%int(self.options["md"]["nsmall"])) 
            lmp.command("unfix            1")

            #now run normal cycle
            lmp.command("fix              1 all npt temp %f %f %f %s %f %f %f"%(self.t, self.t, self.options["md"]["tdamp"], 
                                                self.iso, self.p, self.p,  self.options["md"]["pdamp"]))
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

        lmp = ph.set_double_hybrid_potential(lmp, self.options, self.pair_style, self.pair_coeff)

        #remap the box to get the correct pressure
        lmp = ph.remap_box(lmp, self.lx, self.ly, self.lz)

        ################################     Fixes, computes and constraints     ###############################
        # Integrator & thermostat.
        lmp.command("fix             f1 all nve")                              
        lmp.command("fix             f2 all langevin %f %f %f ${rnd}"%(self.t, self.t, self.options["md"]["tdamp"]))
        lmp.command("variable        rnd equal round(random(0,999999,0))")

        # Compute pair definitions
        #lmp.command("include  %s"%self.calc["compute_file"])

        if self.pair_style[0] == self.pair_style[1]:
            lmp.command("compute         c1 all pair %s 1"%self.pair_style[0])
            lmp.command("compute         c2 all pair %s 2"%self.pair_style[1])
        else:
            lmp.command("compute         c1 all pair %s"%self.pair_style[0])
            lmp.command("compute         c2 all pair %s"%self.pair_style[1])


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

        if self.pair_style[0] == self.pair_style[1]:
            lmp.command("fix             f0 all adapt 0 pair %s:2 scale * * v_zero"%self.pair_style[1])
        else:
            lmp.command("fix             f0 all adapt 0 pair %s scale * * v_zero"%self.pair_style[1])
        lmp.command("run             0")
        lmp.command("unfix           f0")


        # Equilibrate the fluid interacting by Sw potential and switch to UF potential (Forward realization).
        lmp.command("run             %d"%self.options["md"]["te"])

        lmp.command("print           \"${dU1} ${dU2} ${li}\" file forward_%d.dat"%iteration)
        lmp.command("variable        lambda_sw equal ramp(${li},${lf})")                 # Linear lambda protocol from 1 to 0.

        if self.pair_style[0] == self.pair_style[1]:
            lmp.command("fix             f3 all adapt 1 pair %s:1 scale * * v_lambda_sw"%self.pair_style[0])
        else:
            lmp.command("fix             f3 all adapt 1 pair %s scale * * v_lambda_sw"%self.pair_style[0])        
        
        lmp.command("variable        lambda_ufm equal ramp(${lf},${li})")                  # Linear lambda protocol from 0 to 1.
        
        if self.pair_style[0] == self.pair_style[1]:
            lmp.command("fix             f4 all adapt 1 pair %s:2 scale * * v_lambda_ufm"%self.pair_style[1])
        else:
            lmp.command("fix             f4 all adapt 1 pair %s scale * * v_lambda_ufm"%self.pair_style[1])
        
        lmp.command("fix             f5 all print 1 \"${dU1} ${dU2} ${lambda_sw}\" screen no append forward_%d.dat"%iteration)
        lmp.command("run             %d"%self.options["md"]["ts"])

        lmp.command("unfix           f3")
        lmp.command("unfix           f4")
        lmp.command("unfix           f5")

        # Equilibrate the fluid interacting by UF potential and switch to sw potential (Backward realization).
        lmp.command("run             %d"%self.options["md"]["te"])

        lmp.command("print           \"${dU1} ${dU2} ${lf}\" file backward_%d.dat"%iteration)
        lmp.command("variable        lambda_sw equal ramp(${lf},${li})")                 # Linear lambda protocol from 0 to 1.
        
        if self.pair_style[0] == self.pair_style[1]:
            lmp.command("fix             f3 all adapt 1 pair %s:1 scale * * v_lambda_sw"%self.pair_style[0])
        else:
            lmp.command("fix             f3 all adapt 1 pair %s scale * * v_lambda_sw"%self.pair_style[0])        
        
        lmp.command("variable        lambda_ufm equal ramp(${li},${lf})")                  # Linear lambda protocol from 1 to 0.
        
        if self.pair_style[0] == self.pair_style[1]:
            lmp.command("fix             f4 all adapt 1 pair %s:2 scale * * v_lambda_ufm"%self.pair_style[1])
        else:
            lmp.command("fix             f4 all adapt 1 pair %s scale * * v_lambda_ufm"%self.pair_style[1])
        
        lmp.command("fix             f5 all print 1 \"${dU1} ${dU2} ${lambda_sw}\" screen no append backward_%d.dat"%iteration)
        lmp.command("run             %d"%self.options["md"]["ts"])

        lmp.command("unfix           f3")
        lmp.command("unfix           f4")
        lmp.command("unfix           f5")
        #------------------------------------------------------------------------------------------------------#
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
        w, q, qerr = find_w(self.simfolder, nelements=self.options["nelements"], concentration=self.concentration, nsims=self.nsims, 
            full=True, solid=False, alchemy=True)
        
        self.w = w
        self.ferr = qerr

        self.fe = self.w


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
        
        #results
        report["results"] = {}
        report["results"]["free_energy"] = float(self.fe)
        report["results"]["error"] = float(self.ferr)
        report["results"]["work"] = float(self.w)

        reportfile = os.path.join(self.simfolder, "report.yaml")
        with open(reportfile, 'w') as f:
            yaml.dump(report, f)
