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
import calphy.phase as cph

class Liquid(cph.Phase):
    """
    Class for free energy calculation with liquid as the reference state

    Parameters
    ----------
    options : dict
        dict of input options
    
    kernel : int
        the index of the calculation that should be run from
        the list of calculations in the input file

    simfolder : string
        base folder for running calculations

    """
    def __init__(self, options=None, kernel=None, simfolder=None):
        """
        Set up class
        """
        #call base class
        super().__init__(options=options,
        kernel=kernel, simfolder=simfolder)




    def run_averaging(self):
        """
        Run averaging routine

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        Run averaging routine using LAMMPS. Starting from the initial lattice two different routines can
        be followed:
        If pressure is specified, MD simulations are run until the pressure converges within the given
        threshold value.
        If `fix_lattice` option is True, then the input structure is used as it is and the corresponding pressure
        is calculated.
        At the end of the run, the averaged box dimensions are calculated. 
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
        #try with different multiples of thmult until the structure melts.
        #TODO : Add option to skip melting cycle
        
        melted = False
        
        #this is the multiplier for thigh to try melting routines
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
        
        #if melting cycle is over and still not melted, raise error
        if not melted:
            lmp.close()
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
        converged = False

        for i in range(self.options["md"]["ncycles"]):
            lmp.command("run              %d"%int(self.options["md"]["nsmall"]))
            ncount = int(self.options["md"]["nsmall"])//int(self.options["md"]["nevery"]*self.options["md"]["nrepeat"])
            #now we can check if it converted
            file = os.path.join(self.simfolder, "avg.dat")
            lx, ly, lz, ipress = np.loadtxt(file, usecols=(1,2,3,4), unpack=True)

            lxpc = ipress
            mean = np.mean(lxpc)
            std = np.std(lxpc)
            volatom = np.mean((lx*ly*lz)/self.natoms)            
            self.logger.info("At count %d mean pressure is %f with vol/atom %f"%(i+1, mean, volatom))

            if (np.abs(mean - self.p)) < self.options["conv"]["p_tol"]:
                #process other means
                self.lx = np.round(np.mean(lx[-ncount+1:]), decimals=3)
                self.ly = np.round(np.mean(ly[-ncount+1:]), decimals=3)
                self.lz = np.round(np.mean(lz[-ncount+1:]), decimals=3)
                self.volatom = volatom
                self.vol = self.lx*self.ly*self.lz
                self.rho = self.natoms/(self.lx*self.ly*self.lz)

                self.logger.info("finalized vol/atom %f at pressure %f"%(self.volatom, mean))
                self.logger.info("Avg box dimensions x: %f, y: %f, z:%f"%(self.lx, self.ly, self.lz))
                converged = True
                break
            
            laststd = std

        if not converged:
            lmp.close()
            raise ValueError("Pressure did not converge after MD runs, maybe change lattice_constant and try?")

        lmp.command("dump              2 all custom 1 traj.dat id type mass x y z vx vy vz")
        lmp.command("run               0")
        lmp.command("undump            2")

        #finish run and close object
        lmp.close()

        #process the trajectory
        self.process_traj()




    def run_integration(self, iteration=1):
        """
        Run integration routine

        Parameters
        ----------
        iteration : int, optional
            iteration number for running independent iterations

        Returns
        -------
        None

        Notes
        -----
        Run the integration routine where the initial and final systems are connected using
        the lambda parameter. See algorithm 4 in publication.
        """
        lmp = ph.create_object(self.cores, self.simfolder, self.options["md"]["timestep"])

        # Adiabatic switching parameters.
        lmp.command("variable        li       equal   1.0")
        lmp.command("variable        lf       equal   0.0")

        #read in the conf file
        conf = os.path.join(self.simfolder, "conf.dump")
        lmp = ph.read_dump(lmp, conf, species=self.options["nelements"])

        #set hybrid ufm and normal potential
        lmp = ph.set_hybrid_potential(lmp, self.options, self.eps)

        #remap the box to get the correct pressure
        lmp = ph.remap_box(lmp, self.lx, self.ly, self.lz)

        #apply the necessary thermostat
        lmp.command("fix             f1 all nve")                              
        lmp.command("fix             f2 all langevin %f %f %f %d"%(self.t, self.t, self.options["md"]["tdamp"],
            np.random.randint(0, 10000)))

        # Compute the potential energy of each pair style.
        lmp.command("compute         c1 all pair %s"%self.options["md"]["pair_style"])
        lmp.command("compute         c2 all pair ufm")

        # Output variables.
        lmp.command("variable        step equal step")
        lmp.command("variable        dU1 equal c_c1/atoms")             # Driving-force obtained from NEHI procedure.
        lmp.command("variable        dU2 equal c_c2/atoms")

        #force thermo to evaluate variables
        lmp.command("thermo_style    custom step v_dU1 v_dU2")
        lmp.command("thermo          1000")

        #switching completely to potential of interest
        lmp.command("variable        zero equal 0")
        lmp.command("fix             f0 all adapt 0 pair ufm scale * * v_zero")
        lmp.command("run             0")
        lmp.command("unfix           f0")

        #Equilibrate system
        lmp.command("run             %d"%self.options["md"]["te"])

        #print header
        lmp.command("print           \"${dU1} ${dU2} ${li}\" file forward_%d.dat"%iteration)
        
        #set up scaling variables
        lmp.command("variable        lambda_p1 equal ramp(${li},${lf})")
        lmp.command("variable        lambda_p2 equal ramp(${lf},${li})")

        #Forward switching run
        lmp.command("fix             f3 all adapt 1 pair %s scale * * v_lambda_p1"%self.options["md"]["pair_style"])
        lmp.command("fix             f4 all adapt 1 pair ufm scale * * v_lambda_p2")
        lmp.command("fix             f5 all print 1 \"${dU1} ${dU2} ${lambda_p1}\" screen no append forward_%d.dat"%iteration)
        lmp.command("run             %d"%self.options["md"]["ts"])

        #unfix things
        lmp.command("unfix           f3")
        lmp.command("unfix           f4")
        lmp.command("unfix           f5")

        #Equilibriate at UFM potential
        lmp.command("run             %d"%self.options["md"]["te"])

        #print file header
        lmp.command("print           \"${dU1} ${dU2} ${lf}\" file backward_%d.dat"%iteration)
        
        #set up scaling variables
        lmp.command("variable        lambda_p1 equal ramp(${lf},${li})")
        lmp.command("variable        lambda_p2 equal ramp(${li},${lf})")

        #Reverse switching run
        lmp.command("fix             f3 all adapt 1 pair %s scale * * v_lambda_p1"%self.options["md"]["pair_style"])
        lmp.command("fix             f4 all adapt 1 pair ufm scale * * v_lambda_p2")
        lmp.command("fix             f5 all print 1 \"${dU1} ${dU2} ${lambda_p1}\" screen no append backward_%d.dat"%iteration)
        lmp.command("run             %d"%self.options["md"]["ts"])

        #unfix things
        lmp.command("unfix           f3")
        lmp.command("unfix           f4")
        lmp.command("unfix           f5")
        
        #close object
        lmp.close()
    
    def thermodynamic_integration(self):
        """
        Calculate free energy after integration step

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        Calculates the final work, energy dissipation and free energy by
        matching with UFM model
        """
        w, q, qerr = find_w(self.simfolder, nsims=self.nsims, 
            full=True, solid=False)  
        
        #TODO: Hardcoded UFM parameters - enable option to change          
        f1 = get_uhlenbeck_ford_fe(self.t, 
            self.rho, 50, 1.5)

        #Get ideal gas fe
        f2 = get_ideal_gas_fe(self.t, self.rho, 
            self.natoms, self.options["mass"], self.concentration)
        
        self.ferr = qerr
        self.fref = f1
        self.fideal = f2
        self.w = w

        #add pressure contribution if required
        if self.p != 0:
            p = self.p/(10000*160.21766208)
            v = self.vol/self.natoms
            self.pv = p*v
        else:
            self.pv = 0

        #calculate final free energy
        self.fe = self.fideal + self.fref - self.w + self.pv



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
        t0 = self.t
        tf = self.tend
        li = 1
        lf = t0/tf
        pi = self.p
        pf = lf*pi

        #create lammps object
        lmp = ph.create_object(self.cores, self.simfolder, self.options["md"]["timestep"])

        lmp.command("echo              log")
        lmp.command("variable          li equal %f"%li)
        lmp.command("variable          lf equal %f"%lf)

        #read in conf file
        conf = os.path.join(self.simfolder, "conf.dump")
        lmp = ph.read_dump(lmp, conf, species=self.options["nelements"])

        #set up potential
        lmp = ph.set_potential(lmp, self.options)

        #remap the box to get the correct pressure
        lmp = ph.remap_box(lmp, self.lx, self.ly, self.lz)

        #set thermostat and run equilibrium
        lmp.command("fix               f1 all nph iso %f %f %f"%(self.p, self.p, self.options["md"]["pdamp"]))
        lmp.command("fix               f2 all langevin %f %f %f %d zero yes"%(self.t, self.t, self.options["md"]["tdamp"], np.random.randint(0, 10000)))
        lmp.command("run               %d"%self.options["md"]["te"])
        lmp.command("unfix             f1")
        lmp.command("unfix             f2")

        #now fix com
        lmp.command("variable         xcm equal xcm(all,x)")
        lmp.command("variable         ycm equal xcm(all,y)")
        lmp.command("variable         zcm equal xcm(all,z)")
        
        lmp.command("fix              f1 all nph iso %f %f %f fixedpoint ${xcm} ${ycm} ${zcm}"%(self.p, self.p, self.options["md"]["pdamp"]))
        lmp.command("fix              f2 all langevin %f %f %f %d zero yes"%(t0, t0, self.options["md"]["tdamp"], np.random.randint(0, 10000)))
        
        #compute com and modify fix
        lmp.command("compute           tcm all temp/com")
        lmp.command("fix_modify        f1 temp tcm")
        lmp.command("fix_modify        f2 temp tcm")

        lmp.command("variable          step    equal step")
        lmp.command("variable          dU      equal c_thermo_pe/atoms")        
        lmp.command("thermo_style      custom step pe c_tcm press vol")
        lmp.command("thermo            10000")

        #create velocity and equilibriate
        lmp.command("velocity          all create %f %d mom yes rot yes dist gaussian"%(t0, np.random.randint(0, 10000)))   
        lmp.command("run               %d"%self.options["md"]["te"])
        
        #unfix nph
        lmp.command("unfix             f1")

        #define ramp
        lmp.command("variable          lambda equal ramp(${li},${lf})")

        #start scaling over switching time
        lmp.command("fix              f1 all nph iso %f %f %f fixedpoint ${xcm} ${ycm} ${zcm}"%(pi, 
            pf, self.options["md"]["pdamp"]))
        lmp.command("fix_modify        f1 temp tcm")
        lmp.command("fix               f3 all adapt 1 pair %s scale * * v_lambda"%self.options["md"]["pair_style"])
        lmp.command("fix               f4 all print 1 \"${dU} $(press) $(vol) ${lambda}\" screen no file forward_%d.dat"%iteration)
        lmp.command("run               %d"%self.options["md"]["ts"])

        #unfix
        lmp.command("unfix             f3")
        lmp.command("unfix             f4")
        lmp.command("unfix             f1")

         #equilibriate scaled hamiltonian
        lmp.command("fix              f1 all nph iso %f %f %f fixedpoint ${xcm} ${ycm} ${zcm}"%(pf, 
            pf, self.options["md"]["pdamp"]))
        lmp.command("fix_modify        f1 temp tcm")        
        lmp.command("run               %d"%self.options["md"]["te"])
        lmp.command("unfix             f1")

        #reverse lambda ramp
        lmp.command("variable          lambda equal ramp(${lf},${li})")

        #apply fix and perform switching        
        lmp.command("fix              f1 all nph iso %f %f %f fixedpoint ${xcm} ${ycm} ${zcm}"%(pf, 
            pi, self.options["md"]["pdamp"]))
        lmp.command("fix_modify        f1 temp tcm")
        lmp.command("fix               f3 all adapt 1 pair %s scale * * v_lambda"%self.options["md"]["pair_style"])
        lmp.command("fix               f4 all print 1 \"${dU} $(press) $(vol) ${lambda}\" screen no file backward_%d.dat"%iteration)
        lmp.command("run               %d"%self.options["md"]["ts"])
        lmp.command("unfix             f3")
        lmp.command("unfix             f4")
        lmp.command("unfix             f1")
        
        #close the object
        lmp.close()

