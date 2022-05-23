"""
calphy: a Python library and command line interface for automated free
energy calculations.

Copyright 2021  (c) Sarath Menon^1, Yury Lysogorskiy^2, Ralf Drautz^2
^1: Max Planck Institut für Eisenforschung, Dusseldorf, Germany 
^2: Ruhr-University Bochum, Bochum, Germany

calphy is published and distributed under the Academic Software License v1.0 (ASL). 
calphy is distributed in the hope that it will be useful for non-commercial academic research, 
but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
calphy API is published and distributed under the BSD 3-Clause "New" or "Revised" License
See the LICENSE FILE for more details. 

More information about the program can be found in:
Menon, Sarath, Yury Lysogorskiy, Jutta Rogal, and Ralf Drautz.
“Automated Free Energy Calculation from Atomistic Simulations.” Physical Review Materials 5(10), 2021
DOI: 10.1103/PhysRevMaterials.5.103801

For more information contact:
sarath.menon@ruhr-uni-bochum.de/yury.lysogorskiy@icams.rub.de
"""

import numpy as np
import yaml

import pyscal.traj_process as ptp
from calphy.integrators import *
import calphy.lattice as pl
import calphy.helpers as ph
import calphy.phase as cph

class Alchemy(cph.Phase):
    """
    Class for alchemical transformations

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
    def __init__(self, calculation=None, simfolder=None):

        #call base class
        super().__init__(calculation=calculation,
        simfolder=simfolder)


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
        Fix lattice option is not implemented at present.
        At the end of the run, the averaged box dimensions are calculated. 
        """
        lmp = ph.create_object(self.cores, self.simfolder, self.calc.md.timestep)

        #set up structure
        lmp = ph.create_structure(lmp, self.calc)

        #set up potential
        lmp = ph.set_potential(lmp, self.calc)

        #add some computes
        lmp.command("variable         mvol equal vol")
        lmp.command("variable         mlx equal lx")
        lmp.command("variable         mly equal ly")
        lmp.command("variable         mlz equal lz")
        lmp.command("variable         mpress equal press")

        if not self.calc._fix_lattice:
            if self.calc._pressure == 0:
                #This routine should be followed for zero pressure
                lmp.command("velocity         all create %f %d"%(self.calc._temperature, np.random.randint(0, 10000)))
                lmp.command("fix              1 all npt temp %f %f %f %s %f %f %f"%(self.calc._temperature, self.calc._temperature, self.calc.md.thermostat_damping, 
                                                    self.iso, self.calc._pressure, self.calc._pressure, self.calc.md.barostat_damping))
                lmp.command("thermo_style     custom step pe press vol etotal temp lx ly lz")
                lmp.command("thermo           10")
                lmp.command("run              %d"%int(self.calc.md.n_small_steps)) 

            else:
                #Now this routine is for non-zero pressure
                #one has to equilibriate at a low temperature, but high pressure and then increase temp gradually
                #start at 0.25 temp, and increase to 0.50, while keeping high pressure
                lmp.command("velocity         all create %f %d"%(0.25*self.calc._temperature, np.random.randint(0, 10000)))
                lmp.command("fix              1 all npt temp %f %f %f %s %f %f %f"%(0.25*self.calc._temperature, 0.5*self.calc._temperature, self.calc.md.thermostat_damping, 
                                                    self.iso, self.calc._pressure, self.calc._pressure, self.calc.md.barostat_damping))
                lmp.command("thermo_style     custom step pe press vol etotal temp")
                lmp.command("thermo           10")
                lmp.command("run              %d"%int(self.calc.md.n_small_steps)) 
                lmp.command("unfix            1")

                #now heat again
                lmp.command("fix              1 all npt temp %f %f %f %s %f %f %f"%(0.5*self.calc._temperature, self.calc._temperature, self.calc.md.thermostat_damping, 
                                                    self.iso, self.calc._pressure, self.calc._pressure,  self.calc.md.barostat_damping))
                lmp.command("run              %d"%int(self.calc.md.n_small_steps)) 
                lmp.command("unfix            1")

                #now run normal cycle
                lmp.command("fix              1 all npt temp %f %f %f %s %f %f %f"%(self.calc._temperature, self.calc._temperature, self.calc.md.thermostat_damping, 
                                                    self.iso, self.calc._pressure, self.calc._pressure,  self.calc.md.barostat_damping))
                lmp.command("run              %d"%int(self.calc.md.n_small_steps)) 


            #this is when the averaging routine starts
            lmp.command("fix              2 all ave/time %d %d %d v_mlx v_mly v_mlz v_mpress file avg.dat"%(int(self.calc.md.n_every_steps),
                int(self.calc.md.n_repeat_steps), int(self.calc.md.n_every_steps*self.calc.md.n_repeat_steps)))
            
            laststd = 0.00
            converged = False

            for i in range(int(self.calc.md.n_cycles)):
                lmp.command("run              %d"%int(self.calc.md.n_small_steps))
                ncount = int(self.calc.md.n_small_steps)//int(self.calc.md.n_every_steps*self.calc.md.n_repeat_steps)
                #now we can check if it converted
                file = os.path.join(self.simfolder, "avg.dat")
                lx, ly, lz, ipress = np.loadtxt(file, usecols=(1, 2, 3, 4), unpack=True)
                
                lxpc = ipress
                mean = np.mean(lxpc)
                std = np.std(lxpc)
                volatom = np.mean((lx*ly*lz)/self.natoms)
                self.logger.info("At count %d mean pressure is %f with %f vol/atom"%(i+1, mean, volatom))
                
                if (np.abs(mean - self.calc._pressure)) < self.calc.tolerance.pressure:

                    #process other means
                    self.lx = np.round(np.mean(lx[-ncount+1:]), decimals=3)
                    self.ly = np.round(np.mean(ly[-ncount+1:]), decimals=3)
                    self.lz = np.round(np.mean(lz[-ncount+1:]), decimals=3)
                    self.volatom = volatom
                    self.vol = self.lx*self.ly*self.lz
                    self.logger.info("finalized vol/atom %f at pressure %f"%(self.volatom, mean))
                    self.logger.info("Avg box dimensions x: %f, y: %f, z:%f"%(self.lx, self.ly, self.lz))
                    converged = True
                    break
                laststd = std

            if not converged:
                lmp.close()
                raise ValueError("Pressure did not converge after MD runs, maybe change lattice_constant and try?")
                
            #now run for msd
            lmp.command("unfix            1")
            lmp.command("unfix            2")

        else:
            #we should do a small run to eqbrte atom positions
            lmp.command("fix              1 all nvt temp %f %f %f"%(self.calc._temperature, self.calc._temperature, self.calc.md.thermostat_damping))
            lmp.command("velocity         all create %f %d"%(self.calc._temperature, np.random.randint(0, 10000)))
            lmp.command("thermo_style     custom step pe press vol etotal temp lx ly lz")
            lmp.command("thermo           10")            
            lmp.command("fix              2 all ave/time %d %d %d v_mlx v_mly v_mlz v_mpress file avg.dat"%(int(self.calc.md.n_every_steps),
                int(self.calc.md.n_repeat_steps), int(self.calc.md.n_every_steps*self.calc.md.n_repeat_steps)))
            
            lmp.command("run              %d"%int(self.calc.md.n_small_steps))
            lmp.command("unfix            1")
            lmp.command("unfix            2")

            file = os.path.join(self.simfolder, "avg.dat")
            lx, ly, lz, ipress = np.loadtxt(file, usecols=(1, 2, 3, 4), unpack=True)
            mean = np.mean(ipress)
            self.calc._pressure = mean
            volatom = np.mean((lx*ly*lz)/self.natoms)
            self.logger.info("At count 0 mean pressure is %f with %f vol/atom"%(mean, volatom))
            self.lx = np.round(np.mean(lx), decimals=3)
            self.ly = np.round(np.mean(ly), decimals=3)
            self.lz = np.round(np.mean(lz), decimals=3)
            self.volatom = volatom
            self.vol = self.lx*self.ly*self.lz
            self.logger.info("finalized vol/atom %f at pressure %f"%(self.volatom, mean))
            self.logger.info("Avg box dimensions x: %f, y: %f, z:%f"%(self.lx, self.ly, self.lz))


        #dump is common for both
        lmp.command("dump              2 all custom 1 traj.dat id type mass x y z vx vy vz")
        lmp.command("run               0")
        lmp.command("undump            2")


        #close object and process traj
        lmp.close()
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

        #create lammps object
        lmp = ph.create_object(self.cores, self.simfolder, self.calc.md.timestep)
        
        # Adiabatic switching parameters.
        lmp.command("variable        li       equal   1.0")
        lmp.command("variable        lf       equal   0.0")
        
        #read dump file
        conf = os.path.join(self.simfolder, "conf.dump")
        lmp = ph.read_dump(lmp, conf, species=self.calc.n_elements)

        #set up hybrid potential
        #here we only need to set one potential
        lmp = ph.set_potential(lmp, self.calc)
        #lmp = ph.set_double_hybrid_potential(lmp, self.options, self.calc._pressureair_style, self.calc._pressureair_coeff)

        #remap the box to get the correct pressure
        lmp = ph.remap_box(lmp, self.lx, self.ly, self.lz)

        # Integrator & thermostat.
        if self.calc._npt:
            lmp.command("fix             f1 all npt temp %f %f %f %s %f %f %f"%(self.calc._temperature, self.calc._temperature, 
                self.calc.md.thermostat_damping, self.iso, self.calc._pressure, self.calc._pressure, self.calc.md.barostat_damping))        
        else:
            lmp.command("fix             f1 all nvt temp %f %f %f"%(self.calc._temperature, self.calc._temperature, 
                self.calc.md.thermostat_damping))

        lmp.command("thermo_style    custom step pe")
        lmp.command("thermo          1000")
        lmp.command("run             %d"%self.calc.n_equilibration_steps)
        #equilibration run is over
        
        #---------------------------------------------------------------
        # FWD cycle
        #---------------------------------------------------------------
        lmp.command("variable         flambda equal ramp(${li},${lf})")
        lmp.command("variable         blambda equal ramp(${lf},${li})")
    
        #lmp.command("pair_style       hybrid/scaled v_flambda %s v_blambda ufm 7.5"%self.options["md"]["pair_style"])
                    
        # Compute pair definitions
        if self.calc.pair_style[0] == self.calc.pair_style[1]:
            pc =  self.calc.pair_coeff[0]
            pcraw = pc.split()
            pc1 = " ".join([*pcraw[:2], *[self.calc.pair_style[0],], "1", *pcraw[2:]])
            pc =  self.calc.pair_coeff[1]
            pcraw = pc.split()
            pc2 = " ".join([*pcraw[:2], *[self.calc.pair_style[1],], "2", *pcraw[2:]])
        else:
            pc =  self.calc.pair_coeff[0]
            pcraw = pc.split()
            pc1 = " ".join([*pcraw[:2], *[self.calc.pair_style[0],], *pcraw[2:]])
            pc =  self.calc.pair_coeff[1]
            pcraw = pc.split()
            pc2 = " ".join([*pcraw[:2], *[self.calc.pair_style[1],], *pcraw[2:]])


        lmp.command("pair_style       hybrid/scaled v_flambda %s v_blambda %s"%(self.calc.pair_style[0], 
            self.calc.pair_style[1]))
        lmp.command("pair_coeff       %s"%pc1)
        lmp.command("pair_coeff       %s"%pc2)


        #apply pair force commands
        if self.calc.pair_style[0] == self.calc.pair_style[1]:
            lmp.command("compute         c1 all pair %s 1"%self.calc.pair_style[0])
            lmp.command("compute         c2 all pair %s 2"%self.calc.pair_style[1])
        else:
            lmp.command("compute         c1 all pair %s"%self.calc.pair_style[0])
            lmp.command("compute         c2 all pair %s"%self.calc.pair_style[1])


        # Output variables.
        lmp.command("variable        step equal step")
        lmp.command("variable        dU1 equal c_c1/atoms")             # Driving-force obtained from NEHI procedure.
        lmp.command("variable        dU2 equal c_c2/atoms")

        # Thermo output.
        lmp.command("thermo_style    custom step v_dU1 v_dU2")
        lmp.command("thermo          1000")


        #save the necessary items to a file: first step
        lmp.command("fix             f2 all print 1 \"${dU1} ${dU2} ${flambda}\" screen no file forward_%d.dat"%iteration)
        lmp.command("run             %d"%self.calc._n_switching_steps)


        #now equilibrate at the second potential
        lmp.command("unfix           f2")
        lmp.command("uncompute       c1")
        lmp.command("uncompute       c2")


        lmp.command("pair_style      %s"%self.calc.pair_style[1])
        lmp.command("pair_coeff      %s"%self.calc.pair_coeff[1])

        # Thermo output.
        lmp.command("thermo_style    custom step pe")
        lmp.command("thermo          1000")

        #run eqbrm run
        lmp.command("run             %d"%self.calc.n_equilibration_steps)
        
        
        #reverse switching
        lmp.command("variable         flambda equal ramp(${lf},${li})")
        lmp.command("variable         blambda equal ramp(${li},${lf})")
        
        
        lmp.command("pair_style       hybrid/scaled v_flambda %s v_blambda %s"%(self.calc.pair_style[0], 
            self.calc.pair_style[1]))
        lmp.command("pair_coeff       %s"%pc1)
        lmp.command("pair_coeff       %s"%pc2)


        #apply pair force commands
        if self.calc.pair_style[0] == self.calc.pair_style[1]:
            lmp.command("compute         c1 all pair %s 1"%self.calc.pair_style[0])
            lmp.command("compute         c2 all pair %s 2"%self.calc.pair_style[1])
        else:
            lmp.command("compute         c1 all pair %s"%self.calc.pair_style[0])
            lmp.command("compute         c2 all pair %s"%self.calc.pair_style[1])


        # Output variables.
        lmp.command("variable        step equal step")
        lmp.command("variable        dU1 equal c_c1/atoms")             # Driving-force obtained from NEHI procedure.
        lmp.command("variable        dU2 equal c_c2/atoms")

        # Thermo output.
        lmp.command("thermo_style    custom step v_dU1 v_dU2")
        lmp.command("thermo          1000")


        #save the necessary items to a file: first step
        lmp.command("fix             f2 all print 1 \"${dU1} ${dU2} ${flambda}\" screen no file backward_%d.dat"%iteration)
        lmp.command("run             %d"%self.calc._n_switching_steps)


        #now equilibrate at the second potential
        lmp.command("unfix           f2")
        lmp.command("uncompute       c1")
        lmp.command("uncompute       c2")

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
        Calculates the final work, energy dissipation; In alchemical mode, there is reference system,
        the calculated free energy is the same as the work.
        """
        w, q, qerr = find_w(self.simfolder, nelements=self.calc.n_elements, 
            concentration=self.concentration, nsims=self.calc.n_iterations, 
            full=True, solid=False, alchemy=True)
        
        self.w = w
        self.ferr = qerr
        self.fe = self.w


