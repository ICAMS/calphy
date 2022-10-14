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
import copy
import sys

import pyscal.traj_process as ptp
from calphy.integrators import *
import calphy.lattice as pl
import calphy.helpers as ph
import calphy.phase as cph
from calphy.errors import *

class Solid(cph.Phase):
    """
    Class for free energy calculation with solid as the reference state

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

    def run_spring_constant_convergence(self, lmp):
        """
        """
        lmp.command("fix              3 all nvt temp %f %f %f"%(self.calc._temperature, self.calc._temperature, self.calc.md.thermostat_damping[1]))
        
        #apply fix
        lmp = ph.compute_msd(lmp, self.calc)
        
        if ph.check_if_any_is_none(self.calc.spring_constants):
            #similar averaging routine
            laststd = 0.00
            for i in range(self.calc.md.n_cycles):
                lmp.command("run              %d"%int(self.calc.md.n_small_steps))
                ncount = int(self.calc.md.n_small_steps)//int(self.calc.md.n_every_steps*self.calc.md.n_repeat_steps)
                #now we can check if it converted
                file = os.path.join(self.simfolder, "msd.dat")
                quant = np.loadtxt(file, usecols=(1,), unpack=True)[-ncount+1:]
                quant = 3*kb*self.calc._temperature/quant
                #self.logger.info(quant)
                mean = np.mean(quant)
                std = np.std(quant)
                self.logger.info("At count %d mean k is %f std is %f"%(i+1, mean, std))
                if (np.abs(laststd - std) < self.calc.tolerance.spring_constant):
                    #now reevaluate spring constants
                    k = []
                    for i in range(self.calc.n_elements):
                        quant = np.loadtxt(file, usecols=(i+1, ), unpack=True)[-ncount+1:]
                        quant = 3*kb*self.calc._temperature/quant
                        k.append(np.round(np.mean(quant), decimals=2))

                    #first replace any provided values with user values
                    if ph.check_if_any_is_not_none(self.calc.spring_constants):
                        spring_constants = copy.copy(self.calc.spring_constants)
                        k = ph.replace_nones(spring_constants, k, logger=self.logger)
                    
                    #add sanity checks
                    k = ph.validate_spring_constants(k, logger=self.logger)
                    
                    #now save
                    self.k = k

                    self.logger.info("finalized sprint constants")
                    self.logger.info(self.k)
                    break
                laststd = std

        else:
            if not (len(self.calc.spring_constants) == self.calc.n_elements):
                raise ValueError("Spring constant input length should be same as number of elements, spring constant length %d, # elements %d"%(len(self.calc.spring_constants), self.calc.n_elements))

            #still run a small NVT cycle
            lmp.command("run              %d"%int(self.calc.md.n_small_steps))
            self.k = self.calc.spring_constants
            self.logger.info("Used user input sprint constants")
            self.logger.info(self.k)

        lmp.command("unfix         3")

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
        lmp = ph.create_object(self.cores, self.simfolder, self.calc.md.timestep, self.calc.md.cmdargs)

        #set up structure
        lmp = ph.create_structure(lmp, self.calc, species=self.calc.n_elements+self.calc._ghost_element_count)

        #set up potential
        if self.calc.potential_file is None:
            lmp = ph.set_potential(lmp, self.calc)
        else:
            lmp.command("include %s"%self.calc.potential_file)

        #add some computes
        lmp.command("variable         mvol equal vol")
        lmp.command("variable         mlx equal lx")
        lmp.command("variable         mly equal ly")
        lmp.command("variable         mlz equal lz")
        lmp.command("variable         mpress equal press")

        #Run if a constrained lattice is not needed
        if not self.calc._fix_lattice:
            if self.calc._pressure == 0:
                self.run_zero_pressure_equilibration(lmp)
            else:
                self.run_finite_pressure_equilibration(lmp)


            #this is when the averaging routine starts
            self.run_pressure_convergence(lmp)

            #dump snapshot and check if melted
            self.dump_current_snapshot(lmp, "traj.equilibration_stage1.dat")
            self.check_if_melted(lmp, "traj.equilibration_stage1.dat")
        
        #run if a constrained lattice is used
        else:
            #routine in which lattice constant will not varied, but is set to a given fixed value
            self.run_constrained_pressure_convergence(lmp)

        #start MSD calculation routine
        #there two possibilities here - if spring constants are provided, use it. If not, calculate it
        self.run_spring_constant_convergence(lmp)

        #check for melting
        self.dump_current_snapshot(lmp, "traj.equilibration_stage2.dat")
        self.check_if_melted(lmp, "traj.equilibration_stage2.dat")

        #close object and process traj
        lmp.close()
        self.process_traj("traj.equilibration_stage2.dat", "conf.equilibration.dump")



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
        lmp = ph.create_object(self.cores, self.simfolder, self.calc.md.timestep, self.calc.md.cmdargs)

        #read in the conf file
        conf = os.path.join(self.simfolder, "conf.equilibration.dump")
        lmp = ph.read_dump(lmp, conf, species=self.calc.n_elements+self.calc._ghost_element_count)

        #set up potential
        if self.calc.potential_file is None:
            lmp = ph.set_potential(lmp, self.calc)
        else:
            lmp.command("include %s"%self.calc.potential_file)

        #remap the box to get the correct pressure
        lmp = ph.remap_box(lmp, self.lx, self.ly, self.lz)

        #create groups - each species belong to one group
        for i in range(self.calc.n_elements):
            lmp.command("group  g%d type %d"%(i+1, i+1))

        #get counts of each group
        for i in range(self.calc.n_elements):
            lmp.command("variable   count%d equal count(g%d)"%(i+1, i+1))

        #initialise everything
        lmp.command("run               0")

        #apply initial fixes
        lmp.command("fix               f1 all nve")
        
        #apply fix for each spring
        #TODO: Add option to select function
        for i in range(self.calc.n_elements):
            lmp.command("fix               ff%d g%d ti/spring 10.0 100 100 function 2"%(i+1, i+1))
        
        #apply temp fix
        lmp.command("fix               f3 all langevin %f %f %f %d zero yes"%(self.calc._temperature, self.calc._temperature, self.calc.md.thermostat_damping[1], 
                                        np.random.randint(0, 10000)))

        #compute com and apply to fix
        lmp.command("compute           Tcm all temp/com")
        lmp.command("fix_modify        f3 temp Tcm")

        lmp.command("variable          step    equal step")
        lmp.command("variable          dU1      equal pe/atoms")
        for i in range(self.calc.n_elements):
            lmp.command("variable          dU%d      equal f_ff%d/v_count%d"%(i+2, i+1, i+1))
        
        lmp.command("variable          lambda  equal f_ff1[1]")

        #add thermo command to force variable evaluation
        lmp.command("thermo_style      custom step pe c_Tcm")
        lmp.command("thermo            10000")

        #Create velocity
        lmp.command("velocity          all create %f %d mom yes rot yes dist gaussian"%(self.calc._temperature, np.random.randint(0, 10000)))

        #reapply 
        for i in range(self.calc.n_elements):
            lmp.command("fix               ff%d g%d ti/spring %f %d %d function 2"%(i+1, i+1, self.k[i], 
                self.calc._n_switching_steps, self.calc.n_equilibration_steps))

        #Equilibriate structure
        lmp.command("run               %d"%self.calc.n_equilibration_steps)
        
        #write out energy
        str1 = "fix f4 all print 1 \"${dU1} "
        str2 = []
        for i in range(self.calc.n_elements):
            str2.append("${dU%d}"%(i+2))
        str2.append("${lambda}\"")
        str2 = " ".join(str2)
        str3 = " screen no file forward_%d.dat"%iteration
        command = str1 + str2 + str3
        lmp.command(command)

        #Forward switching over ts steps
        lmp.command("run               %d"%self.calc._n_switching_steps)
        lmp.command("unfix             f4")

        #Equilibriate
        lmp.command("run               %d"%self.calc.n_equilibration_steps)

        #write out energy
        str1 = "fix f4 all print 1 \"${dU1} "
        str2 = []
        for i in range(self.calc.n_elements):
            str2.append("${dU%d}"%(i+2))
        str2.append("${lambda}\"")
        str2 = " ".join(str2)
        str3 = " screen no file backward_%d.dat"%iteration
        command = str1 + str2 + str3
        lmp.command(command)

        #Reverse switching over ts steps
        lmp.command("run               %d"%self.calc._n_switching_steps)
        lmp.command("unfix             f4")

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
        matching with Einstein crystal
        """
        f1 = get_einstein_crystal_fe(self.calc._temperature, 
            self.natoms, self.calc.mass, 
            self.vol, self.k, self.concentration)
        w, q, qerr = find_w(self.simfolder, 
            nelements=self.calc.n_elements, 
            concentration=self.concentration, nsims=self.calc.n_iterations, 
            full=True, solid=True)
        
        self.fref = f1
        self.w = w
        self.ferr = qerr

        #add pressure contribution if required
        if self.calc._pressure != 0:
            p = self.calc._pressure/(10000*160.21766208)
            v = self.vol/self.natoms
            self.pv = p*v
        else:
            self.pv = 0 

        #calculate final free energy
        self.fe = self.fref + self.w + self.pv

        