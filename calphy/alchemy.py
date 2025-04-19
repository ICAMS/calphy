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

Notes
-----
- swapping is strictly only performed between types 1 and 2 at the moment; this needs to be refined further
"""

import numpy as np
import yaml

from calphy.integrators import *
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
    def __init__(self, calculation=None, simfolder=None, log_to_screen=False):

        #call base class
        super().__init__(calculation=calculation,
        simfolder=simfolder,
        log_to_screen=log_to_screen)


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
        lmp = ph.create_object(self.cores, self.simfolder, self.calc.md.timestep, 
            self.calc.md.cmdargs, self.calc.md.init_commands)

        lmp.command(f'pair_style {self.calc._pair_style_with_options[0]}')

        #set up structure
        lmp = ph.create_structure(lmp, self.calc)

        #set up potential
        lmp.command(f'pair_coeff {self.calc.pair_coeff[0]}')
        lmp = ph.set_mass(lmp, self.calc)

        #add some computes
        lmp.command("variable         mvol equal vol")
        lmp.command("variable         mlx equal lx")
        lmp.command("variable         mly equal ly")
        lmp.command("variable         mlz equal lz")
        lmp.command("variable         mpress equal press")

        #add some computes
        if not self.calc._fix_lattice:
            if self.calc._pressure == 0:
                self.run_zero_pressure_equilibration(lmp)
            else:
                self.run_finite_pressure_equilibration(lmp)


            #this is when the averaging routine starts
            self.run_pressure_convergence(lmp)
        
        #run if a constrained lattice is used
        else:
            #routine in which lattice constant will not varied, but is set to a given fixed value
            self.run_constrained_pressure_convergence(lmp)    

        #check for melting
        self.dump_current_snapshot(lmp, "traj.equilibration_stage2.dat")
        self.check_if_melted(lmp, "traj.equilibration_stage2.dat")

        #close object and process traj
        lmp = ph.write_data(lmp, "conf.equilibration.data")
        lmp.close()


    

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
        lmp = ph.create_object(self.cores, self.simfolder, self.calc.md.timestep, 
            self.calc.md.cmdargs, self.calc.md.init_commands)
        
        # Adiabatic switching parameters.
        lmp.command("variable        li       equal   1.0")
        lmp.command("variable        lf       equal   0.0")

        lmp.command(f'pair_style {self.calc._pair_style_with_options[0]}')
        
        #read dump file
        #conf = os.path.join(self.simfolder, "conf.equilibration.dump")
        conf = os.path.join(self.simfolder, "conf.equilibration.data")
        lmp = ph.read_data(lmp, conf)

        #set up hybrid potential
        #here we only need to set one potential
        lmp.command(f'pair_coeff {self.calc.pair_coeff[0]}')
        lmp = ph.set_mass(lmp, self.calc)

        #NEW ADDED
        lmp.command("group g1 type 1")
        lmp.command("group g2 type 2")
        #lmp = ph.set_double_hybrid_potential(lmp, self.options, self.calc._pressureair_style, self.calc._pressureair_coeff)

        #remap the box to get the correct pressure
        lmp = ph.remap_box(lmp, self.lx, self.ly, self.lz)

        lmp.command("velocity          all create %f %d mom yes rot yes dist gaussian"%(self.calc._temperature, np.random.randint(1, 10000)))
        # Integrator & thermostat.
        if self.calc.npt:
            lmp.command("fix             f1 all npt temp %f %f %f %s %f %f %f"%(self.calc._temperature, self.calc._temperature, 
                self.calc.md.thermostat_damping[1], self.iso, self.calc._pressure, self.calc._pressure, self.calc.md.barostat_damping[1]))        
        else:
            lmp.command("fix             f1 all nvt temp %f %f %f"%(self.calc._temperature, self.calc._temperature, 
                self.calc.md.thermostat_damping[1]))

        
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
            pc1 = " ".join([*pcraw[:2], *[self.calc._pair_style_names[0],], "1", *pcraw[2:]])
            pc =  self.calc.pair_coeff[1]
            pcraw = pc.split()
            pc2 = " ".join([*pcraw[:2], *[self.calc._pair_style_names[1],], "2", *pcraw[2:]])
        else:
            pc =  self.calc.pair_coeff[0]
            pcraw = pc.split()
            pc1 = " ".join([*pcraw[:2], *[self.calc._pair_style_names[0],], *pcraw[2:]])
            pc =  self.calc.pair_coeff[1]
            pcraw = pc.split()
            pc2 = " ".join([*pcraw[:2], *[self.calc._pair_style_names[1],], *pcraw[2:]])


        lmp.command("pair_style       hybrid/scaled v_flambda %s v_blambda %s"%(
            self.calc._pair_style_with_options[0],
            self.calc._pair_style_with_options[1]
            )
        )
        lmp.command("pair_coeff       %s"%pc1)
        lmp.command("pair_coeff       %s"%pc2)


        #apply pair force commands
        if self.calc._pair_style_names[0] == self.calc._pair_style_names[1]:
            lmp.command("compute         c1 all pair %s 1"%self.calc._pair_style_names[0])
            lmp.command("compute         c2 all pair %s 2"%self.calc._pair_style_names[1])
        else:
            lmp.command("compute         c1 all pair %s"%self.calc._pair_style_names[0])
            lmp.command("compute         c2 all pair %s"%self.calc._pair_style_names[1])


        # Output variables.
        lmp.command("variable        step equal step")
        lmp.command("variable        dU1 equal c_c1/atoms")             # Driving-force obtained from NEHI procedure.
        lmp.command("variable        dU2 equal c_c2/atoms")

        #add swaps if n_swap is > 0
        if self.calc.monte_carlo.n_swaps > 0:
            self.logger.info(f'{self.calc.monte_carlo.n_swaps} swap moves are performed between {self.calc.monte_carlo.swap_types[0]} and {self.calc.monte_carlo.swap_types[1]} every {self.calc.monte_carlo.n_steps}')
            lmp.command("fix  swap all atom/swap %d %d %d %f ke no types %d %d"%(self.calc.monte_carlo.n_steps,
                                                                                self.calc.monte_carlo.n_swaps,
                                                                                np.random.randint(1, 10000),
                                                                                self.calc._temperature,
                                                                                self.calc.monte_carlo.swap_types[0],
                                                                                self.calc.monte_carlo.swap_types[1]))
            lmp.command("variable a equal f_swap[1]")
            lmp.command("variable b equal f_swap[2]")
            lmp.command("fix             swap2 all print 1 \"${a} ${b} ${flambda}\" screen no file swap.forward_%d.dat"%iteration)

        # Thermo output.
        if self.calc.monte_carlo.n_swaps > 0:
            lmp.command("thermo_style    custom step v_dU1 v_dU2 v_a v_b")
        else:
            lmp.command("thermo_style    custom step v_dU1 v_dU2")
        lmp.command("thermo          1000")

        
        #save the necessary items to a file: first step
        lmp.command("fix             f2 all print 1 \"${dU1} ${dU2} ${flambda}\" screen no file forward_%d.dat"%iteration)
        lmp.command("run             %d"%self.calc._n_switching_steps)            

        #now equilibrate at the second potential
        lmp.command("unfix           f2")
        lmp.command("uncompute       c1")
        lmp.command("uncompute       c2")

        #NEW SWAP
        if self.calc.monte_carlo.n_swaps > 0:
            lmp.command("unfix swap")
            lmp.command("unfix swap2")


        lmp.command("pair_style      %s"%self.calc._pair_style_with_options[1])
        lmp.command("pair_coeff      %s"%self.calc.pair_coeff[1])

        # Thermo output.
        lmp.command("thermo_style    custom step pe")
        lmp.command("thermo          1000")
        
        #run eqbrm run
        lmp.command("run             %d"%self.calc.n_equilibration_steps)
        
        
        #reverse switching
        lmp.command("variable         flambda equal ramp(${lf},${li})")
        lmp.command("variable         blambda equal ramp(${li},${lf})")
        
        
        lmp.command("pair_style       hybrid/scaled v_flambda %s v_blambda %s"%(self.calc._pair_style_with_options[0], 
            self.calc._pair_style_with_options[1]))
        lmp.command("pair_coeff       %s"%pc1)
        lmp.command("pair_coeff       %s"%pc2)


        #apply pair force commands
        if self.calc._pair_style_names[0] == self.calc._pair_style_names[1]:
            lmp.command("compute         c1 all pair %s 1"%self.calc._pair_style_names[0])
            lmp.command("compute         c2 all pair %s 2"%self.calc._pair_style_names[1])
        else:
            lmp.command("compute         c1 all pair %s"%self.calc._pair_style_names[0])
            lmp.command("compute         c2 all pair %s"%self.calc._pair_style_names[1])


        # Output variables.
        lmp.command("variable        step equal step")
        lmp.command("variable        dU1 equal c_c1/atoms")             # Driving-force obtained from NEHI procedure.
        lmp.command("variable        dU2 equal c_c2/atoms")

        #add swaps if n_swap is > 0
        if self.calc.monte_carlo.n_swaps > 0:
            if self.calc.monte_carlo.reverse_swap:
                self.logger.info(f'{self.calc.monte_carlo.n_swaps} swap moves are performed between {self.calc.monte_carlo.swap_types[1]} and {self.calc.monte_carlo.swap_types[0]} every {self.calc.monte_carlo.n_steps}')
                lmp.command("fix  swap all atom/swap %d %d %d %f ke no types %d %d"%(self.calc.monte_carlo.n_steps,
                                                                                    self.calc.monte_carlo.n_swaps,
                                                                                    np.random.randint(1, 10000),
                                                                                    self.calc._temperature,
                                                                                    self.calc.monte_carlo.swap_types[1],
                                                                                    self.calc.monte_carlo.swap_types[0]))
            else:
                self.logger.info(f'{self.calc.monte_carlo.n_swaps} swap moves are performed between {self.calc.monte_carlo.swap_types[0]} and {self.calc.monte_carlo.swap_types[1]} every {self.calc.monte_carlo.n_steps}')
                self.logger.info('note that swaps are not reversed')
                lmp.command("fix  swap all atom/swap %d %d %d %f ke no types %d %d"%(self.calc.monte_carlo.n_steps,
                                                                                    self.calc.monte_carlo.n_swaps,
                                                                                    np.random.randint(1, 10000),
                                                                                    self.calc._temperature,
                                                                                    self.calc.monte_carlo.swap_types[0],
                                                                                    self.calc.monte_carlo.swap_types[1]))

            lmp.command("variable a equal f_swap[1]")
            lmp.command("variable b equal f_swap[2]")
            lmp.command("fix             swap2 all print 1 \"${a} ${b} ${blambda}\" screen no file swap.backward_%d.dat"%iteration)

        # Thermo output.
        if self.calc.monte_carlo.n_swaps > 0:
            lmp.command("thermo_style    custom step v_dU1 v_dU2 v_a v_b")
        else:        
            lmp.command("thermo_style    custom step v_dU1 v_dU2")
        lmp.command("thermo          1000")


        
        #save the necessary items to a file: first step
        lmp.command("fix             f2 all print 1 \"${dU1} ${dU2} ${flambda}\" screen no file backward_%d.dat"%iteration)
        lmp.command("run             %d"%self.calc._n_switching_steps)


        #now equilibrate at the second potential
        lmp.command("unfix           f2")
        lmp.command("uncompute       c1")
        lmp.command("uncompute       c2")

        if self.calc.monte_carlo.n_swaps > 0:
            lmp.command("unfix  swap")
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
        w, q, qerr = find_w(self.simfolder, self.calc,
            full=True, solid=False)

        self.w = w
        self.ferr = qerr
        self.fe = self.w
        
        if self.calc.mode == "composition_scaling":
            w_arr, q_arr, qerr_arr, flambda_arr = find_w(self.simfolder, self.calc,
                full=True, solid=False, composition_integration=True)

            #now we need to process the comp scaling
            return flambda_arr, w_arr, q_arr, qerr_arr


    def mass_integration(self, flambda, ref_mass, target_masses, target_counts):
        mcorarr, mcorsum = integrate_mass(flambda, ref_mass, target_masses, target_counts,
    self.calc._temperature, self.natoms)
        return mcorarr, mcorsum      


