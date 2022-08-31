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

from mendeleev import element
import copy
import numpy as np
import os
import time

from calphy.input import read_inputfile
#import calphy.queuekernel as cq
from calphy.errors import *
import calphy.helpers as ph

from calphy.liquid import Liquid
from calphy.solid import Solid
from calphy.composition_transformation import CompositionTransformation

class MeltingTemp:
    """
    Class for automated melting temperature calculation.

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
        self.calc = calculation
        self.simfolder = simfolder
        self.org_tm = 0
        self.dtemp = self.calc.melting_temperature.step
        self.maxattempts = self.calc.melting_temperature.attempts
        self.attempts = 0
        self.exp_tm = self.calc._temperature
        self.calculations = []

        self.get_props(self.calc.element[0])
        self.get_trange()
        self.arg = None
        

        logfile = os.path.join(os.getcwd(), "calphy.log")
        self.logger = ph.prepare_log(logfile)
    
    def prepare_calcs(self):
        """
        Prepare calculations list from given object

        Parameters
        ----------
        None

        Returns
        -------
        None 
        """
        self.calc.temperature = [int(self.tmin), int(self.tmax)]
        self.calc._temperature_stop = int(self.tmax)
        csol = copy.deepcopy(self.calc)
        clqd = copy.deepcopy(self.calc)
        
        csol.lattice = self.lattice.upper()
        clqd.lattice = 'LQD'
        csol.reference_phase = 'solid'
        clqd.reference_phase = 'liquid'
        csol.lattice_constant = self.lattice_constant
        clqd.lattice_constant = self.lattice_constant
        csol._temperature_high = self.tmin
        clqd._temperature_high = 1.5*self.tmax
        csol.mode = 'ts'
        clqd.mode = 'ts'
        
        #csol['directory'] = create_identifier(csol)
        #clqd['directory'] = create_identifier(clqd)
        self.calculations = [csol, clqd]
        
        
    def get_props(self, elem):
        """
        Get properties from mendeleev

        Parameters
        ----------
        elem : string
            Chemical symbol of the element

        Returns
        -------
        None
        """
        chem = element(elem)
        lattice = chem.lattice_structure
        self.lattice_constant = chem.lattice_constant
        self.org_tm = chem.melting_point
        
        if self.exp_tm is None:
            self.exp_tm = chem.melting_point

        if lattice == "HEX":
            lattice = "HCP"
        self.lattice = lattice.lower()
        
    def get_trange(self):
        """
        Get temperature range for calculations

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        tmin = self.exp_tm - self.dtemp
        if tmin < 0:
            tmin = 10
        tmax = self.exp_tm + self.dtemp
        self.tmax = tmax
        self.tmin = tmin
        
        
    def run_jobs(self):
        """
        Run calculations

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        self.prepare_calcs()

        self.soljob = Solid(calculation=self.calculations[0], 
            simfolder=self.calculations[0].create_folders())
        self.lqdjob = Liquid(calculation=self.calculations[1], 
            simfolder=self.calculations[1].create_folders())
        
        self.logger.info("Free energy of %s and %s phases will be calculated"%(self.soljob.calc.lattice, self.lqdjob.calc.lattice))
        self.logger.info("Temperature range of %f-%f"%(self.tmin, self.tmax))
        self.logger.info("STATE: Temperature range of %f-%f K"%(self.tmin, self.tmax))
        self.logger.info('Starting solid fe calculation')
        
        try:
            self.soljob = routine_fe(self.soljob)
        except MeltedError:
            self.logger.info('Solid phase melted')
            return 2
        
        self.logger.info('Starting solid reversible scaling run')
        for i in range(self.soljob.calc.n_iterations):
            try:
                self.soljob.reversible_scaling(iteration=(i+1))
            except MeltedError:
                self.logger.info('Solid system melted during reversible scaling run')
                return 2
            
            self.solres = self.soljob.integrate_reversible_scaling(scale_energy=True,
                                           return_values=True)
        
        self.logger.info('Starting liquid fe calculation')
        try:
            self.lqdjob = routine_fe(self.lqdjob)
        except SolidifiedError:
            self.logger.info('Liquid froze')
            return 3

        self.logger.info('Starting liquid reversible scaling calculation')
        for i in range(self.lqdjob.calc.n_iterations):
            try:
                self.lqdjob.reversible_scaling(iteration=(i+1))
            except SolidifiedError:
                self.logger.info('Liquid froze during reversible scaling calculation')
                return 3

        self.lqdres = self.lqdjob.integrate_reversible_scaling(scale_energy=True,
                                           return_values=True)
    
    def start_calculation(self):
        """
        Start calculation

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        for i in range(100):
            returncode = self.run_jobs()
        
            if returncode == 3:
                self.tmin = self.tmin + self.dtemp
                self.tmax = self.tmax + self.dtemp

            elif returncode == 2:
                self.tmin = self.tmin - self.dtemp
                if self.tmin < 0:
                    self.tmin = 0
                self.tmax = self.tmax - self.dtemp
                
            else:
                return True
            
            self.attempts += 1
            if (self.attempts > self.maxattempts):
                raise ValueError('Maximum number of tries reached')


    def extrapolate_tm(self, arg):
        """
        Extrapolate Tm
        """
        solfit = np.polyfit(self.solres[0], self.solres[1], 1)
        lqdfit = np.polyfit(self.lqdres[0], self.lqdres[1], 1)

        found = False

        for i in range(100):
            if arg==0:
                self.tmin = self.tmin-self.dtemp
                if self.tmin < 0:
                    self.tmin = 0
                new_t = np.linspace(self.tmin, self.tmax, 1000)
            elif arg==999:
                self.tmax = self.tmax + self.dtemp
                new_t = np.linspace(self.tmin, self.tmax, 1000)
            else:
                break
            #evaluate in new range
            new_solfe = np.polyval(solfit, new_t)
            new_lqdfe = np.polyval(lqdfit, new_t)
    
            arg = np.argsort(np.abs(new_solfe-new_lqdfe))[0]
            tpred = new_t[arg]
        else:
            raise ValueError('failed to extrapolate melting temperature')
        
        self.logger.info("Predicted melting temperature from extrapolation: %f"%tpred)
        self.logger.info("STATE: Predicted Tm from extrapolation: %f K"%tpred)
        return tpred    
                
        
    def find_tm(self):
        """
        Find melting temperature

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        for i in range(100):
            arg = np.argsort(np.abs(self.solres[1]-self.lqdres[1]))[0]
            self.arg = arg
        
            if ((arg==0) or (arg==len(self.solres[1])-1)):
                self.logger.info('From calculation, melting temperature is not within the selected range.')
                self.logger.info('STATE: From calculation, Tm is not within range.')
                if arg==len(self.solres[1])-1:
                    arg = 999
                #the above is just a trick to extrapolate
                #now here we need to find a guess value;
                tpred = self.extrapolate_tm(arg)
                #now we have to run calcs again
                self.tmin = tpred - self.dtemp
                if self.tmin < 0:
                    self.tmin = 0
                self.tmax = tpred + self.dtemp
                self.logger.info('Restarting calculation with predicted melting temperature +/- %f'%self.dtemp)
                #self.logger.info('STATE: Restarting calculation with predicted melting temperature +/- %f'%self.dtemp)
                self.start_calculation()
                
            else:
                self.calc_tm = self.solres[0][arg]
                #get errors
                suberr = np.sqrt(self.solres[2][arg]**2 + self.lqdres[2][arg]**2)
                sol_slope = (self.solres[1][arg+50]-self.solres[1][arg-50])/(self.solres[0][arg+50]-self.solres[0][arg-50])
                lqd_slope = (self.lqdres[1][arg+50]-self.lqdres[1][arg-50])/(self.lqdres[0][arg+50]-self.lqdres[0][arg-50])
                slope_diff = sol_slope-lqd_slope
                tmerr = suberr/slope_diff
                self.tmerr = tmerr
                return self.calc_tm, self.tmerr
            
            self.attempts += 1
            self.logger.info('Attempt incremented to %d'%self.attempts)
            if self.attempts>self.maxattempts:
                raise ValueError('Maximum number of tries reached')
    
    def calculate_tm(self):
        #do a first round of calculation
        self.start_calculation()
        tm, tmerr = self.find_tm()
        self.logger.info('Found melting temperature = %.2f +/- %.2f K '%(tm, tmerr))
        self.logger.info('Experimental melting temperature = %.2f K '%(self.org_tm))
        self.logger.info('STATE: Tm = %.2f K +/- %.2f K, Exp. Tm = %.2f K'%(tm, tmerr, self.org_tm))

def routine_fe(job):
    """
    Perform an FE calculation routine
    """
    ts = time.time()
    job.run_averaging()
    te = (time.time() - ts)
    job.logger.info("Averaging routine finished in %f s"%te)

    #now run integration loops
    for i in range(job.calc.n_iterations):
        ts = time.time()
        job.run_integration(iteration=(i+1))
        te = (time.time() - ts)
        job.logger.info("Integration cycle %d finished in %f s"%(i+1, te))

    job.thermodynamic_integration()
    job.submit_report()
    return job

def routine_ts(job):
    """
    Perform ts routine
    """
    routine_fe(job)

    #now do rev scale steps
    for i in range(job.calc.n_iterations):
        ts = time.time()
        job.reversible_scaling(iteration=(i+1))
        te = (time.time() - ts)
        job.logger.info("TS integration cycle %d finished in %f s"%(i+1, te))
    
    job.integrate_reversible_scaling(scale_energy=True)
    return job


def routine_only_ts(job):
    """
    Perform sweep without free energy calculation
    """
    ts = time.time()
    job.run_averaging()
    te = (time.time() - ts)
    job.logger.info("Averaging routine finished in %f s"%te)

    for i in range(job.calc.n_iterations):
        ts = time.time()
        job.reversible_scaling(iteration=(i+1))
        te = (time.time() - ts)
        job.logger.info("TS integration cycle %d finished in %f s"%(i+1, te))
    return job

def routine_tscale(job):
    """
    Perform tscale routine
    """
    routine_fe(job)

    #now do rev scale steps
    for i in range(job.calc.n_iterations):
        ts = time.time()
        job.temperature_scaling(iteration=(i+1))
        te = (time.time() - ts)
        job.logger.info("Temperature scaling cycle %d finished in %f s"%(i+1, te))
    
    job.integrate_reversible_scaling(scale_energy=False)
    return job

def routine_pscale(job):
    """
    Perform pscale routine
    """
    routine_fe(job)

    #now do rev scale steps
    for i in range(job.calc.n_iterations):
        ts = time.time()
        job.pressure_scaling(iteration=(i+1))
        te = (time.time() - ts)
        job.logger.info("Pressure scaling cycle %d finished in %f s"%(i+1, te))
    
    job.integrate_pressure_scaling()
    return job

def routine_alchemy(job):
    """
    Perform an FE calculation routine
    """
    ts = time.time()
    job.run_averaging()
    te = (time.time() - ts)
    job.logger.info("Averaging routine finished in %f s"%te)

    #now run integration loops
    for i in range(job.calc.n_iterations):
        ts = time.time()
        job.run_integration(iteration=(i+1))
        te = (time.time() - ts)
        job.logger.info("Alchemy integration cycle %d finished in %f s"%(i+1, te))

    job.thermodynamic_integration()
    job.submit_report()
    return job 


def routine_composition_scaling(job):
    """
    Perform a compositional scaling routine
    """
    #we set up comp scaling first
    job.logger.info("Calculating composition scaling")
    comp = CompositionTransformation(job.calc.lattice, 
        job.calc.composition_scaling.input_chemical_composition, 
        job.calc.composition_scaling.output_chemical_composition, 
        restrictions=job.calc.composition_scaling.restrictions)

    #update pair styles
    res = comp.update_pair_coeff(job.calc.pair_coeff[0])
    job.calc.pair_style.append(job.calc.pair_style[0])
    job.calc.pair_coeff[0] = res[0]
    job.calc.pair_coeff.append(res[1])
    job.logger.info("Update pair coefficients")
    job.logger.info(f"pair coeff 1: {job.calc.pair_coeff[0]}")
    job.logger.info(f"pair coeff 2: {job.calc.pair_coeff[1]}")
    backup_element = job.calc.element.copy()
    job.calc.element = comp.pair_list_old
    #job.calc._ghost_element_count = len(comp.new_atomtype) - len()

    #write new file out and update lattice
    outfilename = ".".join([job.calc.lattice, "comp", "dump"])
    comp.write_structure(outfilename)
    job.calc.lattice = outfilename
    job.logger.info(f"Modified lattice written to {outfilename}")

    #prepare mass change methods
    #update and backup mass
    job.logger.info(f"Original mass: {job.calc.mass}")
    backup_mass = job.calc.mass.copy()
    mass_dict = {key:val for (key, val) in zip(backup_element, backup_mass)}

    target_masses = []
    target_counts = []

    ref_mass_list = []

    for mdict in comp.transformation_list:
        ref_mass_list.append(mass_dict[mdict["primary_element"]])
        target_masses.append(mass_dict[mdict["secondary_element"]])
        target_counts.append(mdict["count"])

    if len(backup_mass) > 2:
        job.logger.warning("Composition scaling is untested for more than 2 elements!")

    if len(np.unique(ref_mass_list)) > 1:
        job.logger.warning("More than one kind of transformation found! Stopping")
        raise RuntimeError("More than one kind of transformation found! Stopping")

    ref_mass = ref_mass_list[0] 
    
    #now replace mass    
    job.calc.mass = [ref_mass for x in range(len(job.calc.element))] 
    job.logger.info(f"Temporarily replacing mass: {job.calc.mass}")


    #now start cycle
    ts = time.time()
    job.run_averaging()
    te = (time.time() - ts)
    job.logger.info("Averaging routine finished in %f s"%te)

    #now run integration loops
    for i in range(job.calc.n_iterations):
        ts = time.time()
        job.run_integration(iteration=(i+1))
        te = (time.time() - ts)
        job.logger.info("Alchemy integration cycle %d finished in %f s"%(i+1, te))

    flambda_arr, w_arr, q_arr, qerr_arr = job.thermodynamic_integration()

    #read the file
    mcorrarr, mcorsum = job.mass_integration(flambda_arr, ref_mass, target_masses, target_counts)
    netfe = w_arr - mcorrarr

    job.fe = job.fe - mcorsum
    job.submit_report(extra_dict = {"results":{"mass_correction": float(mcorsum)}})

    outfile = os.path.join(job.simfolder, "composition_sweep.dat")
    np.savetxt(outfile, np.column_stack((flambda_arr, netfe, w_arr, mcorrarr)))
    return job