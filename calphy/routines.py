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

from mendeleev import element
import copy
import numpy as np
import os

from calphy.input import read_inputfile, create_identifier
import calphy.queuekernel as cq
from calphy.errors import *
import calphy.helpers as ph


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
    def __init__(self, options=None, kernel=None, simfolder=None):
        self.options = options
        self.kernel = kernel
        self.simfolder = simfolder
        self.org_tm = 0
        self.calc = self.options['calculations'][kernel]
        self.dtemp = self.calc['dtemp']
        self.maxattempts = self.calc['maxattempts']
        self.attempts = 0
        self.exp_tm = self.calc['tguess']

        self.get_props(self.calc['element'][0])
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
        self.calc['temperature'] = int(self.tmin)
        self.calc['temperature_stop'] = int(self.tmax)
        csol = copy.deepcopy(self.calc)
        clqd = copy.deepcopy(self.calc)
        
        csol['lattice'] = self.lattice.upper()
        clqd['lattice'] = 'LQD'
        csol['state'] = 'solid'
        clqd['state'] = 'liquid'
        csol['lattice_constant'] = self.lattice_constant
        clqd['lattice_constant'] = self.lattice_constant
        csol['thigh'] = self.tmin
        clqd['thigh'] = 1.5*self.tmax
        csol['mode'] = 'ts'
        clqd['mode'] = 'ts'
        
        csol['directory'] = create_identifier(csol)
        clqd['directory'] = create_identifier(clqd)
        self.options['calculations'] = [csol, clqd]
        
        
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

        self.soljob = cq.setup_calculation(self.options, 0)
        self.lqdjob = cq.setup_calculation(self.options, 1)
        
        self.logger.info("Free energy of %s and %s phases will be calculated"%(self.soljob.calc['lattice'], self.lqdjob.calc['lattice']))
        self.logger.info("Temperature range of %f-%f"%(self.tmin, self.tmax))
        self.logger.info("STATE: Temperature range of %f-%f K"%(self.tmin, self.tmax))
        self.logger.info('Starting solid fe calculation')
        
        try:
            self.soljob = cq.routine_fe(self.soljob)
        except MeltedError:
            self.logger.info('Solid phase melted')
            return 2
        
        self.logger.info('Starting solid reversible scaling run')
        for i in range(self.soljob.nsims):
            try:
                self.soljob.reversible_scaling(iteration=(i+1))
            except MeltedError:
                self.logger.info('Solid system melted during reversible scaling run')
                return 2
            
            self.solres = self.soljob.integrate_reversible_scaling(scale_energy=True,
                                           return_values=True)
        
        self.logger.info('Starting liquid fe calculation')
        try:
            self.lqdjob = cq.routine_fe(self.lqdjob)
        except SolidifiedError:
            self.logger.info('Liquid froze')
            return 3

        self.logger.info('Starting liquid reversible scaling calculation')
        for i in range(self.lqdjob.nsims):
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
            if (self.attempts>self.maxattempts):
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
        