"""
calphy: a Python library and command line interface for automated free
energy calculations.

Copyright 2021  (c) Sarath Menon^1, Yury Lysogorskiy^2, Ralf Drautz^2
^1: Max Planck Institut für Eisenforschung, Dusseldorf, Germany 
^2: Ruhr-University Bochum, Bochum, Germany

calphy is published and distributed under the Academic Software License v1.0 (ASL). 
calphy is distributed in the hope that it will be useful for non-commercial academic research, 
but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the ASL for more details. 

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

import pyscal.traj_process as ptp
from calphy.integrators import *
import calphy.lattice as pl
import calphy.helpers as ph
from calphy.errors import *

class Phase:
    """
    Class for free energy calculation.

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

        self.options = copy.deepcopy(options)
        self.simfolder = simfolder
        self.kernel = kernel
        
        logfile = os.path.join(self.simfolder, "calphy.log")
        self.logger = ph.prepare_log(logfile)

        self.calc = copy.deepcopy(options["calculations"][kernel])
        self.nsims = self.calc["nsims"]

        self.t = self.calc["temperature"]
        self.tend = self.calc["temperature_stop"]
        self.thigh = self.calc["thigh"] 
        self.p = self.calc["pressure"]
        if "pressure_stop" in self.calc.keys():
            self.pend = self.calc["pressure_stop"]
        else:
            self.pend = self.p
        self.logger.info("Temperature start: %f K, temperature stop: %f K, pressure: %f bar"%(self.t, self.tend, self.p))

        if self.calc["iso"]:
            self.iso = "iso"
        else:
            self.iso = "aniso"
        self.logger.info("Pressure adjusted in %s"%self.iso)

        self.l = None
        self.alat = None
        self.apc = None
        self.vol = None
        self.concentration = None
        self.prepare_lattice()

        #other properties
        self.cores = self.options["queue"]["cores"]
        self.ncells = np.prod(self.calc["repeat"])
        self.natoms = self.ncells*self.apc        
        self.logger.info("%d atoms in %d cells on %d cores"%(self.natoms, self.ncells, self.cores))

        #reference system props; may not be always used
        #TODO : Add option to customize UFM parameters
        self.eps = self.t*50.0*kb

        #properties that will be calculated later
        self.volatom = None
        self.k = None
        self.rho = None

        self.ferr = 0
        self.fref = 0
        self.fideal = 0
        
        self.w = 0
        self.pv = 0
        self.fe = 0

        #box dimensions that need to be stored
        self.lx = None
        self.ly = None
        self.lz = None

        #backup pair styles
        self.pair_style = self.options["md"]["pair_style"]
        self.pair_coeff = self.options["md"]["pair_coeff"]

        #now manually tune pair styles
        self.options["md"]["pair_style"] = self.options["md"]["pair_style"][0]
        self.options["md"]["pair_coeff"] = self.options["md"]["pair_coeff"][0]
        self.logger.info("pair_style: %s"%self.options["md"]["pair_style"])
        self.logger.info("pair_coeff: %s"%self.options["md"]["pair_coeff"])

        #log second pair style
        if len(self.pair_style)>1:
            self.logger.info("second pair_style: %s"%self.pair_style[1])
            self.logger.info("second pair_coeff: %s"%self.pair_coeff[1])

    def __repr__(self):
        """
        String of the class
        """
        data = "%s system with T=%f, P=%f in %s lattice for mode %s"%(self.calc["state"],
            self.t, self.p, self.l, self.calc["mode"]) 
        return data

    def prepare_lattice(self):
        """
        Prepare the lattice for the simulation

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        Calculates the lattic, lattice constant, number of atoms per unit cell
        and concentration of the input system.
        """
        l, alat, apc, conc = pl.prepare_lattice(self.calc)
        self.l = l
        self.alat = alat
        self.apc = apc
        self.concentration = conc
        self.logger.info("Lattice: %s with a=%f"%(self.l, self.alat))
        self.logger.info("%d atoms in the unit cell"%self.apc)
        self.logger.info("concentration:")
        self.logger.info(self.concentration)


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
        report["average"]["vol_atom"] = float(self.volatom)
        
        if self.k is not None:
            report["average"]["spring_constant"] = " ".join(np.array(self.k).astype(str))
        if self.rho is not None:
            report["average"]["density"] = float(self.rho)

        #results
        report["results"] = {}
        report["results"]["free_energy"] = float(self.fe)
        report["results"]["error"] = float(self.ferr)
        report["results"]["reference_system"] = float(self.fref)
        report["results"]["work"] = float(self.w)
        report["results"]["pv"] = float(self.pv)
        self.report = report

        reportfile = os.path.join(self.simfolder, "report.yaml")
        with open(reportfile, 'w') as f:
            yaml.dump(report, f)

        self.logger.info("Report written in %s"%reportfile)

        #now we have to write out the results
        self.logger.info("Please cite the following publications:")
        self.logger.info("- 10.1103/PhysRevMaterials.5.103801")

        if self.calc["mode"] == "fe":
            if self.calc["state"] == "solid":
                self.logger.info("- 10.1016/j.commatsci.2015.10.050")
            else:
                self.logger.info("- 10.1016/j.commatsci.2018.12.029")
                self.logger.info("- 10.1063/1.4967775")
        elif self.calc["mode"] == "ts":
            self.logger.info("- 10.1103/PhysRevLett.83.3973")
        elif self.calc["mode"] == "mts":
            self.logger.info("- 10.1063/1.1420486") 


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
        solid = False
        if self.calc['state'] == 'solid':
            solid = True

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
        lmp.command("fix               f1 all nph %s %f %f %f"%(self.iso, self.p, self.p, self.options["md"]["pdamp"]))
        lmp.command("fix               f2 all langevin %f %f %f %d zero yes"%(self.t, self.t, self.options["md"]["tdamp"], np.random.randint(0, 10000)))
        lmp.command("run               %d"%self.options["md"]["te"])
        lmp.command("unfix             f1")
        lmp.command("unfix             f2")

        #now fix com
        lmp.command("variable         xcm equal xcm(all,x)")
        lmp.command("variable         ycm equal xcm(all,y)")
        lmp.command("variable         zcm equal xcm(all,z)")
        
        lmp.command("fix              f1 all nph %s %f %f %f fixedpoint ${xcm} ${ycm} ${zcm}"%(self.iso, self.p, self.p, self.options["md"]["pdamp"]))
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


        lmp.command("variable         flambda equal ramp(${li},${lf})")
        lmp.command("variable         blambda equal ramp(${lf},${li})")
        lmp.command("variable         fscale equal v_flambda-1.0")
        lmp.command("variable         bscale equal v_blambda-1.0")
        lmp.command("variable         one equal 1.0")

        #set up potential
        pc =  self.options["md"]["pair_coeff"]
        pcraw = pc.split()
        pcnew1 = " ".join([*pcraw[:2], *[self.options["md"]["pair_style"],], "1", *pcraw[2:]])
        pcnew2 = " ".join([*pcraw[:2], *[self.options["md"]["pair_style"],], "2", *pcraw[2:]])

        lmp.command("pair_style       hybrid/scaled v_one %s v_fscale %s"%(self.options["md"]["pair_style"], self.options["md"]["pair_style"]))
        lmp.command("pair_coeff       %s"%pcnew1)
        lmp.command("pair_coeff       %s"%pcnew2)

        #start scaling over switching time
        lmp.command("fix              f1 all nph %s %f %f %f fixedpoint ${xcm} ${ycm} ${zcm}"%(self.iso, pi, 
            pf, self.options["md"]["pdamp"]))
        lmp.command("fix_modify        f1 temp tcm")
        lmp.command("fix               f3 all print 1 \"${dU} $(press) $(vol) ${flambda}\" screen no file forward_%d.dat"%iteration)

        if self.options["md"]["traj_interval"] > 0:
            lmp.command("dump              d1 all custom %d traj.forward_%d.dat id type mass x y z vx vy vz"%(self.options["md"]["traj_interval"],
                iteration))

        lmp.command("run               %d"%self.options["md"]["ts"])

        #unfix
        lmp.command("unfix             f3")
        lmp.command("unfix             f1")

        if self.options["md"]["traj_interval"] > 0:
            lmp.command("undump           d1")

        #switch potential
        lmp = ph.set_potential(lmp, self.options)

        #equilibriate scaled hamiltonian
        lmp.command("fix              f1 all nph %s %f %f %f fixedpoint ${xcm} ${ycm} ${zcm}"%(self.iso, pf, 
            pf, self.options["md"]["pdamp"]))
        lmp.command("fix_modify        f1 temp tcm")        
        lmp.command("run               %d"%self.options["md"]["te"])
        lmp.command("unfix             f1")


        #check melting or freezing
        lmp.command("dump              2 all custom 1 traj.dat id type mass x y z vx vy vz")
        lmp.command("run               0")
        lmp.command("undump            2")
        
        solids = ph.find_solid_fraction(os.path.join(self.simfolder, "traj.dat"))
        if solid:
            if (solids/lmp.natoms < self.options["conv"]["solid_frac"]):
                lmp.close()
                raise MeltedError("System melted, increase size or reduce scaling!")
        else:
            if (solids/lmp.natoms > self.options["conv"]["liquid_frac"]):
                lmp.close()
                raise SolidifiedError('System solidified, increase temperature')


        #reverse scaling
        lmp.command("variable         flambda equal ramp(${li},${lf})")
        lmp.command("variable         blambda equal ramp(${lf},${li})")
        lmp.command("variable         fscale equal v_flambda-1.0")
        lmp.command("variable         bscale equal v_blambda-1.0")
        lmp.command("variable         one equal 1.0")

        lmp.command("pair_style       hybrid/scaled v_one %s v_bscale %s"%(self.options["md"]["pair_style"], self.options["md"]["pair_style"]))
        lmp.command("pair_coeff       %s"%pcnew1)
        lmp.command("pair_coeff       %s"%pcnew2)

        #apply fix and perform switching        
        lmp.command("fix              f1 all nph %s %f %f %f fixedpoint ${xcm} ${ycm} ${zcm}"%(self.iso, pf, 
            pi, self.options["md"]["pdamp"]))
        lmp.command("fix_modify        f1 temp tcm")
        lmp.command("fix               f3 all print 1 \"${dU} $(press) $(vol) ${blambda}\" screen no file backward_%d.dat"%iteration)

        if self.options["md"]["traj_interval"] > 0:
            lmp.command("dump              d1 all custom %d traj.backward_%d.dat id type mass x y z vx vy vz"%(self.options["md"]["traj_interval"],
                iteration))

        lmp.command("run               %d"%self.options["md"]["ts"])
        
        lmp.command("unfix             f3")
        lmp.command("unfix             f1")

        if self.options["md"]["traj_interval"] > 0:
            lmp.command("undump           d1")
        
        #close the object
        lmp.close()

    def integrate_reversible_scaling(self, scale_energy=False, return_values=False):
        """
        Perform integration after reversible scaling

        Parameters
        ----------
        scale_energy : bool, optional
            If True, scale the energy during reversible scaling. 

        return_values : bool, optional
            If True, return integrated values

        Returns
        -------
        res : list of lists of shape 1x3
            Only returned if `return_values` is True.
        """
        res = integrate_rs(self.simfolder, self.fe, self.t, self.natoms, p=self.p,
            nsims=self.nsims, scale_energy=scale_energy, return_values=return_values)

        if return_values:
            return res

    def temperature_scaling(self, iteration=1):
        """
        Perform temperature scaling calculation in NPT
        
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
        p0 = self.p
        pf = lf*p0

        #create lammps object
        lmp = ph.create_object(self.cores, self.simfolder, self.options["md"]["timestep"])

        lmp.command("echo              log")
        lmp.command("variable          li equal %f"%li)
        lmp.command("variable          lf equal %f"%lf)

        #read in conf
        conf = os.path.join(self.simfolder, "conf.dump")
        lmp = ph.read_dump(lmp, conf, species=self.options["nelements"])

        #set up potential
        lmp = ph.set_potential(lmp, self.options)

        #remap the box to get the correct pressure
        lmp = ph.remap_box(lmp, self.lx, self.ly, self.lz)


        #equilibrate first
        lmp.command("fix               1 all npt temp %f %f %f %s %f %f %f"%(t0, t0, self.options["md"]["tdamp"],
                                        self.iso, p0, p0, self.options["md"]["pdamp"]))
        lmp.command("run               %d"%self.options["md"]["te"])
        lmp.command("unfix             1")


        #now scale system to final temp, thereby recording enerfy at every step
        lmp.command("variable          step    equal step")
        lmp.command("variable          dU      equal pe/atoms")
        lmp.command("variable          lambda equal ramp(${li},${lf})")

        lmp.command("fix               f2 all npt temp %f %f %f %s %f %f %f"%(t0, tf, self.options["md"]["tdamp"],
                                        self.iso, p0, pf, self.options["md"]["pdamp"]))
        lmp.command("fix               f3 all print 1 \"${dU} $(press) $(vol) ${lambda}\" screen no file forward_%d.dat"%iteration)
        lmp.command("run               %d"%self.options["md"]["ts"])

        lmp.command("unfix             f2")
        lmp.command("unfix             f3")

        lmp.command("fix               1 all npt temp %f %f %f %s %f %f %f"%(tf, tf, self.options["md"]["tdamp"],
                                        self.iso, pf, pf, self.options["md"]["pdamp"]))
        lmp.command("run               %d"%self.options["md"]["te"])
        lmp.command("unfix             1")

        #start reverse loop
        lmp.command("variable          lambda equal ramp(${lf},${li})")

        lmp.command("fix               f2 all npt temp %f %f %f %s %f %f %f"%(tf, t0, self.options["md"]["tdamp"],
                                        self.iso, pf, p0, self.options["md"]["pdamp"]))
        lmp.command("fix               f3 all print 1 \"${dU} $(press) $(vol) ${lambda}\" screen no file backward_%d.dat"%iteration)
        lmp.command("run               %d"%self.options["md"]["ts"])


    def pressure_scaling(self, iteration=1):
        """
        Perform pressure scaling calculation in NPT
        
        Parameters
        ----------
        iteration : int, optional
            iteration of the calculation. Default 1
        
        Returns
        -------
        None
        """
        t0 = self.t
        li = 1
        lf = self.pend
        p0 = self.p
        pf = self.pend

        #create lammps object
        lmp = ph.create_object(self.cores, self.simfolder, self.options["md"]["timestep"])

        lmp.command("echo              log")
        lmp.command("variable          li equal %f"%li)
        lmp.command("variable          lf equal %f"%lf)

        #read in conf
        conf = os.path.join(self.simfolder, "conf.dump")
        lmp = ph.read_dump(lmp, conf, species=self.options["nelements"])

        #set up potential
        lmp = ph.set_potential(lmp, self.options)

        #remap the box to get the correct pressure
        lmp = ph.remap_box(lmp, self.lx, self.ly, self.lz)

        #equilibrate first
        lmp.command("fix               1 all npt temp %f %f %f %s %f %f %f"%(t0, t0, self.options["md"]["tdamp"],
                                        self.iso, p0, p0, self.options["md"]["pdamp"]))
        lmp.command("run               %d"%self.options["md"]["te"])
        lmp.command("unfix             1")


        #now scale system to final temp, thereby recording enerfy at every step
        lmp.command("variable          step    equal step")
        lmp.command("variable          dU      equal pe/atoms")
        lmp.command("variable          lambda equal ramp(${li},${lf})")

        lmp.command("fix               f2 all npt temp %f %f %f %s %f %f %f"%(t0, t0, self.options["md"]["tdamp"],
                                        self.iso, p0, pf, self.options["md"]["pdamp"]))
        lmp.command("fix               f3 all print 1 \"${dU} $(press) $(vol) ${lambda}\" screen no file forward_%d.dat"%iteration)
        lmp.command("run               %d"%self.options["md"]["ts"])

        lmp.command("unfix             f2")
        lmp.command("unfix             f3")


        lmp.command("fix               1 all npt temp %f %f %f %s %f %f %f"%(t0, t0, self.options["md"]["tdamp"],
                                        self.iso, pf, pf, self.options["md"]["pdamp"]))
        lmp.command("run               %d"%self.options["md"]["te"])
        lmp.command("unfix             1")

        #start reverse loop
        lmp.command("variable          lambda equal ramp(${lf},${li})")

        lmp.command("fix               f2 all npt temp %f %f %f %s %f %f %f"%(t0, t0, self.options["md"]["tdamp"],
                                        self.iso, pf, p0, self.options["md"]["pdamp"]))
        lmp.command("fix               f3 all print 1 \"${dU} $(press) $(vol) ${lambda}\" screen no file backward_%d.dat"%iteration)
        lmp.command("run               %d"%self.options["md"]["ts"])

    
    
    def integrate_pressure_scaling(self, return_values=False):
        """
        Perform integration after reversible scaling
        
        Parameters
        ----------
        scale_energy : bool, optional
            If True, scale the energy during reversible scaling. 
        return_values : bool, optional
            If True, return integrated values
        
        Returns
        -------
        res : list of lists of shape 1x3
            Only returned if `return_values` is True.
        """
        res = integrate_ps(self.simfolder, self.fe, self.natoms,
            nsims=self.nsims, return_values=return_values)

        if return_values:
            return res