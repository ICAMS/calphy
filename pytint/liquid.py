"""
Contains methods for liquid
"""

import numpy as np
import yaml
from pytint.integrators import *
import pyscal.core as pc
import pyscal.traj_process as ptp

class Liquid:
    """
    Liquid class
    """
    def __init__(self, t=None, p=None, l=None, apc=None,
                    alat=None, c=None, options=None, simfolder=None,
                    thigh=None):
        """
        Set up class
        """
        self.t = t
        self.p = p
        self.l = l
        self.apc = apc
        self.alat = alat
        self.c = c
        self.options = options
        self.simfolder = simfolder
        self.thigh = thigh

    def write_average_script(self, mdscriptfile):
        """
        Write averagin script for solid
        """
        with open(mdscriptfile, 'w') as fout:
            fout.write("units            metal\n")
            fout.write("boundary         p p p\n")
            fout.write("atom_style       atomic\n")
            fout.write("timestep         %f\n"%self.options["md"]["timestep"])

            fout.write("lattice          %s %f\n"%(self.l, self.alat))
            fout.write("region           box block 0 %d 0 %d 0 %d\n"%(self.options["md"]["nx"], self.options["md"]["ny"], self.options["md"]["nz"]))
            fout.write("create_box       1 box\n")
            fout.write("create_atoms     1 box\n")

            fout.write("pair_style       %s\n"%self.options["md"]["pair_style"])
            fout.write("pair_coeff       %s\n"%self.options["md"]["pair_coeff"])
            fout.write("mass             * %f\n"%self.options["md"]["mass"])

            fout.write("velocity         all create %f %d\n"%(self.t, np.random.randint(0, 10000)))
            fout.write("fix              1 all npt temp %f %f %f iso %f %f %f\n"%(self.thigh, self.thigh, self.options["md"]["tdamp"], 
                                                    0, self.p, self.options["md"]["pdamp"]))
            fout.write("thermo_style     custom step pe press vol etotal temp\n")
            fout.write("thermo           10\n")
            fout.write("run              %d\n"%int(self.options["md"]["nsmall"])) 
            fout.write("unfix            1\n")

            fout.write("velocity         all create %f %d\n"%(self.thigh, np.random.randint(0, 10000)))
            fout.write("fix              1 all npt temp %f %f %f iso %f %f %f\n"%(self.thigh, self.thigh, self.options["md"]["tdamp"], 
                                                    self.p, self.p, self.options["md"]["pdamp"]))
            fout.write("dump             2 all custom %d traj.dat id type mass x y z vx vy vz\n"%(int(self.options["md"]["nsmall"])/10))
            fout.write("run              %d\n"%int(self.options["md"]["nsmall"])) 
            fout.write("unfix            1\n")

            fout.write("velocity         all create %f %d\n"%(self.t, np.random.randint(0, 10000)))
            fout.write("fix              1 all npt temp %f %f %f iso %f %f %f\n"%(self.t, self.t, self.options["md"]["tdamp"], self.p, self.p, self.options["md"]["pdamp"]))
            fout.write("run              %d\n"%int(self.options["md"]["nsmall"])) 
            fout.write("fix              2 all print 10 \"$(step) $(press) $(vol) $(temp)\" file avg.dat\n")
            fout.write("run              %d\n"%int(self.options["md"]["nlarge"]))

    def gather_average_data(self):
        """
        Gather average data
        """
        avgfile = os.path.join(self.simfolder, "avg.dat")
        vol = np.loadtxt(avgfile, usecols=(2,), unpack=True)
        avgvol = np.mean(vol[-100:])
        ncells = self.options["md"]["nx"]*self.options["md"]["ny"]*self.options["md"]["nz"]
        self.natoms = ncells*self.apc
        self.rho = self.natoms/avgvol
        #WARNING: hard coded ufm parameter
        self.eps = self.t*50.0*kb


    def process_traj(self):
        """
        Copy conf
        """
        trajfile = os.path.join(self.simfolder, "traj.dat")
        files = ptp.split_trajectory(trajfile)
        conf = os.path.join(self.simfolder, "conf.dump")

        sys = pc.System()
        sys.read_inputfile(files[-1], customkeys=["vx", "vy", "vz", "mass"])
        sys.to_file(conf, customkeys=["vx", "vy", "vz", "mass"])

        os.remove(trajfile)
        for file in files:
            os.remove(file)


    def write_integrate_script(self, mdscriptfile):
        """
        Write TI integrate script
        """

        self.process_traj()

        with open(mdscriptfile, 'w') as fout:
            fout.write("label RESTART\n")

            fout.write("variable        rnd      equal   round(random(0,999999,%d))\n"%np.random.randint(0, 10000))


            fout.write("variable        dt       equal   %f\n"%self.options["md"]["timestep"])             # Timestep (ps).

            # Adiabatic switching parameters.
            fout.write("variable        li       equal   1.0\n")               # Initial lambda.
            fout.write("variable        lf       equal   0.0\n")               # Final lambda.
            fout.write("variable        N_sim    loop    %d\n"%self.options["main"]["nsims"])                # Number of independent simulations.
            #------------------------------------------------------------------------------------------------------#


            ########################################     Atomic setup     ##########################################
            # Defines the style of atoms, units and boundary conditions.
            fout.write("units            metal\n")
            fout.write("boundary         p p p\n")
            fout.write("atom_style       atomic\n")
            fout.write("timestep         %f\n"%self.options["md"]["timestep"])

            # Read atoms positions, velocities and box parameters.
            fout.write("lattice          %s %f\n"%(self.l, self.alat))
            fout.write("region           box block 0 %d 0 %d 0 %d\n"%(self.options["md"]["nx"], self.options["md"]["ny"], self.options["md"]["nz"]))
            fout.write("create_box       1 box\n")

            conf = os.path.join(self.simfolder, "conf.dump")
            fout.write("read_dump        %s 0 x y z vx vy vz scaled no box yes add keep\n"%conf)

            fout.write("neigh_modify    delay 0\n")

            # Define MEAM and UF potentials parameters.
            fout.write("pair_style       hybrid/overlay %s ufm 7.5\n"%self.options["md"]["pair_style"])
            
            #modify pair style
            pc =  self.options["md"]["pair_coeff"]
            pcraw = pc.split()
            #now add style
            pcnew = " ".join([*pcraw[:2], *[self.options["md"]["pair_style"],], *pcraw[2:]])

            fout.write("pair_coeff       %s\n"%pcnew)
            fout.write("pair_coeff       * * ufm %f 1.5\n"%self.eps) 
            fout.write("mass             * %f\n"%self.options["md"]["mass"])

            #------------------------------------------------------------------------------------------------------#


            ################################     Fixes, computes and constraints     ###############################
            # Integrator & thermostat.
            fout.write("fix             f1 all nve\n")                              
            fout.write("fix             f2 all langevin %f %f %f ${rnd}\n"%(self.t, self.t, self.options["md"]["tdamp"]))
            fout.write("variable        rnd equal round(random(0,999999,0))\n")

            # Compute the potential energy of each pair style.
            fout.write("compute         c1 all pair %s\n"%self.options["md"]["pair_style"])
            fout.write("compute         c2 all pair ufm\n")
            #------------------------------------------------------------------------------------------------------#


            ##########################################     Output setup     ########################################
            # Output variables.
            fout.write("variable        step equal step\n")
            fout.write("variable        dU equal (c_c1-c_c2)/atoms\n")             # Driving-force obtained from NEHI procedure.

            # Thermo output.
            fout.write("thermo_style    custom step v_dU\n")
            fout.write("thermo          1000\n")
            #------------------------------------------------------------------------------------------------------#


            ##########################################     Run simulation     ######################################
            # Turn UF potential off (completely) to equilibrate the Sw potential.
            fout.write("variable        zero equal 0\n")
            fout.write("fix             f0 all adapt 0 pair ufm scale * * v_zero\n")
            fout.write("run             0\n")
            fout.write("unfix           f0\n")

            # Equilibrate the fluid interacting by Sw potential and switch to UF potential (Forward realization).
            fout.write("run             %d\n"%self.options["md"]["te"])

            fout.write("print           \"${dU} ${li}\" file forward_${N_sim}.dat\n")
            fout.write("variable        lambda_sw equal ramp(${li},${lf})\n")                 # Linear lambda protocol from 1 to 0.
            fout.write("fix             f3 all adapt 1 pair %s scale * * v_lambda_sw\n"%self.options["md"]["pair_style"])
            fout.write("variable        lambda_ufm equal ramp(${lf},${li})\n")                  # Linear lambda protocol from 0 to 1.
            fout.write("fix             f4 all adapt 1 pair ufm scale * * v_lambda_ufm\n")
            fout.write("fix             f5 all print 1 \"${dU} ${lambda_sw}\" screen no append forward_${N_sim}.dat\n")
            fout.write("run             %d\n"%self.options["md"]["ts"])

            fout.write("unfix           f3\n")
            fout.write("unfix           f4\n")
            fout.write("unfix           f5\n")

            # Equilibrate the fluid interacting by UF potential and switch to sw potential (Backward realization).
            fout.write("run             %d\n"%self.options["md"]["te"])

            fout.write("print           \"${dU} ${lf}\" file backward_${N_sim}.dat\n")
            fout.write("variable        lambda_sw equal ramp(${lf},${li})\n")                 # Linear lambda protocol from 0 to 1.
            fout.write("fix             f3 all adapt 1 pair %s scale * * v_lambda_sw\n"%self.options["md"]["pair_style"])
            fout.write("variable        lambda_ufm equal ramp(${li},${lf})\n")                  # Linear lambda protocol from 1 to 0.
            fout.write("fix             f4 all adapt 1 pair ufm scale * * v_lambda_ufm\n")
            fout.write("fix             f5 all print 1 \"${dU} ${lambda_sw}\" screen no append backward_${N_sim}.dat\n")
            fout.write("run             %d\n"%self.options["md"]["ts"])

            fout.write("unfix           f3\n")
            fout.write("unfix           f4\n")
            fout.write("unfix           f5\n")
            #------------------------------------------------------------------------------------------------------#


            ##########################################     Loop procedure     ######################################
            fout.write("next N_sim\n")
            fout.write("clear\n")
            fout.write("jump %s RESTART\n"%mdscriptfile)
            #------------------------------------------------------------------------------------------------------#

    
    def thermodynamic_integration(self):
        """
        Perform thermodynamic integration
        """
        w, q, qerr = find_w(self.simfolder, nsims=self.options["main"]["nsims"], 
            full=True)  
        #WARNING: hardcoded UFM parameters           
        f1 = get_uhlenbeck_ford_fe(self.t, 
            self.rho, 50, 1.5)
        f2 = get_ideal_gas_fe(self.t, self.rho, 
            self.natoms, self.options["md"]["mass"], xa=(1-self.c), 
            xb=self.c)
        self.fe = f2 + f1 - w
        self.ferr = qerr

    def submit_report(self):
        report = {}
        report["temperature"] = int(self.t)
        report["pressure"] = float(self.p)
        report["lattice"] = str(self.l)
        report["concentration"] = float(self.c)
        report["rho"] = float(self.rho)
        report["fe"] = float(self.fe)
        report["fe_err"] = float(self.ferr)

        reportfile = os.path.join(self.simfolder, "report.yaml")
        with open(reportfile, 'w') as f:
            yaml.dump(report, f)

    def write_rs_script(self, mdscriptfile):
        """
        Write TI integrate script
        """
        #rev scale needs tstart and tend; here self.t will be start
        #tend will be the final temp
        
        t0 = self.t
        tf = self.options["main"]["temperature"][-1]
        li = 1
        lf = t0/tf

        with open(mdscriptfile, 'w') as fout:
            fout.write("echo              log\n")

            fout.write("variable          T0 equal %f\n"%t0)  # Initial temperature.
            fout.write("variable          te equal %d\n"%self.options["md"]["te"])   # Equilibration time (steps).
            fout.write("variable          ts equal %d\n"%self.options["md"]["ts"])  # Switching time (steps).
            fout.write("variable          li equal %f\n"%li)
            fout.write("variable          lf equal %f\n"%lf)
            fout.write("variable          rand equal %d\n"%np.random.randint(0, 1000))


        #-------------------------- Atomic Setup --------------------------------------#  
            fout.write("units            metal\n")
            fout.write("boundary         p p p\n")
            fout.write("atom_style       atomic\n")

            fout.write("lattice          %s %f\n"%(self.l, self.alat))
            fout.write("region           box block 0 %d 0 %d 0 %d\n"%(self.options["md"]["nx"], self.options["md"]["ny"], self.options["md"]["nz"]))
            fout.write("create_box       1 box\n")
            conf = os.path.join(self.simfolder, "conf.dump")
            fout.write("read_dump        %s 0 x y z vx vy vz scaled no box yes add keep\n"%conf)


            fout.write("neigh_modify    every 1 delay 0 check yes once no\n")

            fout.write("pair_style       %s\n"%self.options["md"]["pair_style"])
            fout.write("pair_coeff       %s\n"%self.options["md"]["pair_coeff"])
            fout.write("mass             * %f\n"%self.options["md"]["mass"])

        #---------------------- Thermostat & Barostat ---------------------------------#
            fout.write("fix               f1 all nph aniso %f %f %f\n"%(self.p, self.p, self.options["md"]["pdamp"]))
            fout.write("fix               f2 all langevin ${T0} ${T0} %f %d zero yes\n"%(self.options["md"]["tdamp"], np.random.randint(0, 10000)))
            fout.write("run               ${te}\n")
            fout.write("unfix             f1\n")
            fout.write("unfix             f2\n")

            fout.write("variable         xcm equal xcm(all,x)\n")
            fout.write("variable         ycm equal xcm(all,y)\n")
            fout.write("variable         zcm equal xcm(all,z)\n")
            
            fout.write("fix              f1 all nph aniso %f %f %f fixedpoint ${xcm} ${ycm} ${zcm}\n"%(self.p, self.p, self.options["md"]["pdamp"]))
            fout.write("fix              f2 all langevin ${T0} ${T0} %f %d zero yes\n"%(self.options["md"]["tdamp"], np.random.randint(0, 10000)))
            
        #------------------ Computes, variables & modifications -----------------------#
            fout.write("compute           tcm all temp/com\n")
            fout.write("fix_modify        f1 temp tcm\n")
            fout.write("fix_modify        f2 temp tcm\n")

            fout.write("variable          step    equal step\n")
            fout.write("variable          dU      equal c_thermo_pe/atoms\n")
            fout.write("variable          te_run  equal ${te}-1\n")
            fout.write("variable          ts_run  equal ${ts}+1\n")
            fout.write("thermo_style      custom step pe c_tcm\n")
            fout.write("timestep          %f\n"%self.options["md"]["timestep"])
            fout.write("thermo            10000\n")
            

            fout.write("velocity          all create ${T0} ${rand} mom yes rot yes dist gaussian\n")   
            fout.write("variable          i loop %d\n"%self.options["main"]["nsims"])
            fout.write("label repetitions\n")
            fout.write("    run               ${te}\n")
            fout.write("    variable          lambda equal ramp(${li},${lf})\n")

            #we need to similar to liquid here

            fout.write("    fix               f3 all adapt 1 pair %s scale * * v_lambda\n"%self.options["md"]["pair_style"])
            fout.write("    fix               f4 all print 1 \"${dU} ${lambda}\" screen no file forward_$i.dat 
\n")
            fout.write("    run               ${ts}\n")
            fout.write("    unfix             f3\n")
            fout.write("    unfix             f4\n")
            fout.write("    run               ${te}\n")
            fout.write("    variable          lambda equal ramp(${lf},${li})\n")
            fout.write("    fix               f3 all adapt 1 pair %s scale * * v_lambda\n"%self.options["md"]["pair_style"])
            fout.write("    fix               f4 all print 1 \"${dU} ${lambda}\" screen no file backward_$i.dat\n")
            fout.write("    run               ${ts}\n")
            fout.write("    unfix             f3\n")
            fout.write("    unfix             f4\n")
            
            fout.write("next i\n")
            fout.write("jump SELF repetitions\n")