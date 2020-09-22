import numpy as np

def write_script_solid(mdscriptfile, temp, k, lat, options):
    """
    Write the md script file for submission of job for
    solid
    """
    with open(mdscriptfile, 'w') as fout:
        fout.write("echo              log\n")

        fout.write("variable          T0 equal 0.7*%f\n"%temp)  # Initial temperature.
        fout.write("variable          te equal %d\n"%options["md"]["te"])   # Equilibration time (steps).
        fout.write("variable          ts equal %d\n"%options["md"]["ts"])  # Switching time (steps).
        fout.write("variable          k equal %f\n"%k)
        fout.write("variable          rand equal %d\n"%np.random.randint(0, 1000))


    #-------------------------- Atomic Setup --------------------------------------#  
        fout.write("units            metal\n")
        fout.write("boundary         p p p\n")
        fout.write("atom_style       atomic\n")

        fout.write("lattice          %s %f\n"%(options["md"]["lattice"], lat))
        fout.write("region           box block 0 %d 0 %d 0 %d\n"%(options["md"]["nx"], options["md"]["ny"], options["md"]["nz"]))
        fout.write("create_box       1 box\n")
        fout.write("create_atoms     1 box\n")


        fout.write("neigh_modify    every 1 delay 0 check yes once no\n")

        fout.write("pair_style       %s\n"%options["md"]["pair_style"])
        fout.write("pair_coeff       %s\n"%options["md"]["pair_coeff"])
        fout.write("mass             * %f\n"%options["md"]["mass"])

    #---------------------- Thermostat & Barostat ---------------------------------#
        fout.write("fix               f1 all nve\n")
        fout.write("fix               f2 all ti/spring 10.0 1000 1000 function 2\n")
        fout.write("fix               f3 all langevin %f %f %f %d zero yes\n"%(temp, temp, options["md"]["tdamp"], 
                                        np.random.randint(0, 10000)))
    
    #------------------ Computes, variables & modifications -----------------------#
        fout.write("compute           Tcm all temp/com\n")
        fout.write("fix_modify        f3 temp Tcm\n")

        fout.write("variable          step    equal step\n")
        fout.write("variable          dU      equal pe/atoms-f_f2/atoms\n")
        fout.write("variable          lambda  equal f_f2[1]\n")

        fout.write("variable          te_run  equal %d-1\n"%options["md"]["te"]) # Print correctly on fix print.
        fout.write("variable          ts_run  equal %d+1\n"%options["md"]["ts"]) # Print correctly on fix print.

    #------------------------- Thermo stuff ---------------------------------------#
        fout.write("thermo_style      custom step pe c_Tcm\n")
        fout.write("timestep          0.001\n")
        fout.write("thermo            10000\n")

    #------------------------- Running the Simulation -----------------------------#
        fout.write("velocity          all create ${T0} ${rand} mom yes rot yes dist gaussian\n")


        fout.write("variable          i loop %d\n"%options["main"]["nsims"])
        fout.write("label repetitions\n")
        fout.write("fix               f2 all ti/spring %f %d %d function 2\n"%(k, options["md"]["ts"], options["md"]["te"]))

        # Forward. 
        fout.write("run               ${te_run}\n")
        fout.write("fix               f4 all print 1 \"${dU} ${lambda}\" screen no file forward_$i.dat \n")
        fout.write("run               ${ts_run}\n")
        fout.write("unfix             f4\n")

        # Backward. 
        fout.write("run               ${te_run}\n")
        fout.write("fix               f4 all print 1 \"${dU} ${lambda}\" screen no file backward_$i.dat\n")
        fout.write("run               ${ts_run}\n")
        fout.write("unfix             f4\n")

        fout.write("next i\n")
        fout.write("jump SELF repetitions\n")



def write_script_liquid(mdscriptfile, temp, epsilon, dumpfile, options):
    """
    Write the md script file for submission of job for
    solid
    """
    with open(mdscriptfile, 'w') as fout:
        fout.write("label RESTART\n")

        fout.write("variable        rnd      equal   %d\n"%np.random.randint(0, 10000))


        fout.write("variable        dt       equal   0.001\n")             # Timestep (ps).

        # Adiabatic switching parameters.
        fout.write("variable        li       equal   1.0\n")               # Initial lambda.
        fout.write("variable        lf       equal   0.0\n")               # Final lambda.
        fout.write("variable        N_sim    loop    %d\n"%options["main"]["nsims"])                # Number of independent simulations.
        #------------------------------------------------------------------------------------------------------#


        ########################################     Atomic setup     ##########################################
        # Defines the style of atoms, units and boundary conditions.
        fout.write("units            metal\n")
        fout.write("boundary         p p p\n")
        fout.write("atom_style       atomic\n")
        fout.write("timestep         %f\n"%options["md"]["timestep"])

        # Read atoms positions, velocities and box parameters.
        fout.write("lattice          %s %f\n"%(options["md"]["lattice"], options["md"]["lattice_constant"]))
        fout.write("region           box block 0 %d 0 %d 0 %d\n"%(options["md"]["nx"], options["md"]["ny"], options["md"]["nz"]))
        fout.write("create_box       1 box\n")
        fout.write("read_dump        %s 0 x y z vx vy vz scaled no box yes add keep\n"%dumpfile)

        fout.write("neigh_modify    delay 0\n")

        # Define MEAM and UF potentials parameters.
        fout.write("pair_style       hybrid/overlay %s ufm 7.5\n"%options["md"]["pair_style"])
        
        #modify pair style
        pc =  options["md"]["pair_coeff"]
        pcraw = pc.split()
        #now add style
        pcnew = " ".join([*pcraw[:2], *[options["md"]["pair_style"],], *pcraw[2:]])

        fout.write("pair_coeff       %s\n"%pcnew)
        fout.write("pair_coeff       * * ufm %f 1.5\n"%epsilon) 
        fout.write("mass             * %f\n"%options["md"]["mass"])

        #------------------------------------------------------------------------------------------------------#


        ################################     Fixes, computes and constraints     ###############################
        # Integrator & thermostat.
        fout.write("fix             f1 all nve\n")                              
        fout.write("fix             f2 all langevin %f %f %f %d zero yes\n"%(temp, temp, options["md"]["tdamp"], np.random.randint(0, 10000)))

        # Compute the potential energy of each pair style.
        fout.write("compute         c1 all pair %s\n"%options["md"]["pair_style"])
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
        fout.write("unfix           f0\n")

        # Equilibrate the fluid interacting by Sw potential and switch to UF potential (Forward realization).
        fout.write("run             %d\n"%options["md"]["te"])

        fout.write("print           \"${dU} ${li}\" file forward_${N_sim}.dat\n")
        fout.write("variable        lambda_sw equal ramp(${li},${lf})\n")                 # Linear lambda protocol from 1 to 0.
        fout.write("fix             f3 all adapt 1 pair %s scale * * v_lambda_sw\n"%options["md"]["pair_style"])
        fout.write("variable        lambda_ufm equal ramp(${lf},${li})\n")                  # Linear lambda protocol from 0 to 1.
        fout.write("fix             f4 all adapt 1 pair ufm scale * * v_lambda_ufm\n")
        fout.write("fix             f5 all print 1 \"${dU} ${lambda_sw}\" screen no append forward_${N_sim}.dat\n")
        fout.write("run             $%d\n"%options["md"]["ts"])

        fout.write("unfix           f3\n")
        fout.write("unfix           f4\n")
        fout.write("unfix           f5\n")

        # Equilibrate the fluid interacting by UF potential and switch to sw potential (Backward realization).
        fout.write("run             %d\n"%options["md"]["te"])

        fout.write("print           \"${dU} ${lf}\" file backward_${N_sim}.dat\n")
        fout.write("variable        lambda_sw equal ramp(${lf},${li})\n")                 # Linear lambda protocol from 0 to 1.
        fout.write("fix             f3 all adapt 1 pair %s scale * * v_lambda_sw\n"%options["md"]["pair_style"])
        fout.write("variable        lambda_ufm equal ramp(${li},${lf})\n")                  # Linear lambda protocol from 1 to 0.
        fout.write("fix             f4 all adapt 1 pair ufm scale * * v_lambda_ufm\n")
        fout.write("fix             f5 all print 1 \"${dU} ${lambda_sw}\" screen no append backward_${N_sim}.dat\n")
        fout.write("run             %d\n"%options["md"]["ts"])

        fout.write("unfix           f3\n")
        fout.write("unfix           f4\n")
        fout.write("unfix           f5\n")
        #------------------------------------------------------------------------------------------------------#


        ##########################################     Loop procedure     ######################################
        fout.write("next N_sim\n")
        fout.write("clear\n")
        fout.write("jump %s RESTART\n"%mdscriptfile)
        #------------------------------------------------------------------------------------------------------#

