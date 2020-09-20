import numpy as np

def write_script_solid(mdscriptfile, temp, press, k, lat, options):
	"""
	Write the md script file for submission of job for
	solid
	"""
	with open(mdscriptfile, 'w') as fout:
		fout.write("echo              log\n")

	  	fout.write("variable          T0 equal 0.7*%f\n"%temp)  # Initial temperature.
	  	fout.write("variable          te equal %d\n"%options["te"])   # Equilibration time (steps).
	  	fout.write("variable          ts equal %d\n"%options["ts"])  # Switching time (steps).
	  	fout.write("variable          k equal %f\n"%k)
	  	fout.write("variable          rand equal %d\n"%np.random.randint(0, 1000))


	#-------------------------- Atomic Setup --------------------------------------#  
		fout.write("units            metal\n")
		fout.write("boundary         p p p\n")
		fout.write("atom_style       atomic\n")

		fout.write("lattice          %s %f\n"%(options["lattice"], lat))
		fout.write("region           box block 0 %d 0 %d 0 %d\n"%(options["nx"], options["ny"], options["nz"]))
		fout.write("create_box       1 box\n")
		fout.write("create_atoms     1 box\n")


		fout.write("neigh_modify    every 1 delay 0 check yes once no\n")

		fout.write("pair_style       %s\n"%options["pair_style"])
		fout.write("pair_coeff       %s\n"%options["pair_coeff"])
		fout.write("mass             * %f\n"%options["mass"])

	#---------------------- Thermostat & Barostat ---------------------------------#
		fout.write("fix               f1 all nve\n")
		fout.write("fix               f2 all ti/spring 10.0 1000 1000 function 2\n")
		fout.write("fix               f3 all langevin %f %f %f %d zero yes\n"%(temp, temp, tdamp, np.random.randint(0, 10000)))
	
	#------------------ Computes, variables & modifications -----------------------#
		fout.write("compute           Tcm all temp/com\n")
		fout.write("fix_modify        f3 temp Tcm\n")

		fout.write("variable          step    equal step\n")
		fout.write("variable          dU      equal pe/atoms-f_f2/atoms\n")
		fout.write("variable          lambda  equal f_f2[1]\n")

		fout.write("variable          te_run  equal %d-1\n"%options["te"]) # Print correctly on fix print.
		fout.write("variable          ts_run  equal %d+1\n"%options["ts"]) # Print correctly on fix print.

	#------------------------- Thermo stuff ---------------------------------------#
		fout.write("thermo_style      custom step pe c_Tcm\n")
		fout.write("timestep          0.001\n")
		fout.write("thermo            10000\n")

	#------------------------- Running the Simulation -----------------------------#
		fout.write("velocity          all create $\{T0\} $\{rand\} mom yes rot yes dist gaussian\n")


		fout.write("variable          i loop 5\n")
		fout.write("label repetitions\n")
			fout.write("fix               f2 all ti/spring %f %d %d function 2\n"%(k, options["ts"], options["te"]))

			# Forward. 
			fout.write("run               $\{te_run\}\n")
			fout.write("fix               f4 all print 1 \"$\{dU\} $\{lambda\}\" screen no file forward_$i.dat \n")
			fout.write("run               $\{ts_run\}\n")
			fout.write("unfix             f4\n")

			# Backward. 
			fout.write("run               $\{te_run\}\n")
			fout.write("fix               f4 all print 1 \"$\{dU\} $\{lambda\}\" screen no file backward_$i.dat\n")
			fout.write("run               $\{\ts_run\}\n")
			fout.write("unfix             f4\n")

			fout.write("next i\n")
		fout.write("jump SELF repetitions\n")




