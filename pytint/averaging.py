"""
Initial module which is used for calculation of
average values
"""
import numpy as np

def write_script_solid(mdscriptfile, temp, press, options):
	"""
	Write the md script file for submission of job for
	solid
	"""
	with open(mdscriptfile, 'w') as fout:
		fout.write("units            metal\n")
		fout.write("boundary         p p p\n")
		fout.write("atom_style       atomic\n")
		fout.write("timestep         %f\n"%options["timestep"])

		fout.write("lattice          %s %f\n"%(options["lattice"], options["lattice_constant"]))
		fout.write("region           box block 0 10 0 10 0 10\n")
		fout.write("create_box       1 box\n")
		fout.write("create_atoms     1 box\n")

		fout.write("pair_style       %s\n"%options["pair_style"])
		fout.write("pair_coeff       %s\n"%options["pair_coeff"])

		fout.write("mass             * %f\n"%options["mass"])

		fout.write("velocity         all create %f %d\n"%(temp, np.random.randint(0, 10000)))
		fout.write("fix              1 all npt temp %f %f %f iso %f %f %f\n"%(temp, temp, tdamp, press, press, pdamp))
		fout.write("thermo_style     custom step pe press vol etotal temp\n")
		fout.write("thermo           10\n")
		fout.write("run              %d\n"%int(options["nsmall"])) 

		fout.write("fix              2 all print 10 \"$(step) $(press) $(vol) $(temp)\" file avg.dat\n")
		fout.write("run              %d\n"%int(options["nlarge"]))

		#now run for msd
		fout.write("unfix            1\n")
		fout.write("unfix            2\n")

		fout.write("fix              3 all nvt temp %f %f %f\n"%(temp, temp, tdamp))
		fout.write("compute          1 all msd com yes\n")

		fout.write("variable         msd equal c_1[4]\n")

		fout.write("fix              4 all print 10 \"$(step) ${msd}\" file msd.dat\n")
		fout.write("run              %d\n"%int(options["nlarge"]))

def write_script_liquid(mdscriptfile, temp, thigh, press, options):
	with open(mdscriptfile, 'w') as fout:
		fout.write("units            metal\n")
		fout.write("boundary         p p p\n")
		fout.write("atom_style       atomic\n")
		fout.write("timestep         %f\n"%options["timestep"])

		fout.write("lattice          %s %f\n"%(options["lattice"], options["lattice_constant"]))
		fout.write("region           box block 0 10 0 10 0 10\n")
		fout.write("create_box       1 box\n")
		fout.write("create_atoms     1 box\n")

		fout.write("pair_style       %s\n"%options["pair_style"])
		fout.write("pair_coeff       %s\n"%options["pair_coeff"])
		fout.write("mass             * %f\n"%options["mass"])

		fout.write("velocity         all create %f %d\n"%(thigh, np.random.randint(0, 10000)))
		fout.write("fix              1 all npt temp %f %f %f iso %f %f %f\n"%(thigh, thigh, tdamp, press, press, pdamp))
		fout.write("thermo_style     custom step pe press vol etotal temp\n")
		fout.write("thermo           10\n")
		fout.write("dump             2 all custom 1000 %s id type mass x y z vx vy vz\n"%options["trajfile"])
		fout.write("run              %d\n"%int(options["nsmall"])) 
		fout.write("unfix            1\n")

		fout.write("velocity         all create %f %d\n"%(temp, np.random.randint(0, 10000)))
		fout.write("fix              1 all npt temp %f %f %f iso %f %f %f\n"%(temp, temp, tdamp, press, press, pdamp))
		fout.write("run              %d\n"%int(options["nsmall"])) 
		fout.write("fix              2 all print 10 \"$(step) $(press) $(vol) $(temp)\" file avg.dat\n")
		fout.write("run              %d\n"%int(options["nlarge"]))
