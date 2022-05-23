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

import subprocess as sub
import os
import stat

class Local:
    """
    Local submission script
    """
    def __init__(self, options, cores=1, directory=os.getcwd()):
        self.queueoptions = {"scheduler": "local",
                             "jobname": "tis",
                             "walltime": None,
                             "queuename": None,
                             "memory": None,
                             "cores": cores,
                             "hint": None,
                             "directory": directory,
                             "options": [],
                             "commands": [],
                             "modules": [],
                             "header": "#!/bin/bash"

                            }
        for (key, val) in options.items():
            if key in self.queueoptions.keys():
                if val is not None:
                    self.queueoptions[key] = val
        self.maincommand = ""

    def write_script(self, outfile):
        """
        Write the script file
        """
        jobout = ".".join([outfile, "out"])
        joberr = ".".join([outfile, "err"])

        with open(outfile, "w") as fout:
            fout.write(self.queueoptions["header"])
            fout.write("\n")

            #now write modules
            for module in self.queueoptions["modules"]:
                fout.write("module load %s\n"   %module)

            #now finally commands
            for command in self.queueoptions["commands"]:
                fout.write("%s\n"   %command)
            fout.write("%s > %s 2> %s\n"   %(self.maincommand, jobout, joberr))
        self.script = outfile

    def submit(self):
        """
        Submit the job
        """
        st = os.stat(self.script)
        os.chmod(self.script, st.st_mode | stat.S_IEXEC)
        cmd = [self.script]
        proc = sub.Popen(cmd, stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)
        return proc

class SLURM:
    """
    Slurm class for writing submission script
    """
    def __init__(self, options, cores=1, directory=os.getcwd()):
        """
        Create class
        """
        self.queueoptions = {"scheduler": "slurm",
                             "jobname": "tis",
                             "walltime": "23:59:00",
                             "queuename": "shorttime",
                             "memory": "3GB",
                             "cores": cores,
                             "hint": "nomultithread",
                             "directory": directory,
                             "options": [],
                             "commands": [ "uss=$(whoami)",
                                           "find /dev/shm/ -user $uss -type f -mmin +30 -delete",
                                         ],
                             "modules": [],
                             "header": "#!/bin/bash"

                            }
        for (key, val) in options.items():
            if key in self.queueoptions.keys():
                if val is not None:
                    self.queueoptions[key] = val
        self.maincommand = ""


    def write_script(self, outfile):
        """
        Write the script file
        """
        jobout = ".".join([outfile, "out"])
        joberr = ".".join([outfile, "err"])

        with open(outfile, "w") as fout:
            fout.write(self.queueoptions["header"])
            fout.write("\n")

            #write the main header options
            fout.write("#SBATCH --job-name=%s\n" %self.queueoptions["jobname"])
            fout.write("#SBATCH --time=%s\n"     %self.queueoptions["walltime"])
            fout.write("#SBATCH --partition=%s\n"%self.queueoptions["queuename"])
            fout.write("#SBATCH --ntasks=%s\n"   %str(self.queueoptions["cores"]))
            fout.write("#SBATCH --mem-per-cpu=%s\n"%self.queueoptions["memory"])
            fout.write("#SBATCH --hint=%s\n"     %self.queueoptions["hint"])
            fout.write("#SBATCH --chdir=%s\n"    %self.queueoptions["directory"])

            #now write extra options
            for option in self.queueoptions["options"]:
                fout.write("#SBATCH %s\n"   %option)

            #now write modules
            for module in self.queueoptions["modules"]:
                fout.write("module load %s\n"   %module)

            #now finally commands
            for command in self.queueoptions["commands"]:
                fout.write("%s\n"   %command)
            fout.write("%s > %s 2> %s\n"   %(self.maincommand, jobout, joberr))

        self.script = outfile

    def submit(self):
        """
        Submit the job
        """
        cmd = ['sbatch', self.script]
        proc = sub.Popen(cmd, stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)
        return proc



class SGE:
    """
    Slurm class for writing submission script
    """
    def __init__(self, options, cores=1, directory=os.getcwd()):
        """
        Create class
        """
        self.queueoptions = {"scheduler": "sge",
                             "jobname": "tis",
                             "walltime": "23:59:00",
                             "queuename": None,
                             "memory": "3GB",
                             "system": "smp",
                             "commands": [],
                             "modules": [], 
                             "options": ["-j y",
                                         "-R y",
                                         "-P ams.p",
                                        ],
                             "cores": cores,
                             "hint": None,
                             "directory": directory,
                             "header": "#!/bin/bash"
                            }
        for (key, val) in options.items():
            if key in self.queueoptions.keys():
                if val is not None:
                    self.queueoptions[key] = val
        self.maincommand = ""

    def write_script(self, outfile):
        """
        Write the script file
        """
        with open(outfile, "w") as fout:
            fout.write(self.queueoptions["header"])
            fout.write("\n")

            #write the main header options
            fout.write("#$ -N %s\n" %self.queueoptions["jobname"])
            fout.write("#$ -l h_rt=%s\n"     %self.queueoptions["walltime"])
            fout.write("#$ -l qname=%s\n"%self.queueoptions["queuename"])
            fout.write("#$ -pe %s %s\n"   %( self.queueoptions["system"], str(self.queueoptions["cores"])))
            fout.write("#$ -l h_vmem=%s\n"%self.queueoptions["memory"])
            fout.write("#$ -cwd %s\n"    %self.queueoptions["directory"])

            #now write extra options
            for option in self.queueoptions["options"]:
                fout.write("#$ %s\n"   %option)

            #now write modules
            for module in self.queueoptions["modules"]:
                fout.write("module load %s\n"   %module)

            #now finally commands
            for command in self.queueoptions["commands"]:
                fout.write("%s\n"   %command)

            fout.write("%s > %s 2> %s\n"   %(self.maincommand, jobout, joberr))
        
        self.script = outfile   


    def submit(self):
        """
        Submit the job
        """
        cmd = ['qsub', self.script]
        proc = sub.Popen(cmd, stdin=sub.PIPE,stdout=sub.PIPE,stderr=sub.PIPE)
        return proc
    
