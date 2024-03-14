from calphy.input import read_inputfile
import shutil
from ase.io import read
import numpy as np
from tqdm.notebook import trange
import os
import re

try:
    from pyiron_atomistics import Project
    from pyiron_atomistics.atomistics.structure.atoms import ase_to_pyiron
except ImportError:
    raise ImportError('This feature needs pyiron_atomistics installed')

def create_job_from_inputfile(pr, inputfile, potential, kernel=None):
    """
    Create a pyiron job from calphy input file
    
    Parameters
    ----------
    pr: pyiron Project
        project to which the job is to be added
    
    inputfile: string
        calphy input file
    
    potential: string
        name of the potential as present in pyiron
    
    kernel: int, optional
        the index of the calculation to be read in. If None, all calculations are read in
    """
    calcs = read_inputfile(inputfile)
    
    if kernel is None:
        kernel = [x for x in range(len(calcs))]
    kernel = np.atleast_1d(kernel)
    
    for i in trange(len(kernel)):
        calc = calcs[kernel[i]]
        basedir = calc.create_identifier()
        basedir_path = os.path.join(os.path.dirname(inputfile), basedir)
        
        if os.path.exists(basedir_path):
            #create job and copy files
            try:
                #make sure that the report file exists
                reportfile = os.path.join(basedir_path, 'report.yaml')
                if os.path.exists(reportfile):
                    job = pr.create.job.Calphy(basedir.replace('-', '_'))
                    job._job_id = pr.db.add_item_dict(job.db_entry())
                    job.refresh_job_status()
                    shutil.copytree(basedir_path, job.working_directory, dirs_exist_ok=True)

                    #read in structure, assign potential
                    Z_of_type = dict([(count+1, calc._element_dict[element]['atomic_number']) for count, element in enumerate(calc.element)])
                    structure = read(calc.lattice, format='lammps-data', style='atomic', Z_of_type=Z_of_type)
                    job.structure = ase_to_pyiron(structure)
                    job.potential = potential
                    pr.db.item_update({"ChemicalFormula": job.structure.get_chemical_formula()}, job._job_id)

                    #collect output
                    job.input.mode = calc.mode
                    job.status.collect = True
                    job.collect_output()

                    #populate inputs
                    calcdict = calc.model_dump()
                    #temporary fix for comp scaling until its introduced in pyiron
                    del calcdict['composition_scaling']
                    job.input.update(calcdict)
                    job._create_calc()
                    job.to_hdf()
                    job.status.finished = True
                else:
                    print(f'parsing {basedir_path} failed, skipping')   
            except:
                #delete job
                pr.remove_job(basedir.replace('-', '_'))
                print(f'parsing {basedir_path} failed, skipping')
        else:
            print(f'could not find {basedir_path}, skipping')


def get_free_energy(job):
    return job["output/energy_free"]

def get_temperature(job):
    return job["output/temperature"]

def get_phase(job):
    raw = job.name.split('_')
    if raw[-3] == 'liquid':
        phase = 'liquid'
    else:
        phase = raw[1]
    return phase

def get_composition(job):
    chem = job.project.db.get_item_by_id(job.id)['chemicalformula']
    comp_split = re.split('(\d+)', chem)[:-1]
    if len(comp_split) == 2:
        if comp_split[0] == 'Al':
            return 0.00
        else:
            return 1.00
    else:
        return int(comp_split[3])/(int(comp_split[1])+int(comp_split[3]))