"""
Module to handle input files
"""
import os
import yaml

def create_calculation_dict():
    """
    Create a template for calculation dict

    Parameters
    ----------
    None

    Returns
    -------
    dict: dict
        template dict
    """
    cdict = {}
    cdict["mode"] = "fe"
    cdict["temperature"] = None
    cdict["pressure"] = None
    cdict["lattice"] = None
    cdict["repeat"] = [7, 7, 7]
    cdict["state"] = None
    cdict["nsims"] = 1
    return cdict

def read_yamlfile(file):
    """
    Read a yaml input file
    Parameters
    ----------
    inputfile: string
        name of inout yaml file
    Returns
    -------
    dict: dict
        the read input dict options
    """
    #there are three blocks - main, md and queue
    #main block has subblocks of calculations

    #we need to set up def options
    options = {}
    
    #main dictionary
    options["main"] = {
        "element": None,
        "mass": 1.00,
    }

    #create a list for calculations
    options["main"]["calculations"] = []

    #options for md
    options["md"] = {
        #pair elements
        "pair_style": None,
        "pair_coeff": None,
        #time related properties
        "timestep": 0.001,
        "nsmall": 10000,
        "nevery": 10,
        "nrepeat": 10,
        "ncycles": 100,
        #ensemble properties
        "tdamp": 0.1,
        "pdamp": 0.1,
        #eqbr and switching time
        "te": 25000,
        "ts": 50000
    }

    #queue properties
    options["queue"] = {
        "scheduler": "local",
        "cores": 1,
        "jobname": "ti",
        "walltime": "23:50:00",
        "queuename": None,
        "memory": "3GB",
        "commands": None,
        "modules": None,
        "options": None
    }

    #convergence factors that can be set if required
    options["conv"] = {
        "alat_tol": 0.0002,
        "k_tol": 0.01,
        "solid_frac": 0.7,
        "liquid_frac": 0.05
    }


    if os.path.exists(file):
        with open(file) as file:
            indata = yaml.load(file, Loader=yaml.FullLoader)
            #now read keys
            for okey in options.keys():
                if okey in indata.keys():
                    for key, val in indata[okey].items():
                        options[okey][key] = indata[okey][key] 
    else:
        raise FileNotFoundError('%s input file not found'% file)

    return options