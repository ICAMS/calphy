"""
Module to handle input files
"""
import os
import yaml

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

    #we need to set up def options
    options = {}
    
    options["main"] = {
        "temperature": None,
        "pressure": 0,
        "element": None,
        "lattice": None,
        "concentration": 0,
        "nsims": 3
    }

    options["md"] = {
        "timestep": 0.001,
        "pair_style": None,
        "pair_coeff": None,
        "nsmall": 10000,
        "nevery": 10,
        "nrepeat": 10,
        "nfreq": 100,
        "mass": None,
        "tdamp": 0.1,
        "pdamp": 0.1,
        "nx": 7,
        "ny": 7,
        "nz": 7,
        "te": 25000,
        "ts": 50000
    }

    options["queue"] = {
        "scheduler": "local",
        "cores": 1,
        "jobname": "ti",
        "walltime": "23:50:00",
        "queuename": None,
        "memory": "3GB"
        "commands": None,
        "modules": None,
        "options": None
    }


    if os.path.exists(file):
        with open(file) as file:
            indata = yaml.load(file, Loader=yaml.FullLoader)
            #now read keys
            for okey in options.keys():
                if okey in indata.keys():
                    for key, val in indata[okey].items():
                        options["okey"]["key"] = indata["okey"]["key"] 
    else:
        raise FileNotFoundError('%s input file not found'% file)

    return options