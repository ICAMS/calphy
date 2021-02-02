"""
Module to handle input files
"""
import os
import yaml
import warnings

def check_and_convert_to_list(data):
    """
    Check if the given item is a list, if not convert to a single item list

    Parameters
    ----------
    data : single value or list

    Returns
    -------
    data : list
    """
    if not isinstance(data, list):
        return [data]
    else:
        return data

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
    options["element"]: None
    options["mass"]: 1.00

    #create a list for calculations
    options["calculations"] = []

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

    #keys that need to be read in directly
    directkeys = ["md", "queue", "conv"]

    #now read the file
    if os.path.exists(file):
        with open(file) as file:
            indata = yaml.load(file, Loader=yaml.FullLoader)
    else:
        raise FileNotFoundError('%s input file not found'% file)


    #now read keys
    for okey in directkeys:
        if okey in indata.keys():
            for key, val in indata[okey].items():
                options[okey][key] = indata[okey][key] 

    options["element"] = check_and_convert_to_list(indata["element"])
    options["mass"] = check_and_convert_to_list(indata["mass"])

    if not len(options["element"]) == len(options["mass"]):
        raise ValueError("length of elements and mass should be same!")
    options["nelements"] = len(options["element"])

    #now we need to process calculation keys
    #loop over calculations
    if "calculations" in indata.keys():
        #if the key is present
        #Loop 0: over each calc block
        #Loop 1: over lattice
        #Loop 2: over pressure
        #Loop 3: over temperature if needed - depends on mode
        for calc in indata["calculations"]:
            #check and convert items to lists if needed
            lattice = check_and_convert_to_list(calc["lattice"])
            state = check_and_convert_to_list(calc["state"])
            pressure = check_and_convert_to_list(calc["pressure"])
            temperature = check_and_convert_to_list(calc["temperature"])
            mode = calc["mode"]
            

            #now start looping
            for i, lat in enumerate(lattice):
                for press in pressure:
                    if mode == "ts":
                        cdict = {}
                        cdict["mode"] = calc["mode"]
                        cdict["temperature"] = temperature[0]
                        cdict["pressure"] = press
                        cdict["lattice"] = lat
                        cdict["state"] = state[i]
                        cdict["temperature_stop"] = temperature[-1]
                        cdict["nelements"] = options["nelements"]
                        cdict["element"] = options["element"]

                        if options["nelements"] == 1:
                            concentration = [1,]
                        else:
                            concentration = check_and_convert_to_list(calc["concentration"])
                            if not len(concentration) == options["nelements"]:
                                raise ValueError("concentration has to be same length as elements")
                        cdict["concentration"] = concentration

                        #optional keys
                        if "repeat" in calc.keys():
                            cdict["repeat"] = calc["repeat"]
                            if not (cdict["repeat"][0] == cdict["repeat"][1] == cdict["repeat"][2]):
                                raise ValueError("For LAMMPS structure creation, use nx=ny=nz")
                        else:
                            cdict["repeat"] = [1, 1, 1]
                        if "nsims" in calc.keys():
                            cdict["nsims"] = calc["nsims"]
                        else:
                            cdict["nsims"] = 1
                        if "thigh" in calc.keys():
                            cdict["thigh"] = calc["thigh"]
                        else:
                            cdict["thigh"] = 2.0*cdict["temperature_stop"]

                        options["calculations"].append(cdict)
                    else:
                        for temp in temperature:
                            cdict = {}
                            cdict["mode"] = calc["mode"]
                            cdict["temperature"] = temp
                            cdict["pressure"] = press
                            cdict["lattice"] = lat
                            cdict["state"] = state[i]
                            cdict["temperature_stop"] = temp
                            cdict["nelements"] = options["nelements"]
                            cdict["element"] = options["element"]

                            if options["nelements"] == 1:
                                concentration = [1,]
                            else:
                                concentration = check_and_convert_to_list(calc["concentration"])
                                if not len(concentration) == options["nelements"]:
                                    raise ValueError("concentration has to be same length as elements")
                            cdict["concentration"] = concentration

                            #optional keys
                            if "repeat" in calc.keys():
                                cdict["repeat"] = calc["repeat"]
                                if not (cdict["repeat"][0] == cdict["repeat"][1] == cdict["repeat"][2]):
                                    raise ValueError("For LAMMPS structure creation, use nx=ny=nz")
                            else:
                                cdict["repeat"] = [1, 1, 1]
                            if "nsims" in calc.keys():
                                cdict["nsims"] = calc["nsims"]
                            else:
                                cdict["nsims"] = 1
                            if "thigh" in calc.keys():
                                cdict["thigh"] = calc["thigh"]
                            else:
                                cdict["thigh"] = 2.0*cdict["temperature_stop"]

                            options["calculations"].append(cdict)

    return options

def create_identifier(calc):
    """
    Generate an identifier

    Parameters
    ----------
    calc: dict
        a calculation dict

    Returns
    -------
    identistring: string
        unique identification string
    """
    #lattice processed
    ts = int(calc["temperature"])
    ps = int(calc["pressure"])

    l = calc["lattice"]
    l = l.split('/')
    l = l[-1]
    
    #print(calc.keys())
    prefix = calc["mode"]

    identistring = "-".join([prefix, l, str(ts), str(ps)])
    return identistring