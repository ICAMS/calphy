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
    if os.path.exists(file):
        with open(file) as file:
            input = yaml.load(file, Loader=yaml.FullLoader)
    else:
        raise FileNotFoundError('%s input file not found'% file)

    return input