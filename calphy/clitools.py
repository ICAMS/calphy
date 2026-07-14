import os
import numpy as np
import shutil
import argparse as ap
import subprocess
import yaml
import time
import datetime

from calphy.input import read_inputfile, _convert_legacy_inputfile
from calphy.phase_diagram import prepare_inputs_for_phase_diagram


def convert_legacy_inputfile():
    arg = ap.ArgumentParser()
    arg.add_argument("-i", "--input", required=True, type=str,
    help="name of the input file")
    arg.add_argument("-s", "--split", required=False, type=bool, 
    help="split each calculation into new file.", default=True)
    arg.add_argument("-o", "--output", required=False, type=str, 
    help="output file string, calculations will be named <outputstring>.*.yaml", 
        default='input')    
    args = vars(arg.parse_args())
    calculations = _convert_legacy_inputfile(args['input'], return_calcs=True)
    outputstr = args['output']

    if args['split']:
        #now we have to write this out to file
        for count, calc in enumerate(calculations):
            data = {}
            data['calculations'] = [calc]
            outfile = ".".join([outputstr, str(count+1), 'yaml'])
            with open(outfile, 'w') as fout:
                yaml.safe_dump(data, fout)
    else:
        data = {}        
        data['calculations'] = calculations
        outfile = ".".join([outputstr, 'yaml'])
        with open(outfile, 'w') as fout:
            yaml.safe_dump(data, fout)


def phase_diagram():
    arg = ap.ArgumentParser()
    arg.add_argument("-i", "--input", required=True, type=str,
    help="name of the input file")
    args = vars(arg.parse_args())
    prepare_inputs_for_phase_diagram(args['input'])
