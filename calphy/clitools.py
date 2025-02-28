import os
import numpy as np
import shutil
import argparse as ap
import subprocess
import yaml
import time
import datetime

from calphy.input import read_inputfile, load_job, save_job, _convert_legacy_inputfile
from calphy.liquid import Liquid
from calphy.solid import Solid
from calphy.alchemy import Alchemy
from calphy.phase_diagram import prepare_inputs_for_phase_diagram

def _generate_job(calc, simfolder):
    if calc.mode == "alchemy" or calc.mode == "composition_scaling":
        job = Alchemy(calculation=calc, simfolder=simfolder)
        return job
    else:
        if calc.reference_phase == "liquid":
            job = Liquid(calculation=calc, simfolder=simfolder)
            return job
        else:
            job = Solid(calculation=calc, simfolder=simfolder)
            return job


def run_averaging():
    arg = ap.ArgumentParser()
    arg.add_argument("-i", "--input", required=True, type=str,
    help="name of the input file")
    arg.add_argument("-k", "--kernel", required=True, type=int, 
    help="kernel number of the calculation to be run.")
    args = vars(arg.parse_args())
    kernel = args["kernel"]
    calculations = read_inputfile(args["input"])
    calc = calculations[kernel]

    simfolder = calc.create_folders()
    job = _generate_job(calc, simfolder)
    os.chdir(simfolder)
    

    job.run_averaging()
    save_job(job)


def process_averaging():
    arg = ap.ArgumentParser()
    arg.add_argument("-i", "--input", required=True, type=str,
    help="name of the input file")
    arg.add_argument("-k", "--kernel", required=True, type=int, 
    help="kernel number of the calculation to be run.")
    args = vars(arg.parse_args())
    kernel = args["kernel"]
    calculations = read_inputfile(args["input"])
    calc = calculations[kernel]

    job = load_job(calc.savefile)
    job.process_averaging_results()
    save_job(job)

def run_integration():
    arg = ap.ArgumentParser()
    arg.add_argument("-i", "--input", required=True, type=str,
    help="name of the input file")
    arg.add_argument("-k", "--kernel", required=True, type=int, 
    help="kernel number of the calculation to be run.")
    args = vars(arg.parse_args())
    kernel = args["kernel"]
    calculations = read_inputfile(args["input"])
    calc = calculations[kernel]

    job = load_job(calc.savefile)
    job.run_integration()
    save_job(job)

def process_integration():
    arg = ap.ArgumentParser()
    arg.add_argument("-i", "--input", required=True, type=str,
    help="name of the input file")
    arg.add_argument("-k", "--kernel", required=True, type=int, 
    help="kernel number of the calculation to be run.")
    args = vars(arg.parse_args())
    kernel = args["kernel"]
    calculations = read_inputfile(args["input"])
    calc = calculations[kernel]

    job = load_job(calc.savefile)
    job.thermodynamic_integration()
    job.submit_report()
    save_job(job)

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
    prepare_inputs_for_phase_diagram(args['input'])
