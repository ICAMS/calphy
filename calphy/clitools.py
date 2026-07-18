import os
import numpy as np
import shutil
import argparse as ap
import subprocess
import yaml
import time
import datetime

from calphy.input import read_inputfile
from calphy.phase_diagram import prepare_inputs_for_phase_diagram


def phase_diagram():
    arg = ap.ArgumentParser()
    arg.add_argument("-i", "--input", required=True, type=str,
    help="name of the input file")
    args = vars(arg.parse_args())
    prepare_inputs_for_phase_diagram(args['input'])
