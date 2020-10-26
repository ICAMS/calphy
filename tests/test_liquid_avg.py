import pytest
from pytint.input import read_yamlfile
import pytint.lattice as pl
from pytint.liquid import Liquid
import os
import numpy as np

def test_liquid_averaging():
    options = read_yamlfile("tests/input.yaml")
    lqd = Liquid(t=1000, p=0, l="fcc", apc=4,
                    alat=4.05, c=0.0, options=options, simfolder=os.getcwd(),
                    thigh=1500)
    lqd.write_average_script()
    assert os.path.exists("traj.dat") == True
