import pytest
from pytint.input import read_yamlfile
import pytint.lattice as pl
from pytint.liquid import Liquid
import os
import numpy as np

def test_liquid_averaging():
    options = read_yamlfile("tests/input.yaml")
    sol = Liquid(t=1000, p=0, l="fcc", apc=4,
                    alat=4.05, c=0.0, options=options, simfolder=os.getcwd(),
                    thigh=1500)


    sol.reversible_scaling(iteration=1)
    assert os.path.exists("forward_1.dat") == True
    assert os.path.exists("backward_1.dat") == True

    sol.reversible_scaling(iteration=2)
    assert os.path.exists("forward_2.dat") == True
    assert os.path.exists("backward_2.dat") == True 

    sol.reversible_scaling(iteration=3)
    assert os.path.exists("forward_3.dat") == True
    assert os.path.exists("backward_3.dat") == True

    sol.integrate_reversible_scaling(f0=-3.2)
    assert os.path.exists("reversible_scaling.dat") == True
