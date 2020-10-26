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
    lqd.run_averaging()
    assert os.path.exists("traj.dat") == True

    lqd.gather_average_data()

    assert lqd.natoms == 500
    assert np.abs(lqd.rho - 0.05283928377836297) < 1E-2
    assert np.abs(lqd.eps - 4.308666631) < 1E-2

    lqd.run_integration(iteration=1)
    assert os.path.exists("forward_1.dat") == True
    assert os.path.exists("backward_1.dat") == True

    lqd.run_integration(iteration=2)
    assert os.path.exists("forward_2.dat") == True
    assert os.path.exists("backward_2.dat") == True 

    lqd.run_integration(iteration=3)
    assert os.path.exists("forward_3.dat") == True
    assert os.path.exists("backward_3.dat") == True

    lqd.thermodynamic_integration()
    assert np.abs(lqd.fe - -2.2387655392976304) < 1E-2

    lqd.submit_report()
    assert os.path.exists("report.yaml") == True