import pytest
from calphy.input import read_yamlfile
import calphy.lattice as pl
from calphy.liquid import Liquid
import os
import numpy as np

def test_liquid_averaging():
    options = read_yamlfile(os.path.join(os.getcwd(), "tests/input.yaml"))
    lqd = Liquid(options=options, kernel=0, simfolder=os.getcwd())
    lqd.run_averaging()
    assert os.path.exists("traj.melt") == True

    assert lqd.natoms == 256
    assert np.abs(lqd.rho - 0.075) < 1E-2

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
    assert np.abs(lqd.fe - -4.17) < 0.5

    lqd.submit_report()
    assert os.path.exists("report.yaml") == True