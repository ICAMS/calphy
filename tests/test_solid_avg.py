import pytest
from calphy.input import read_yamlfile
import calphy.lattice as pl
from calphy.solid import Solid
import os
import numpy as np

def test_solid_averaging():
    options = read_yamlfile(os.path.join(os.getcwd(), "tests/input.yaml"))
    sol = Solid(options=options, kernel=0, simfolder=os.getcwd())
    sol.run_averaging()
    assert os.path.exists("msd.dat") == True

    assert sol.natoms == 256
    assert np.abs(sol.alat - 3.704) < 1E-1
    assert sol.k[0] > 0

    sol.run_integration(iteration=1)
    assert os.path.exists("forward_1.dat") == True
    assert os.path.exists("backward_1.dat") == True

    sol.run_integration(iteration=2)
    assert os.path.exists("forward_2.dat") == True
    assert os.path.exists("backward_2.dat") == True 

    sol.run_integration(iteration=3)
    assert os.path.exists("forward_3.dat") == True
    assert os.path.exists("backward_3.dat") == True

    sol.thermodynamic_integration()
    assert np.abs(sol.fe - -4.00) < 0.5

    sol.submit_report()
    assert os.path.exists("report.yaml") == True