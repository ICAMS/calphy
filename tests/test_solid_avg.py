import pytest
from pytint.input import read_yamlfile
import pytint.lattice as pl
from pytint.solid import Solid
import os
import numpy as np

def test_solid_averaging():
    options = read_yamlfile("tests/input.yaml")
    sol = Solid(t=1300, p=0, l="fcc", apc=4,
                    alat=3.704, c=0.0, options=options, simfolder=os.getcwd(),
                    )
    sol.run_averaging()
    assert os.path.exists("msd.dat") == True

    sol.gather_average_data()

    assert sol.natoms == 256
    assert np.abs(sol.alat - 3.704) < 1E-1
    assert np.abs(sol.k - 1.31) < 1E-1

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
    assert np.abs(sol.fe - -4.00) < 1E-1

    sol.submit_report()
    assert os.path.exists("report.yaml") == True