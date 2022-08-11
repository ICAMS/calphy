import pytest                                                                                                        
from calphy.input import read_inputfile                                                                              
import calphy.lattice as pl                                                                                          
from calphy.solid import Solid                                                                                       
import os                                                                                                            
import numpy as np                                                                                                   
                                                                                                                     
def test_solid_averaging_nose_hoover():                                                                                          
    calculations = read_inputfile(os.path.join(os.getcwd(), "tests/inp1.yaml"))                                     
    sol = Solid(calculation=calculations[0], simfolder=os.getcwd())                                                  
    sol.run_averaging()                                                                                              
    assert os.path.exists("msd.dat") == True

def test_solid_averaging_berendsen_finite_pressure(): 
    calculations = read_inputfile(os.path.join(os.getcwd(), "tests/inp2.yaml"))                                     
    sol = Solid(calculation=calculations[0], simfolder=os.getcwd())                                                  
    sol.run_averaging()                                                                                              
    assert os.path.exists("msd.dat") == True

def test_solid_averaging_nose_hoover_finite_pressure(): 
    calculations = read_inputfile(os.path.join(os.getcwd(), "tests/inp3.yaml"))                                     
    sol = Solid(calculation=calculations[0], simfolder=os.getcwd())                                                  
    sol.run_averaging()                                                                                              
    assert os.path.exists("msd.dat") == True 

