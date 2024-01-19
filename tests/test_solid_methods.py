import pytest                                                                                                        
from calphy.input import read_inputfile                                                                              
from calphy.solid import Solid
from calphy.liquid import Liquid                                                                                       
import os                                                                                                            
import numpy as np                                                                                                   
                                                                                                                     
def test_solid_averaging_berendsen_finite_pressure(): 
    calculations = read_inputfile(os.path.join(os.getcwd(), "tests/inp2.yaml"))                                     
    sol = Solid(calculation=calculations[0], simfolder=os.getcwd())
    assert sol.calc.equilibration_control == "nose-hoover"                                                  

def test_solid_averaging_nose_hoover_finite_pressure(): 
    calculations = read_inputfile(os.path.join(os.getcwd(), "tests/inp3.yaml"))                                     
    sol = Liquid(calculation=calculations[0], simfolder=os.getcwd())
    assert sol.calc.equilibration_control == "nose-hoover"                                                  

