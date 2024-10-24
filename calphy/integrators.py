"""
calphy: a Python library and command line interface for automated free
energy calculations.

Copyright 2021  (c) Sarath Menon^1, Yury Lysogorskiy^2, Ralf Drautz^2
^1: Max Planck Institut für Eisenforschung, Dusseldorf, Germany 
^2: Ruhr-University Bochum, Bochum, Germany

calphy is published and distributed under the Academic Software License v1.0 (ASL). 
calphy is distributed in the hope that it will be useful for non-commercial academic research, 
but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
calphy API is published and distributed under the BSD 3-Clause "New" or "Revised" License
See the LICENSE FILE for more details. 

More information about the program can be found in:
Menon, Sarath, Yury Lysogorskiy, Jutta Rogal, and Ralf Drautz.
“Automated Free Energy Calculation from Atomistic Simulations.” Physical Review Materials 5(10), 2021
DOI: 10.1103/PhysRevMaterials.5.103801

For more information contact:
sarath.menon@ruhr-uni-bochum.de/yury.lysogorskiy@icams.rub.de
"""

import numpy as np
import scipy.constants as const
import sys
import math
import os
import warnings
from calphy.splines import splines, sum_spline1, sum_spline25, sum_spline50, sum_spline75, sum_spline100
try:
    from scipy.integrate import cumtrapz
except ImportError:
    from scipy.integrate import cumulative_trapezoid as cumtrapz
from tqdm import tqdm
import pyscal3.core as pc
from ase.io import read

#Constants
h = const.physical_constants["Planck constant in eV/Hz"][0]
hJ = const.physical_constants["Planck constant"][0]
hbar = h/(2*np.pi)
kb = const.physical_constants["Boltzmann constant in eV/K"][0]
kbJ = const.physical_constants["Boltzmann constant"][0]
Na = const.physical_constants["Avogadro constant"][0]
eV2J = const.eV
J2eV = 6.242E18

#--------------------------------------------------------------------
#             TI PATH INTEGRATION ROUTINES
#--------------------------------------------------------------------

def integrate_path(calc,
    fwdfilename, 
    bkdfilename,  
    solid=True,
    composition_integration=False):
    """
    Get a filename with columns du and dlambda and integrate

    Parameters
    ----------
    fwdfilename: string
        name of fwd integration file

    bkdfilename: string
        name of bkd integration file

    usecols : list
        column numbers to be used from input file

    Returns
    -------
    w : float
        irreversible work in switching the system

    q : float
        heat dissipation during switching of system
    """
    natoms = np.array([calc._element_dict[x]['count'] for x in calc.element])
    concentration = np.array([calc._element_dict[x]['composition'] for x in calc.element])
    fdata = np.loadtxt(fwdfilename, unpack=True, comments="#")
    bdata = np.loadtxt(bkdfilename, unpack=True, comments="#")

    if solid:
        fdui = fdata[0]
        bdui = bdata[0]

        fdur = np.zeros(len(fdui))
        bdur = np.zeros(len(bdui))

        for i in range(calc.n_elements):
            if natoms[i] > 0:
                fdur += concentration[i]*fdata[i+1]/natoms[i]
                bdur += concentration[i]*bdata[i+1]/natoms[i]

        flambda = fdata[calc.n_elements+1]
        blambda = bdata[calc.n_elements+1]

    else:
        fdui = fdata[0]
        bdui = bdata[0]

        fdur = fdata[1]
        bdur = bdata[1]

        flambda = fdata[2]
        blambda = bdata[2]


    fdu = fdui - fdur
    bdu = bdui - bdur
    
    if composition_integration:
        fw = cumtrapz(fdu, flambda, initial=0)
        bw = cumtrapz(bdu, blambda, initial=0)
    else:
        fw = np.trapz(fdu, flambda)
        bw = np.trapz(bdu, blambda)

    w = 0.5*(fw - bw)
    q = 0.5*(fw + bw)

    return w, q, flambda


def find_w(mainfolder,
    calc, 
    full=False, 
    solid=True,
    composition_integration=False):
    """
    Integrate the irreversible work and dissipation for independent simulations

    Parameters
    ----------
    mainfolder: string
        main simulation folder

    nsims : int, optional
        number of independent simulations, default 5

    full : bool, optional
        If True return error values, default False

    usecols : tuple, optional
        Columns to read in from data file. Default (0, 1)

    Returns
    -------
    ws : float
        average irreversible work

    qs : float
        average energy dissipation, only returned if full is True

    err : float
        Error in free energy, only returned if full is True
    """
    ws = []
    qs = []

    for i in range(calc.n_iterations):
        fwdfilestring = 'forward_%d.dat' % (i+1)
        fwdfilename = os.path.join(mainfolder,fwdfilestring)
        
        bkdfilestring = 'backward_%d.dat' % (i+1)
        bkdfilename = os.path.join(mainfolder,bkdfilestring)
        
        w, q, flambda = integrate_path(calc,
            fwdfilename, 
            bkdfilename, 
            solid=solid,
            composition_integration=composition_integration)
        
        ws.append(w)
        qs.append(q)
    
    if composition_integration:
        wsmean = np.mean(ws, axis=0)
        qsmean = np.mean(qs, axis=0)
        wsstd = np.std(ws, axis=0)
        return wsmean, qsmean, wsstd, flambda

    wsmean = np.mean(ws)
    qsmean = np.mean(qs)
    wsstd = np.std(ws)

    if full:
        return wsmean, qsmean, wsstd
    else:
        return wsmean


def integrate_rs(simfolder, f0, t, 
    natoms, p=0, nsims=5, 
    scale_energy=False, 
    return_values=False):
    """
    Carry out the reversible scaling integration

    Parameters
    ----------
    simfolder : string
        main simulation folder

    f0 : float
        initial free energy for integration

    t : float
        initial temperature

    nsims : int, optional
        number of independent switching

    scale_energy: bool, optional
        if True, scale energy with switching parameter

    Returns
    -------
    None

    Notes
    -----
    Writes the output in a file reversible_scaling.dat
    
    """
    ws = []
    p = p/(10000*160.21766208)
    
    for i in range(1, nsims+1):
        fdx, fp, fvol, flambda = np.loadtxt(os.path.join(simfolder, "ts.forward_%d.dat"%i), unpack=True, comments="#")
        bdx, bp, bvol, blambda = np.loadtxt(os.path.join(simfolder, "ts.backward_%d.dat"%i), unpack=True, comments="#")
        
        if scale_energy:
            fdx /= flambda
            bdx /= blambda

        #add pressure contribution
        fvol = fvol/natoms
        bvol = bvol/natoms
        fdx = fdx + p*fvol
        bdx = bdx + p*bvol
        
        wf = cumtrapz(fdx, flambda,initial=0)
        wb = cumtrapz(bdx[::-1], blambda[::-1],initial=0)
        w = (wf + wb) / (2*flambda)
        ws.append(w)
    

    wmean = np.mean(ws, axis=0)
    werr = np.std(ws, axis=0)
    temp = t/flambda

    f = f0/flambda + 1.5*kb*temp*np.log(flambda) + wmean

    if not return_values:
        outfile = os.path.join(simfolder, "temperature_sweep.dat")
        np.savetxt(outfile, np.column_stack((temp, f, werr)))
    else:
        return (temp, f, werr)


def integrate_ps(simfolder, f0, natoms, pi, pf, nsims=1, 
    return_values=False):
    """
    Carry out the reversible scaling integration
    
    Parameters
    ----------
    simfolder : string
        main simulation folder
    f0 : float
        initial free energy for integration
    nsims : int, optional
        number of independent switching
    
    Returns
    -------
    None
    
    Notes
    -----
    Writes the output in a file pressure_sweep.dat
    
    """

    ws = []

    for i in range(1, nsims+1):
        _, fp, fvol, _ = np.loadtxt(os.path.join(simfolder, "ps.forward_%d.dat"%i), unpack=True, comments="#")
        _, bp, bvol, _ = np.loadtxt(os.path.join(simfolder, "ps.backward_%d.dat"%i), unpack=True, comments="#")
        
        fvol = fvol/natoms
        bvol = bvol/natoms
        
        fp = fp/(10000*160.21766208)
        bp = bp/(10000*160.21766208)

        wf = cumtrapz(fvol, fp, initial=0)
        wb = cumtrapz(bvol[::-1], bp[::-1], initial=0)

        w = (wf + wb)/2
        ws.append(w)

    wmean = np.mean(ws, axis=0)
    werr = np.std(ws, axis=0)
    
    press = np.linspace(pi, pf, len(wmean))

    f = f0 + wmean

    if not return_values:
        outfile = os.path.join(simfolder, "pressure_sweep.dat")
        np.savetxt(outfile, np.column_stack((press, f, werr)))
    else:
        return (press, f, werr)


def integrate_mass(flambda, ref_mass, target_masses, target_counts,
    temperature, natoms):
    
    mcorarr = np.zeros(len(flambda))
    mcorsum = 0

    for i in range(len(target_masses)):
        mcorarr += 1.5*kb*temperature*(flambda*(target_counts[i]/natoms)*np.log(target_masses[i]/ref_mass))
        mcorsum += 1.5*kb*temperature*(1.0*(target_counts[i]/natoms)*np.log(target_masses[i]/ref_mass))

    return mcorarr, mcorsum

def remove_steps(w, stdscale):
    peak  = np.abs(w-np.roll(w, shift=-1))
    peak = np.where(peak> stdscale*np.std(peak), peak, 0)
    args = np.nonzero(peak)[0]
    print(f'No of peaks #{len(args)}')
    diff = [w[x]-w[x-1] if x in args else 0 for x in range(len(w[:-1]))]
    diff.append(0)
    cum_diff = np.cumsum(diff)
    w = w-cum_diff
    return w

def remove_peaks(w, stdscale):
    peak = np.minimum(np.abs(w-np.roll(w, shift=-1)), np.abs(np.roll(w, shift=1)-w))
    args = np.argwhere(peak > stdscale*np.std(peak))
    k = [(w[x-1]+w[x+1])/2 if x in args else w[x] for x in range(len(w[:-1]))]
    k.append(w[-1])
    return k

def integrate_dcc(folder1, folder2, nsims=1, scale_energy=True, 
                  full=False, stdscale=0.25, fit_order=None):
    """
    Integrate Dynamic Clausius-Clapeyron equation

    Parameters
    ----------
    folder1: string
        Calculation folder for calculation of the first phase

    folder2: string
        Calculation folder for calculation of the second phase

    nsims: int
        Number of iterations, default 1

    full: bool
        If True, return error, default False

    stdscale: float
        Function used to smooth the integrand. If a value is provided,
        values in the greater than stdscale*std_dev and less than -stdscale*std_dev
        are identified as steps/peaks and smoothened. Default 0.25

    fit_order: None
        If specified, a final fitting is done to remove kinks in the values
        Optional, default None.

    """

    #get number of atoms
    sys = pc.System(os.path.join(folder1, "traj.equilibration_stage1.dat"))
    natoms1 = int(sys.natoms)

    sys = pc.System(os.path.join(folder2, "traj.equilibration_stage1.dat"))
    natoms2 = int(sys.natoms)

    #get temp and pressure
    f1_raw = os.path.basename(folder1).split("-")
    pressure = float(f1_raw[-1])
    temperature = float(f1_raw[-2])

    #recheck for folder2
    f2_raw = os.path.basename(folder2).split("-")
    if pressure != float(f2_raw[-1]):
        raise ValueError("Pressure for both calculations are different!")
    if temperature != float(f2_raw[-2]):
        raise ValueError("Temperature for both calculations are different!")
    if pressure == 0:
        raise ValueError("A non-zero pressure is needed for dcc")

    ws = []
    for i in range(nsims):
        fsu, fsp, fsv, fsl = np.loadtxt(os.path.join(folder1, "ts.forward_%d.dat"%(i+1)), unpack=True)
        bsu, bsp, bsv, bsl = np.loadtxt(os.path.join(folder1, "ts.backward_%d.dat"%(i+1)), unpack=True)
        flu, flp, flv, fll = np.loadtxt(os.path.join(folder2, "ts.forward_%d.dat"%(i+1)), unpack=True)
        blu, blp, blv, bll = np.loadtxt(os.path.join(folder2, "ts.backward_%d.dat"%(i+1)), unpack=True)

        if scale_energy:
            fsu = fsu/fsl
            bsu = bsu/bsl
            flu = flu/fll
            blu = blu/bll

        #now convert the pressure units
        fsp = fsp/(10000*160.21766208)
        bsp = bsp/(10000*160.21766208)
        flp = flp/(10000*160.21766208)
        blp = blp/(10000*160.21766208)

        #scale volume per number of atoms
        fsv = fsv/natoms1
        bsv = bsv/natoms1
        flv = flv/natoms2
        blv = blv/natoms2

        #get the integrand
        fx = (fsu-flu)/(fsv-flv)
        bx = (bsu-blu)/(bsv-blv)

        if stdscale > 0:
            fx = remove_peaks(fx, stdscale=stdscale)
            bx = remove_peaks(bx, stdscale=stdscale)

        wf = cumtrapz(fx, fsl, initial=0)
        wb = cumtrapz(bx[::-1], bsl[::-1], initial=0)

        w = (wf + wb) / (2*fsl)
        q = (wf - wb) / (2*fsl)
        
        if stdscale > 0:
            w = remove_steps(w, stdscale=stdscale)
        
        ws.append(w)
    
    #return ws
    wmean = np.mean(ws, axis=0)
    werr = np.std(ws, axis=0)
    xp = fsl*(pressure/(10000*160.21766208)) - wmean
    temp = temperature/fsl

    #convert back
    xp = xp*160.217*10000

    #fit if needed
    if fit_order is not None:
        fit = np.polyfit(temp, xp, fit_order)
        xp = np.polyval(fit, temp)

    #return the values
    if not full:
        return xp, temp
    else:
        return xp, temp, werr

#--------------------------------------------------------------------
#             REF. STATE ROUTINES: SOLID
#--------------------------------------------------------------------

def get_einstein_crystal_fe(
    calc,
    vol, 
    k, 
    cm_correction=True,
    return_contributions=False):
    """
    Get the free energy of einstein crystal

    Parameters
    ----------
    calc : Calculation object
        contains all input parameters

    vol : float
        converged volume per atom

    k : spring constant, float
        units - eV/Angstrom^2

    cm_correction : bool, optional, default - True
        add the centre of mass correction to free energy

    return_contributions: bool, optional, default - True
        If True, return individual contributions to the reference free energy.

    Returns
    -------
    F_tot : float
        total free energy of reference crystal

    F_e : float
        Free energy of Einstein crystal without centre of mass correction. Only if `return_contributions` is True.

    F_cm : float
        centre of mass correction. Only if `return_contributions` is True.

    Notes
    -----
    The equations for free energy of Einstein crystal and centre of mass correction are from https://doi.org/10.1063/5.0044833.

    """
    #temperature
    temp = calc._temperature

    #natoms
    natoms = np.sum([calc._element_dict[x]['count'] for x in calc.element])

    #convert a to m3
    vol = vol*1E-30

    #whats the beta
    beta = (1/(kbJ*temp))   

    #create an array of mass
    mass = []
    for x in calc.element:
        for count in range(calc._element_dict[x]['count']):
            mass.append(calc._element_dict[x]['mass'])
    mass = np.array(mass)
    
    #convert mass to kg
    mass = (mass/Na)*1E-3

    #create an array of k as well
    karr = []
    for c, x in enumerate(calc.element):
        for count in range(calc._element_dict[x]['count']):
            karr.append(k[c])
    k = np.array(karr)
    #convert k from ev/A2 to J/m2
    k = k*(eV2J/1E-20)

    #fe of Einstein crystal
    Z_e = ((beta**2*k*hJ**2)/(4*np.pi**2*mass))**1.5
    F_e = np.log(Z_e)
    F_e = kb*temp*np.sum(F_e)/natoms #*J2eV #convert back to eV

    #now get the cm correction
    if cm_correction:
        mass_sum = np.sum(mass)
        mu = mass/mass_sum
        mu2_over_k = mu**2/k
        mu2_over_k_sum = np.sum(mu2_over_k) 
        prefactor = vol
        F_cm = np.log(prefactor*(beta/(2*np.pi*mu2_over_k_sum))**1.5)
        F_cm = kb*temp*F_cm/natoms #convert to eV
    else:
        F_cm = 0
    
    F_tot = F_e - F_cm
    if return_contributions:
        return F_e, -F_cm
    return F_tot

#--------------------------------------------------------------------
#             REF. STATE ROUTINES: LIQUID
#--------------------------------------------------------------------

def get_ideal_gas_fe(temp, rho, natoms, mass, concentration):
    """
    Get the free energy of an single/binary ideal gas

    Parameters
    ----------
    temp : temperature, float
        the reference temperature in K

    rho : number density, float
        units - no of atoms/ angstrom^3

    natoms: int
        total number of atoms

    mass : atomic mass, float
        units - g/mol

    xa : concentration of species a, float, optional
        default 1

    xb : concentration of species b, float, optional
        default 0

    Returns
    -------
    fe : float
        free energy/atom of ideal gas system

    """
    #find mass of one particle
    mass = np.array(mass)/Na
    beta = (1/(kb*temp)) #units - eV

    #omega needs to be in m
    omega = (beta*h*h/(2*np.pi*mass))**0.5
    #convert omega
    omega = omega*(const.eV/1E-3)**0.5
    #the above is in metres - change to Angstrom
    omega = omega*1E10
    prefactor = 1/beta

    fe = 0
    for count, conc in enumerate(concentration):
        if concentration[count] > 0:
            fe += conc*(3*np.log(omega[count]) + np.log(rho) -1 + np.log(conc))

    #return prefactor*(ta + tb + (1/(2*natoms))*np.log(2*np.pi*natoms))
    return prefactor*fe


def get_uhlenbeck_ford_fe(temp, rho, p, sigma):
    """
    Get the excess free energy of Uhlenbeck-Ford model

    Parameters
    ----------
    temp : temperature, float
        units - K

    rho : density, float
        units - no of atoms/ angstrom^3

    p : uf scale, float
    sigma : uf length scale, float

    Returns
    -------
    fe : float
        excess free energy/atom of uf system
    """
    x = (0.5*(np.pi*sigma*sigma)**1.5)*rho
    _, fe = find_fe(p, x)
    beta = (1/(kb*temp))
    fe = fe/beta
    return fe

def press(x, coef):
    """
    Find pressure of system

    Parameters
    ----------
    x : float
        x value for UF system
    
    coef : list of floats
        coefficients

    Returns
    -------
    result : float, optional
        result pressure

    """
    result = coef[0]*(x**3) + coef[1]*(x**2) + coef[2]*x + coef[3]
    return result

def fe(x, coef, sum_spline, index):
    """
    Fe inbuilt method
    """
    if x < 0.0025:
        result = coef[0]*(x**2)/2.0 + coef[1]*x
        return result

    elif x < 0.1:
        if x*10000%25 == 0:
            return sum_spline[index-1]
        else:
            x_0 = 0.0025*int(x*400)
    
    elif x < 1:
        if x*1000%25 == 0:
            return sum_spline[index-1]
        else:
            x_0 = 0.025*int(x*40)

    elif x < 4:
        if x*100%10 == 0:
            return sum_spline[index-1]
        else:
            x_0 = 0.1*int(x*10)
    else:
        return sum_spline[index]

    result =  sum_spline[index-1] + coef[0]*(x**2.0 - x_0**2.0)/2.0 + coef[1]*(x - x_0) + (coef[2] - 1.0)*math.log(x/x_0) - coef[3]*(1.0/x - 1.0/x_0)
    return result


def find_fe(p, x):
    """
    Find free energy of UF system

    Parameters
    ----------
    x : float
        x value of system

    coef : list of floats
        Coefficients of system

    Returns
    -------
    fe : float
        free energy of UF system

    """
    if not p in splines:
        raise ValueError('Invalid p. Valid numbers are: 1, 25, 50, 75, and 100.')

    if (x <= 0.0) or (x > 4.0):
        raise ValueError('Invalid x. Valid numbers are 0.0 < x <= 4.0')

    table1 = splines[p]

    if x < 0.1:
        index =  0 + int(x*400)
    elif x < 1:
        index = 40 + int((x*40 - 4))
    elif x < 4:
        index = 76 + int((x*10 - 10))
    else:
        index = 105

    coef = table1[index]

    pressure = press(x, coef)

    if p==1:
        sum_spline = sum_spline1
    elif p==25:
        sum_spline = sum_spline25
    elif p==50:
        sum_spline = sum_spline50
    elif p==75:
        sum_spline = sum_spline75
    else:
        sum_spline = sum_spline100

    free_energy = fe(x,coef,sum_spline,index)

    return pressure, free_energy

#--------------------------------------------------------------------
#             PHASE DIAGRAM ROUTINES
#--------------------------------------------------------------------

def calculate_entropy_mix(conc):
    """
    Calculate the entropy of mixing

    Parameters
    ----------
    conc : float
        concentration

    Returns
    -------
    s: float
        entropy
    """
    s = -kb*(conc*np.log(conc) + (1-conc)*np.log(1-conc))
    return s


def calculate_fe_impurity(temp, natoms, fepure, feimpure):
    """
    Calculate energy change of mixing, imput energies are
    in eV/atom

    Parameters
    ----------
    temp : float
        temperature

    natoms : int
        number of atoms

    fepure : float
        free energy of pure phase

    feimpure : float
        free energy of impure phase

    Returns
    -------
    dg : float
        entropy of mixing
    """
    dg = feimpure*natoms - fepure*natoms + kb*temp*np.log(natoms)
    return dg
    
def calculate_fe_mix(temp, fepure, feimpure, concs, natoms=4000):
    """
    Calculate energy of mixing

    Parameters
    ----------
    temp : float
        temperature

    fepure : float
        free energy of the pure phase

    feimpure : float
        energy due to impurity

    concs : list of floats
        concentration array

    natoms : int
        number of atoms

    Returns
    -------
    fes : list of floats
        free energy with concentration

    """
    if concs[0] == 0:
        print("zero is autodone")
        concs = concs[1:]
    fes = [fepure]
    
    s = calculate_entropy_mix(concs)
    dg = calculate_fe_impurity(temp, natoms, fepure, feimpure)
    fe_conc = fepure + concs*dg - temp*s
    
    for f in fe_conc:
        fes.append(f)    
    return fes