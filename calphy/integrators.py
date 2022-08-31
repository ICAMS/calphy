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
from scipy.integrate import cumtrapz
from tqdm import tqdm

#Constants
h = const.physical_constants["Planck constant in eV/Hz"][0]
hbar = h/(2*np.pi)
kb = const.physical_constants["Boltzmann constant in eV/K"][0]
kbJ = const.physical_constants["Boltzmann constant"][0]
Na = const.physical_constants["Avogadro constant"][0]
eV2J = const.eV

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


def get_einstein_crystal_fe(temp, natoms, mass, vol, k, concentration, cm_correction=True):
    """
    Get the free energy of einstein crystal

    Parameters
    ----------
    temp : temperature, float
        units - K

    natoms : int
        no of atoms in the system

    mass : float
        units - g/mol

    a : lattice constant, float
        units - Angstrom

    k : spring constant, float
        units - eV/Angstrom^2

    cm_correction : bool, optional, default - True
        add the centre of mass correction to free energy

    Returns
    -------
    fe : float
        free energy of Einstein crystal

    """
    #convert mass first for single particle in kg
    mass = (np.array(mass)/Na)*1E-3

    #convert k from ev/A2 to J/m2
    k = np.array(k)*(eV2J/1E-20)
    omega = np.sqrt(k/mass)

    #convert a to m3
    vol = vol*1E-30

    F_harm = 0
    F_cm = 0

    for count, om in enumerate(omega):
        F_harm += concentration[count]*np.log((hbar*om)/(kb*temp))
        if cm_correction:
            F_cm += np.log((natoms*concentration[count]/vol)*(2*np.pi*kbJ*temp/(natoms*concentration[count]*k[count]))**1.5)
            #F_cm = 0
    F_harm = 3*kb*temp*F_harm
    F_cm = (kb*temp/natoms)*F_cm

    F_harm = F_harm + F_cm

    return F_harm

def integrate_path(fwdfilename, bkdfilename, nelements=1, concentration=[1,], usecols=(0, 1, 2), solid=True,
    alchemy=False, composition_integration=False):
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
    if solid:
        fdui = np.loadtxt(fwdfilename, unpack=True, comments="#", usecols=(0,))
        bdui = np.loadtxt(bkdfilename, unpack=True, comments="#", usecols=(0,))

        fdur = np.zeros(len(fdui))
        bdur = np.zeros(len(bdui))

        for i in range(nelements):
            fdur += concentration[i]*np.loadtxt(fwdfilename, unpack=True, comments="#", usecols=(i+1,))
            bdur += concentration[i]*np.loadtxt(bkdfilename, unpack=True, comments="#", usecols=(i+1,))

        flambda = np.loadtxt(fwdfilename, unpack=True, comments="#", usecols=(nelements+1,))
        blambda = np.loadtxt(bkdfilename, unpack=True, comments="#", usecols=(nelements+1,))

    else:
        fdui, fdur, flambda = np.loadtxt(fwdfilename, unpack=True, comments="#", usecols=usecols)
        bdui, bdur, blambda = np.loadtxt(bkdfilename, unpack=True, comments="#", usecols=usecols)

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

def find_w(mainfolder, nelements=1, concentration=[1,], nsims=5, full=False, usecols=(0,1,2), solid=True,
    alchemy=False, composition_integration=False):
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

    for i in range(nsims):
        fwdfilestring = 'forward_%d.dat' % (i+1)
        fwdfilename = os.path.join(mainfolder,fwdfilestring)
        bkdfilestring = 'backward_%d.dat' % (i+1)
        bkdfilename = os.path.join(mainfolder,bkdfilestring)
        w, q, flambda = integrate_path(fwdfilename, bkdfilename, nelements=nelements, concentration=concentration, usecols=usecols, solid=solid,
            alchemy=alchemy, composition_integration=composition_integration)
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


def integrate_mts(folder1, folder2, natoms1, natoms2, 
    pressure, temperature, nsims=5, scale_energy=True, 
                  full=False, stdscale=0):

    ws = []
    for i in tqdm(range(nsims)):
        fsu, fsp, fsv, fsl = np.loadtxt(os.path.join(folder1, "forward_%d.dat"%(i+1)), unpack=True)
        bsu, bsp, bsv, bsl = np.loadtxt(os.path.join(folder1, "backward_%d.dat"%(i+1)), unpack=True)
        flu, flp, flv, fll = np.loadtxt(os.path.join(folder2, "forward_%d.dat"%(i+1)), unpack=True)
        blu, blp, blv, bll = np.loadtxt(os.path.join(folder2, "backward_%d.dat"%(i+1)), unpack=True)

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

        wf = cumtrapz(fx, fsl, initial=0)
        wb = cumtrapz(bx[::-1], bsl[::-1], initial=0)

        w = (wf + wb) / (2*fsl)
        q = (wf - wb) / (2*fsl)
        
        if stdscale > 0:
            peak  = (w-np.roll(w, shift=-1))
            stdcut = np.std(peak[:-1])
            for i in range(len(peak[:-1])):
                if not (-stdscale*stdcut < peak[i] < stdscale*stdcut):
                    w[i:] = w[i:] + peak[i]
        ws.append(w)
    
    #return ws
    wmean = np.mean(ws, axis=0)
    werr = np.std(ws, axis=0)
    xp = fsl*(pressure/(10000*160.21766208)) - wmean
    temp = temperature/fsl

    #convert back
    xp = xp*160.217*10000

    #return the values
    if not full:
        return xp, temp
    else:
        return xp, temp, werr


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


def integrate_rs(simfolder, f0, t, natoms, p=0, nsims=5, 
    scale_energy=False, return_values=False):
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
        mcorarr += 1.5*kb*temperature*(flambda*(target_counts[i]/natoms)*np.log(ref_mass/target_masses[i]))
        mcorsum += 1.5*kb*temperature*(1.0*(target_counts[i]/natoms)*np.log(ref_mass/target_masses[i]))

    return mcorarr, mcorsum