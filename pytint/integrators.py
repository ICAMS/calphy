"""
Integrators for pytint
"""
import numpy as np
import scipy.constants as const
import sys
import math
import os
import warnings
from pytint.splines import splines, sum_spline1, sum_spline25, sum_spline50, sum_spline75, sum_spline100
from scipy.integrate import cumtrapz

#Constants
h = const.physical_constants["Planck constant in eV/Hz"][0]
hbar = h/(2*np.pi)
kb = const.physical_constants["Boltzmann constant in eV/K"][0]
kbJ = const.physical_constants["Boltzmann constant"][0]
Na = const.physical_constants["Avogadro constant"][0]
eV2J = const.eV

def get_ideal_gas_fe(temp, rho, natoms, mass, xa=1.0, xb=0.0):
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
    mass = mass/Na
    beta = (1/(kb*temp)) #units - eV
    #omega needs to be in m
    omega = (beta*h*h/(2*np.pi*mass))**0.5
    #convert omega
    omega = omega*(const.eV/1E-3)**0.5
    #the above is in metres - change to Angstrom
    omega = omega*1E10
    prefactor = 1/beta

    if xa> 0:
        ta = xa*(3*np.log(omega) + np.log(rho) -1 + np.log(xa))
    else:
        ta = 0
    if xb> 0:
        tb = xb*(3*np.log(omega) + np.log(rho) -1 + np.log(xb))
    else:
        tb = 0


    return prefactor*(ta + tb + (1/(2*natoms))*np.log(2*np.pi*natoms))


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


def get_einstein_crystal_fe(temp, natoms, mass, a, k, atoms_per_cell, cm_correction=True):
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

    atoms_per_cell : int
        number of atoms per unit cell

    cm_correction : bool, optional, default - True
        add the centre of mass correction to free energy

    Returns
    -------
    fe : float
        free energy of Einstein crystal

    """
    #convert mass first for single particle in kg
    mass = (mass/Na)*1E-3
    #convert k from ev/A2 to J/m2
    k = k*(eV2J/1E-20)
    omega = np.sqrt(k/mass)
    F_harm = -3*kb*temp*np.log((kb*temp)/(hbar*omega))

    #convert a to m
    a = a*1E-10
    vol = (natoms/atoms_per_cell) * (a**3)

    if cm_correction:
        F_cm = (kb*temp/natoms)*np.log((natoms/vol)*(2*np.pi*kbJ*temp/(natoms*k))**1.5)
        F_harm = (F_harm + F_cm)

    return F_harm

def integrate_path(fwdfilename, bkdfilename, usecols=(0, 1, 2), solid=True):
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
    fdui, fdur, flambda = np.loadtxt(fwdfilename, unpack=True, comments="#", usecols=usecols)
    bdui, bdur, blambda = np.loadtxt(bkdfilename, unpack=True, comments="#", usecols=usecols)

    #SOLID HAS NO ISSUES - NO SCALING NEEDED
    #THIS IS TEMPORARY
    #UFM ENERGY IS NOT SCALED IN LAMMPS-THIS IS WRONG! BUT UNTIL THEN, WE KEEP THIS
    if not solid:

        #now scale with lambda
        for i in range(len(fdui)):
            if flambda[i] !=0:
                fdui[i] = fdui[i]/flambda[i]
        for i in range(len(bdui)):
            if blambda[i] !=0:
                bdui[i] = bdui[i]/blambda[i]

        """
        for i in range(len(fdur)):
            if flambda[i] !=0:
                fdur[i] = fdur[i]/flambda[i]
        for i in range(len(bdur)):
            if blambda[i] !=0:
                bdur[i] = bdur[i]/blambda[i]
        """

    fdu = fdui - fdur
    bdu = bdui - bdur
    
    fw = np.trapz(fdu, flambda)
    bw = np.trapz(bdu, blambda)

    w = 0.5*(fw - bw)
    q = 0.5*(fw + bw)

    return w, q


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

def find_w(mainfolder, nsims=5, full=False, usecols=(0,1,2), solid=True):
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
        w, q = integrate_path(fwdfilename, bkdfilename, usecols=usecols, solid=solid)
        ws.append(w)
        qs.append(q)
        
    if full:
        return np.mean(ws), np.mean(qs), np.std(qs)
    else:
        return np.mean(ws)

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


def integrate_rs(simfolder, f0, t, nsims=5, scale_energy=False, 
    correction=True, natoms=None):
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
    for i in range(1, nsims+1):
        fdx, flambda = np.loadtxt(os.path.join(simfolder, "forward_%d.dat"%i), unpack=True, comments="#")
        bdx, blambda = np.loadtxt(os.path.join(simfolder, "backward_%d.dat"%i), unpack=True, comments="#")
        
        if scale_energy:
            fdx /= flambda
            bdx /= blambda
        wf = cumtrapz(fdx, flambda,initial=0)
        wb = cumtrapz(bdx[::-1], blambda[::-1],initial=0)
        w = (wf + wb) / (2*flambda)
        ws.append(w)
    
    if correction:
        cterm = kb*t*np.log(2*np.pi*natoms)/(2*natoms)
    else:
        cterm = 0

    wmean = np.mean(ws, axis=0)
    werr = np.std(ws, axis=0)
    temp = t/flambda

    f = (f0-cterm)/flambda + 1.5*kb*temp*np.log(flambda) + wmean
    warnings.warn("Correcting free energy: org: %f, corr: %f"%(f0, cterm))

    outfile = os.path.join(simfolder, "reversible_scaling.dat")
    np.savetxt(outfile, np.column_stack((temp, f, werr)))