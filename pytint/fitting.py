import numpy as np
import pandas as pd

clqd = np.arange(0, 1.2, 0.2)
cfcc = np.arange(0, 0.12, 0.02)
cdia = 1 - cfcc
cfit = np.arange(0, 1.01, 0.01)

def solve_tangent(ps, pl, x1=1):
    fs = np.poly1d([ps(x1)])
    newpoly = pl - np.poly1d([ps(x1)]) - np.poly1d([1, -x1])*np.polyder(pl)
    r = newpoly.r
    realr = r[np.isreal(r)]
    return [x.real for x in realr if (0<=x<=1)]
    
def get_line(fl):
    """
    Create a line
    """
    y2 = fl[-1]
    y1 = fl[0]
    x2 = 1
    x2 = 0
    slope = (y2-y1)
    c = y1
    return slope, c

def normalise_fe(fl=None, ff=None, fd=None):
    #get line
    nfl = nff = nfd = None
    
    if fl is not None:
        slope, c = get_line(fl)
        flline = slope*clqd + c
    
        #now modify originals
        nfl = fl - flline
    
    #now the others
    if ff is not None:
        fccline = slope*cfcc + c
        nff = ff - fccline
    
    if fd is not None:
        dialine = slope*cdia + c
        nfd = fd - dialine
    
    return nfl, nff, nfd


"""
Legacy functions for conversion
"""

def convert_fe_file(infile, outfile, comments=""):
    df = pd.read_hdf(infile)
    tempdict = {}
    temps = df["temp"].values
    w = df["w"].values
    q = df["q"].values
    err = df["err"].values
    
    for count, temp in enumerate(temps):
        tempdict[str(temp)] = {"w": w[count], "q": q[count], "err": err[count]}
    
    tempdict["info"] = comments
    np.save(outfile, tempdict)
    return tempdict

def convert_solid_prep_data(infile, outfile, comments=""):
    df = pd.read_hdf(infile)
    tempdict = {}
    temps = df["temp"].values
    k = df["k"].values
    lat = df["lat"].values

    for count, temp in enumerate(temps):
        tempdict[str(temp)] = {"lat": lat[count], "k": k[count]}
    
    tempdict["info"] = comments
    np.save(outfile, tempdict)
    return tempdict

def convert_liquid_prep_data(infile, outfile, comments=""):
    df = pd.read_hdf(infile)
    tempdict = {}
    temps = df["temp"].values
    conc = df["conc"].values
    lat = df["lat"].values
    dens = df["dens"].values
    
    utemps = np.unique(temps)
    tempdict = {str(t):{} for t in utemps}
    
    for count, temp in enumerate(temps):
        tempdict[str(temp)][str(np.round(conc[count], decimals=1))] = {"lat": lat[count], "rho": dens[count]}
    
    tempdict["info"] = comments
    np.save(outfile, tempdict)
    return tempdict

def convert_data_file(infile, outfile, liquid=False, comments=""):
    if liquid:
        dd = convert_liquid_prep_data(infile, outfile, comments=comments)
    else:
        dd = convert_solid_prep_data(infile, outfile, comments=comments)
    return dd
    
