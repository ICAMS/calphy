import numpy as np

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