from pytint.integrators import *
from scipy.integrate import cumtrapz
import numpy as np

def rev_scale(folder, reps=5, t0=None, f0=None, scale_pe=False):
    """
    Perform the reversible scaling method
    """
    ws = []

    for i in tqdm(range(1, reps+1)):
        fdx, flambda = np.loadtxt(os.path.join(folder, "forward_%d.dat"%i), unpack=True, comments="#")
        bdx, blambda = np.loadtxt(os.path.join(folder, "backward_%d.dat"%i), unpack=True, comments="#")

        if scale_pe:
            fdx /= flambda
            bdx /= blambda
        wf = cumtrapz(fdx, flambda,initial=0)
        wb = cumtrapz(bdx[::-1], blambda[::-1],initial=0)
        w = (wf + wb) / (2*flambda)
        ws.append(w)
 
    wmean = np.mean(ws, axis=0)
    werr = np.std(ws, axis=0)
    temp = t0/flambda
    f = f0/flambda + 1.5*kb*temp*np.log(flambda) + wmean
    return temp, f, werr
