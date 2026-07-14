try:
    from calphy.input import Calculation
    from calphy.liquid import Liquid
    from calphy.solid import Solid
    from calphy.alchemy import Alchemy
    from calphy.routines import MeltingTemp
except ImportError as _exc:
    # calphy.runner and calphy.errors must remain importable without a LAMMPS
    # backend (pylammpsmpi today).  Defer only the backend-dependent driver
    # imports when that backend is absent; re-raise anything else so real
    # import bugs still surface.
    if "pylammpsmpi" not in str(_exc) and "lammps" not in str(_exc):
        raise

__version__ = "1.8.4"

def addtest(a,b):
    return a+b
