import gi
from gi.repository import GLib
gi.require_version('Hkl', '5.0')
from gi.repository import Hkl
import numpy as np

def transform_motor2xyz(geom, coords):
    '''
    input::geom - the diffractometer geometry
        one of E4CV/E4CH/K4C/E6C/K6C
    input::coords - a numpy array of motor rotations
        E4CV (omega, phi, chi, tth)
        E4CH (omega, phi, chi, tth)
        K4C [...]
        E6C [...]
        K6C [...]

    output::xyz - a numpy array of 3D coordinates (x, y, z)
    '''
    xyz = np.
    return xyz

def sweep(wlen_start, wlen_end, wlen_incr, res):
    '''
    
    '''
    

if __name__ == "__main__":


