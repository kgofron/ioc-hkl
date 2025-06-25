import numpy as np


def businglevy(wlen, lattice, refl1, refl2):

    return UB

def simplex(wlen, lattice, refl_list):
    return lattice, UB

def motor2hkl(wlen, lattice, omega, tth):
    return h, k, l

def hkl2motor(wlen, lattice, h, k, l):
    return omega, tth

def TAX_motor2hkl(wlen, lattice, omega, tth, mono, analyzer):
    tthplusmono = np.sin(mono...) ... # input to diffractometer is tth including mono angle, same thing for analyzer?
    motor2hkl(wlen, lattice, omega, tthplusmono)
    return h, k, l

def TAX_hkl2motor(wlen, lattice, h, k, l):
    return omega, tth, mono, analyzer
