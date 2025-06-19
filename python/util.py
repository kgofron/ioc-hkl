import numpy as np
import math

e = 1.6021766300e-19 # [C]
h = 6.6260701500e-34 # [m^2*kg/s]
c = 299792458 # [m/s^2]
eV2J = 1.602176634e-19 # 1eV = ... J
m2A = 1e10 # meter to angstrom
m_neutron = 1.6749274710e-27 #[kg]
m_proton = 1.6726219200e-27 #[kg]
m_electron = 9.1093837e-27 #[kg]
velocity = 0. # sqrt(2E/m) #TODO

def sci2dec(sci):
    dec = np.round(np.float64(sci), 1)
    return dec

def energy2wavelength_neutron(energy):
    # E = p^2/2m = h^2/(2m*lambda^2) 
    # lambda = sqrt(h^2/2*E*m_neutron) = h*sqrt(1/2*E*m_neutron) * m2A #converting to A
    if energy > 0:
        energy = eV2J*energy/1000 # input is meV not eV
        coeff = h*m2A
        wavelength = coeff*math.sqrt(1/(2*energy*m_neutron))
    else:
        wavelength = 0
    return wavelength

def energy2wavelength_xrays(energy):
    # E [kev] = hc/lambda [m^2kgs^-2]
    coeff = h*c #12.39841987 # hc [kev*A]
    if energy > 0:
        wavelength = coeff/energy # hc/energy
    else:
        wavelength = 0 
    return wavelength
