import numpy as np
from astropy import constants

PI = np.pi
C = constants.c.value
H = constants.h.value
K = constants.k_B.value
SIGMA = constants.sigma_sb.value
R_SUN = constants.R_sun.value
AU = constants.au.value

# planck's law for specific intensity at specified frequencies        
def planck(freq, temp):
    coeff = 2*H*np.power(freq,3)/C**2
    exp = 1/(np.exp(np.divide(H*freq,(K*temp))) - 1)
    return np.multiply(coeff,exp)

# the blackbody emission spectrum, integrated over solid angle
def spectrum(temp):
    def _function(freq):
        return PI * planck(freq,temp)
    return _function

# converts floats to integers if they're integers
def ftoi(value):
    return int(value) if value.is_integer() else value

# returns 'num' with 'n' significant figures
def sig_figs(num, n=1):
    numstr = ("{0:.%ie}" % (n-1)).format(num)
    return float(numstr)