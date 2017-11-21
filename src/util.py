import numpy as np
from astropy import constants

PI = np.pi
C = constants.c.value
G = constants.G.value
H = constants.h.value
K = constants.k_B.value
SIGMA = constants.sigma_sb.value
R_SUN = constants.R_sun.value
M_SUN = constants.M_sun.value
AU = constants.au.value
PC = constants.pc.value
# 1 Jy = 1e-26 W/(m^2*Hz)
JY = 1e-26

# planck's law for specific intensity at specified frequencies        
def planck(freq, temp):
    # suppress numpy errors for overflow in exponentiation
    np.seterr(all='ignore')
    coeff = 2*H*np.power(freq,3) / C**2
    exp = 1 / (np.exp(np.divide(H*freq,(K*temp))) - 1)
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

# binary search
def bin_search(func, l, r, value):
    mid = (l + r) / 2
    val = func(mid)
    dif = abs(val - value)
    avg = abs(val + value) / 2
    # if percent difference is less than 1/1000000, accept as ~ equal
    if (dif/avg) < 1/1000000:
        return mid
    elif val > value:
        return bin_search(func,l,mid,value)
    else:
        return bin_search(func,mid,r,value)