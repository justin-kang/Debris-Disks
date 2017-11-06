import numpy as np
from scipy.integrate import simps, quad
from decimal import Decimal
from util import *
from star import Star

class Dust:
    _radius = 0
    _a = 0
    # the q factors and corresponding wavelengths for non-perfectly absorbers
    _wavelength = None
    _q = None

    # 'radius' is in microns
    def __init__(self, radius, a):
        # set this to precise numbers, because we need to get file names
        self._radius = float(Decimal(radius) * Decimal(10)**Decimal(-6))
        self._a = float(Decimal(a))

    def radius(self):
        return float(Decimal(self._radius))

    def semimajor_axis(self):
        return float(Decimal(self._a))

    # the total power in on the dust
    # 'perfect' should be not None if the dust is not a perfect absorber
    def power(self, star, perfect=True):
        # if the dust is not a perfect absorber, adjust the spectrum to match 
        # the frequency range provided by the file
        if not perfect:
            # read in the Q as a function of wavelength from a file
            w = []
            q = []
            radius = ftoi(self._radius * 10**6)
            with open(str(radius)+'micron.txt') as input_file:
                for line in input_file:
                    line = line.strip().split()
                    w.append(float(line[0]))
                    q.append(float(line[1]))
            self._wavelength = w
            self._q = q
            # calculate the power in
            freq = np.divide(C, self._wavelength)
            power_in = star.spectrum(self._a)
            return PI*self._radius**2 * \
                simps(np.multiply(power_in(freq),self._q),freq)
        # the dust is a perfect absorber
        # suppress numpy errors for overflow in exponentiation
        np.seterr(all='ignore')
        power_in = star.spectrum(self._a);
        return PI*self._radius**2 * quad(power_in,0,3*10**18)[0]

    # the equilibrium temperature of the dust
    # 'star' should be not None if the dust is not a perfect absorber
    def temp(self, power, star):
        # the dust is a perfect absorber
        if self._q is None:
            return (power/(4*PI*self._radius**2*SIGMA))**(1/4)
        # the dust is not a perfect absorber
        def _p_out(temp):
            area = 4*PI*self._radius**2
            freq = np.divide(C, self._wavelength)
            return area*PI*simps(np.multiply(self._q,planck(freq,temp)),freq)
        # numerically solve for temperature using binary search
        def _bin_search(func, l, r, power):
            mid = l + (r-l)/2
            value = func(mid)
            dif = abs(value - power)
            avg = abs(value + power) / 2
            #if value == power:
            if (dif/avg) < 0.0001:
                return mid
            elif value > power:
                return _bin_search(func,l,mid,power)
            else:
                return _bin_search(func,mid,r,power)
        return _bin_search(_p_out,0,star.temp(),power)

    # the emission spectrum of the dust
    # 'star' should be not None if the dust is not a perfect absorber
    def spectrum(self, power, star):
        temp = self.temp(power,star)
        # the dust is a perfect absorber
        if self._q is None:
            reference = spectrum(temp)
            def _function(freq):
                return reference(freq)
            return _function
        # the dust is not a perfect absorber
        freq = np.divide(C,self._wavelength)
        def _function(unused):
            return PI * simps(np.multiply(self._q,planck(freq,temp)),freq)
        return _function