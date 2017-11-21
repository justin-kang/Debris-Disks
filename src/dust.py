import numpy as np
from scipy.integrate import simps, quad
from decimal import Decimal
from util import *
from star import Star

class Dust:
    _mass = 0
    _radius = 0
    _a = 0
    # the q factors and corresponding wavelengths for non-perfectly absorbers
    _wavelength = None
    _q = None
    # the star this dust is orbiting around, power/temp of dust due to star
    _star = None
    _power = None
    _temp = None

    def __init__(self, radius, star, a):
        # set this to precise numbers, because we need to get file names
        self._radius = float(Decimal(radius) * Decimal(10)**Decimal(-6))
        # assume that dust has density ~2 g/cc
        self._mass = 2000 * (4/3*PI*self._radius**3)
        self._star = star
        self._a = float(Decimal(a))

    def mass(self):
        return self._mass

    def radius(self):
        return float(Decimal(self._radius))

    def axis(self):
        return float(Decimal(self._a))

    def freq(self):
        if self._wavelength is None:
            return None
        return np.divide(C, self._wavelength)

    def star(self):
        return self._star

    # the total power in on the dust
    # 'perfect' should be True if the dust is a perfect absorber
    def power(self, a=None, perfect=False, flag=0):
        if self._power is not None and (a is None or a == self._a):
            return self._power
        if a is None:
            a = self._a
        # if the dust is not a perfect absorber, adjust the spectrum to 
        # match the frequency range provided by the file
        if not perfect:
            # read in the Q as a function of wavelength from a file
            if self._wavelength is None:
                w = []
                q = []
                radius = ftoi(self._radius * 10**6)
                with open(str(radius)+'micron.txt') as input_file:
                    for line in input_file:
                        line = line.strip().split()
                        w.append(float(Decimal(line[0]))*1e-6)
                        q.append(float(Decimal(line[1])))
                self._wavelength = w
                self._q = q
            # calculate the power in
            freq = np.divide(C, self._wavelength)
            flux_in = self._star.spectrum(a)
            power = PI*self._radius**2 * \
                simps(np.multiply(flux_in(freq),self._q),freq) * JY
            if flag:
                return power
            self._power = power
        # the dust is a perfect absorber
        else:
            # suppress numpy errors for overflow in exponentiation
            np.seterr(all='ignore')
            flux_in = self._star.spectrum(a)
            power = PI*self._radius**2 * \
                quad(flux_in,0,3*10**18)[0] * JY
            if flag:
                return power
            self._power = power
        return self._power

    # the equilibrium temperature of the dust
    def temp(self, a=None, perfect=False):
        if self._temp is not None and a is None:
            return self._temp
        power = self.power(a,perfect)
        # the dust is a perfect absorber
        if self._wavelength is None:
            return (power/(4*PI*self._radius**2*SIGMA))**(1/4)
        # the dust is not a perfect absorber
        def _p_out(temp):
            area = 4*PI*self._radius**2
            freq = np.divide(C, self._wavelength)
            return area*PI*simps(np.multiply(planck(freq,temp),self._q),freq)
        # numerically solve for temperature using binary search
        self._temp = bin_search(_p_out,0,self._star.temp(),power)
        return self._temp

    # the emission spectrum of the dust
    def spectrum(self):
        temp = self.temp()
        # the dust is a perfect absorber
        if self._wavelength is None:
            reference = spectrum(temp)
            def _function(freq):
                return reference(freq) / JY
            return _function
        # the dust is not a perfect absorber
        def _function(freq):
            return 4*PI * np.multiply(planck(freq, temp), self._q) / JY
        return _function

    # returns the force exerted on the dust from Poynting-Robertson drag
    def pr_drag(self, a=None):
        flag = 1
        if a is None or a == self._a:
            a = self._a
            flag = 0
        v = np.sqrt(G*self._star.mass()/a)
        if self._wavelength is not None:
            return (v / C**2) * self.power(a,flag=flag)
        return (v / C**2) * self.power(a,True,flag)

    # returns the force exerted on the dust from radiation pressure
    def radiation_pressure(self, a=None):
        flag = 1
        if a is None or a == self._a:
            a = self._a
            flag = 0
        if self._wavelength is not None:
            return self.power(a,flag=flag) / C
        return self.power(a,True,flag) / C