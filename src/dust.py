import numpy as np
from scipy.integrate import simps, quad
from decimal import Decimal
from util import *
from star import Star

# assume that dust has density ~2 g/cc = 2000 kg/m^3
_DENSITY = 2000

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
        self._mass = _DENSITY * (4/3*PI*self._radius**3)
        self._star = star
        self._a = float(Decimal(a))

    def mass(self):
        return self._mass

    def radius(self):
        return float(Decimal(self._radius))

    def semimajor_axis(self):
        return float(Decimal(self._a))

    def freq(self):
        if self._wavelength is None:
            return None
        return np.divide(C, self._wavelength)

    def star(self):
        return self._star

    # the total power incident from the star
    # 'perfect' should be True if the dust is a perfect absorber
    def power(self, perfect=False):
        if self._power is not None:
            return self._power
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
            flux_in = self._star.spectrum(self._a)
            power = PI*self._radius**2 * \
                simps(np.multiply(flux_in(freq),self._q),freq) * JY
            self._power = power
        # the dust is a perfect absorber
        else:
            # suppress numpy errors for overflow in exponentiation
            np.seterr(all='ignore')
            flux_in = self._star.spectrum(self._a)
            power = PI*self._radius**2 * \
                quad(flux_in,0,3*10**18)[0] * JY
            self._power = power
        return self._power

    # the equilibrium temperature of the dust
    def temp(self, perfect=False):
        if self._temp is not None:
            return self._temp
        power = self.power(perfect)
        # the dust is a perfect absorber, treat it as a blackbody
        if self._wavelength is None:
            return (power/(4*PI*self._radius**2*SIGMA))**(1/4)
        # the dust is not a perfect absorber, account for Q
        def _p_out(temp):
            area = 4*PI*self._radius**2
            freq = np.divide(C, self._wavelength)
            return area*PI*simps(np.multiply(planck(freq,temp),self._q),freq)
        # numerically solve for temperature using binary search
        self._temp = bin_search(_p_out,0,self._star.temp(),power)
        return self._temp

    # returns the flux density (Jy) from dust at an orbital distance 'd' 
    # as a function of frequency (aka the emission spectrum)
    def spectrum(self, d):
        temp = self.temp()
        # the dust is a perfect absorber, treat as a perfect blackbody
        if self._wavelength is None:
            reference = spectrum(temp)
            def _function(freq):
                return reference(freq) * \
                    float((Decimal(self._radius)/Decimal(d))**Decimal(2)) / JY
            return _function
        # the dust is not a perfect absorber
        def _function(freq):
            return PI * np.multiply(planck(freq, temp), self._q) * \
                float((Decimal(self._radius)/Decimal(d))**Decimal(2)) / JY
        return _function

    # returns the force exerted on the dust from Poynting-Robertson drag
    def pr_drag(self):
        v = np.sqrt(G*self._star.mass()/self._a)
        if self._wavelength is not None:
            return (v / C**2) * self.power()
        return (v / C**2) * self.power(True)

    # returns the force exerted on the dust from radiation pressure
    def radiation_pressure(self):
        if self._wavelength is not None:
            return self.power() / C
        return self.power(True) / C

    # returns the time needed to remove the dust from the system
    def removal_time(self):
        # beta is the ratio of radiation pressure's and gravity's forces
        rp = self.radiation_pressure()
        fg = G * self._star.mass() * self._mass / self._a**2
        beta = rp / fg
        # if beta > 1, radiation pressure blows it out of the system
        if beta > 1:
            return 0
        # if beta <= 1, poynting-robertson may drag it into the star
        return 1/4 * C/(G*self._star.mass()*beta) \
            * (self._a**2 - self._star.radius()**2)