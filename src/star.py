from decimal import Decimal
from util import *

class Star:
    _name = None
    _mass = 0
    _radius = 0
    _temp = 0

    def __init__(self, name, mass, radius, temp):
        self._name = name
        self._mass = mass
        self._radius = float(Decimal(radius))
        self._temp = float(Decimal(temp))

    def name(self):
        return self._name

    def mass(self):
        return self._mass

    def radius(self):
        return float(Decimal(self._radius))

    def temp(self):
        return float(Decimal(self._temp))

    # returns the flux density (Jy) from a star at an orbital distance 'a' 
    # as a function of frequency
    def spectrum(self, a):
        reference = spectrum(self._temp)
        def _function(freq):
            return reference(freq) * float((Decimal(self._radius)/Decimal(a))
                **Decimal(2)) / JY
        return _function