from decimal import Decimal
from util import *

class Star:
    _name = None
    _radius = 0
    _temp = 0

    def __init__(self, name, radius, temp):
        self._name = name
        self._radius = float(Decimal(radius))
        self._temp = float(Decimal(temp))

    def name(self):
        return self._name

    def radius(self):
        return float(Decimal(self._radius))

    def temp(self):
        return float(Decimal(self._temp))

    # returns the flux density (Jy) from a star at an orbital distance 'a' 
    # as a function of frequency
    def spectrum(self, a):
        reference = spectrum(self._temp)
        def _function(freq):
            return reference(freq) * (self._radius/a)**2
        return _function