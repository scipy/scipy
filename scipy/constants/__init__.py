"""
Various useful constants and conversion formulae

Modules
-------

.. autosummary::
   :toctree: generated/

   codata - CODATA Recommended Values of Fundamental Physical Const (2006)
   constants - Collection of physical constants and conversion factors

Functions
---------

.. autosummary::
   :toctree: generated/

   C2F - Convert Celsius to Fahrenheit
   C2K - Convert Celsius to Kelvin
   F2C - Convert Fahrenheit to Celsius
   F2K - Convert Fahrenheit to Kelvin
   K2C - Convert Kelvin to Celsius
   K2F - Convert Kelvin to Fahrenheit
   find - Find the codata.physical_constant keys containing a given string
   lambda2nu - Convert wavelength to optical frequency
   nu2lambda - Convert optical frequency to wavelength
   precision - Relative precision in physical_constants indexed by key
   unit - Unit in physical_constants indexed by key
   value - Value in physical_constants indexed by key

"""

# Modules contributed by BasSw (wegwerp@gmail.com)
from codata import *
from constants import *

__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import Tester
test = Tester().test
