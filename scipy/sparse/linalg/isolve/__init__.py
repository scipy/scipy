"Iterative Solvers for Sparse Linear Systems"

from __future__ import division, print_function, absolute_import

#from info import __doc__
from .iterative import *
from .minres import minres
from .lgmres import lgmres
from .lsqr import lsqr
from .lsmr import lsmr

__all__ = [s for s in dir() if not s.startswith('_')]
from numpy.testing import Tester
test = Tester().test
