"Iterative Solvers for Sparse Linear Systems"

#from info import __doc__
from iterative import *
from minres import minres
from lgmres import lgmres
from lsqr import lsqr

__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import Tester
test = Tester().test
bench = Tester().bench
