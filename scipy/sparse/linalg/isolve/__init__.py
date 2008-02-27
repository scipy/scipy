"Iterative Solvers for Sparse Linear Systems"

#from info import __doc__
from iterative import *
from minres import minres

__all__ = filter(lambda s:not s.startswith('_'),dir())
from scipy.testing.pkgtester import Tester
test = Tester().test
bench = Tester().bench

