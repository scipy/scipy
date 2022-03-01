"Iterative Solvers for Sparse Linear Systems"

#from info import __doc__
from .iterative import *
from .minres import minres
from .lgmres import lgmres
from .lsqr import lsqr
from .lsmr import lsmr
from ._gcrotmk import gcrotmk
from .tfqmr import tfqmr
from .gmresk import gmresk

__all__ = [
    'bicg', 'bicgstab', 'cg', 'cgs', 'gcrotmk', 'gmres',
    'lgmres', 'lsmr', 'lsqr',
    'minres', 'qmr', 'tfqmr', 'gmresk'
]

from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
