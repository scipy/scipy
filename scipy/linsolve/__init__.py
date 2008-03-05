
from warnings import warn

warn('scipy.linsolve has moved to scipy.sparse.linalg.dsolve', DeprecationWarning)

from scipy.sparse.linalg.dsolve import *
