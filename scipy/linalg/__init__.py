#
# linalg - Dense Linear Algebra routines
#

from info import __doc__
from linalg_version import linalg_version as __version__

from basic import *
from decomp import *
from matfuncs import *
from blas import *

from iterative import *


__all__ = filter(lambda s:not s.startswith('_'),dir())

from numpy.dual import register_func
for k in ['norm', 'inv', 'svd', 'solve', 'det', 'eig', 'eigh', 'eigvals',
          'eigvalsh', 'lstsq', 'cholesky']:
    try:
        register_func(k, eval(k))
    except ValueError:
        pass

try:
    register_func('pinv', pinv2)
except ValueError:
    pass

del k, register_func

from numpy.testing import Tester
test = Tester().test
bench = Tester().bench
