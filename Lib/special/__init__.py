#
# special - Special Functions
#

from info import __doc__
#from special_version import special_version as __version__

from basic import *
import specfun
import orthogonal
from orthogonal import legendre, chebyt, chebyu, chebyc, chebys, \
     jacobi, laguerre, genlaguerre, hermite, hermitenorm, gegenbauer, \
     sh_legendre, sh_chebyt, sh_chebyu, sh_jacobi, poch

__all__ = filter(lambda s:not s.startswith('_'),dir())
from numpy.testing import ScipyTest 
test = ScipyTest().test
