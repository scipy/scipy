#
# special - Special Functions
#

from pre___init__ import __doc__
from special_version import special_version as __version__

from basic import *
import specfun
import orthogonal
from orthogonal import legendre, chebyt, chebyu, chebyc, chebys, \
     jacobi, laguerre, genlaguerre, hermite, hermitenorm, gegenbauer, \
     sh_legendre, sh_chebyt, sh_chebyu, sh_jacobi


################## test functions #########################

def test(level=10):
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite(level))
    return runner

def test_suite(level=1):
    import scipy_test.testing
    exec 'import %s as this_mod' % (__name__)
    return scipy_test.testing.harvest_test_suites(this_mod,level=level)
