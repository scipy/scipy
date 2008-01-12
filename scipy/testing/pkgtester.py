''' Define test function for scipy package

Module tests for presence of nose.  If present returns NoseTester,
otherwise returns a placeholder test routine reporting lack of nose
and inability to run tests.  Typical use is in module __init__:

from scipy.testing.pkgtester import Tester
test = Tester().test

See nosetester module for test implementation

'''
try:
    import nose
except ImportError:
    from scipy.testing.nulltester import NullTester as Tester
else:
    from scipy.testing.nosetester import NoseTester as Tester

