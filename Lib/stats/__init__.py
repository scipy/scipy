""" Statistical Functions
"""
_modules = ['pstat']
_namespaces = ['rv','stats']

__all__ = []
import scipy
scipy.modules2all(__all__, _modules, globals())
scipy.names2all(__all__, _namespaces, globals())
del scipy

#---- testing ----#
def test(level=10):
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite(level=level))
    return runner

def test_suite(level=1):
    import scipy_test
    import scipy.stats
    this_mod = scipy.stats
    return scipy_test.harvest_test_suites(this_mod,level=level)
