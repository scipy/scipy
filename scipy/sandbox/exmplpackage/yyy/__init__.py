"""
yyy - Subpackage of Scipy module exmplpackage
"""

__all__ = ['fun', 'test']

# Import testing rig, allowing scipy.examplpackage.yyy.test()
from numpy.testing import Tester
test = Tester().test

def fun():
    return 'Hello from yyy.fun'
