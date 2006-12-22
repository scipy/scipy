import unittest

from scipy.sandbox.models.tests import test_formula
from scipy.sandbox.models.tests import test_regression
from scipy.sandbox.models.tests import test_utils

def suite():
    return unittest.TestSuite([test_formula.suite(),
                               test_regression.suite(),
                               test_utils.suite()])
