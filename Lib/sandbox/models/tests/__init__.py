import test_formula
import test_regression
import test_utils
import unittest

def suite():
    return unittest.TestSuite([test_formula.suite(),
                               test_regression.suite(),
                               test_utils.suite()])
