import model
import formula
import regression
import robust
import family
from glm import model as glm
from rlm import model as rlm


import unittest
def suite():
    return unittest.TestSuite([tests.suite()])


from numpy.testing import ScipyTest
test = ScipyTest().test
