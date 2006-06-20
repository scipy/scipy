import model
import formula
import regression
import robust
import family
from glm import Model as glm
from rlm import Model as rlm

import unittest
def suite():
    return unittest.TestSuite([tests.suite()])

