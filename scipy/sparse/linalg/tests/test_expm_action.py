
import numpy as np
from numpy.testing import (TestCase, run_module_suite, assert_allclose)

import scipy.sparse.linalg._expm_action

class TestTheta(TestCase):

    def test_theta(self):
        for m in range(1, 55+1):
            theta = scipy.sparse.linalg._expm_action._theta(m)
            print m, theta
        raise Exception


if __name__ == '__main__':
    run_module_suite()

