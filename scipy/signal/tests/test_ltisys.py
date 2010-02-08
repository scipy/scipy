import numpy as np
from numpy.testing import *

from scipy.signal.ltisys import ss2tf

class TestSS2TF:
    def tst_matrix_shapes(self, p, q, r):
        ss2tf(np.zeros((p, p)),
              np.zeros((p, q)),
              np.zeros((r, p)),
              np.zeros((r, q)), 0)

    def test_basic(self):
        for p, q, r in [
            (3, 3, 3),
            (0, 3, 3),
            (1, 1, 1)]:
            yield self.tst_matrix_shapes, p, q, r

if __name__ == "__main__":
    run_module_suite()
