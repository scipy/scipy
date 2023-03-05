import numpy as np
from scipy.linalg import _qr_solve

class TestQRSolver(unittest.TestCase):
    def test_is_no_soultion(self):
        m,n = (100,150)
        A = np.random.rand( m, n )
        y = np.random.rand( m )
        x_0 = _qr_solve.qr_solve(A, y, silent=True)
        return np.allclose(y, A@x_0)
    
    def test_is_solution(self):
        m,n = (100,150)
        A = np.random.rand( m, n )
        x = np.random.rand( n )
        y = A@x
        x_0 = _qr_solve.qr_solve(A, y, silent=True)
        return np.allclose(y, A@x_0)