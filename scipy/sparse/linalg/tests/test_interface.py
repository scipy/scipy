"""Test functions for the sparse.linalg._interface module
"""

import numpy as np
import scipy as sp

from scipy.sparse.linalg import aslinearoperator, LinearOperator


def define_matricies(d):
    # several ways to construct diagonal matrix as linear operator
    
    A = np.diag(d)
    A_linop = []
    
    # Numpy
    A_linop.append( aslinearoperator(np.diag(d)) )
    
    # Scipy sparse
    A_linop.append( aslinearoperator(sp.sparse.diags_array(d)) )
    
    # Abstract functions
    matmul = lambda x: (d * x.T).T
    A_linop.append( LinearOperator((5,5), matmul, tmatmul=matmul, dtype=d.dtype) )
    
    # Explicit conversion via AsLinearOperatorDunderProtocol
    class ExplicitConversion:
        def __init__(self, d):
            self.d = d
        def __aslinearoperator__(self):
            matmul = lambda x: (self.d * x.T).T
            shape = (len(self.d), len(self.d))
            return LinearOperator(shape, matmul, tmatmul=matmul, dtype=self.d.dtype)
    
    A_linop.append( aslinearoperator(ExplicitConversion(d)) )
    
    # Implicit conversion via MatrixProtocol
    class ImplicitConversion:
        ndim = 2
        def __init__(self, d):
            self.d = d
            self.shape = (len(d), len(d))
            self.dtype = d.dtype
            
        def __matmul__(self, other):
            return (self.d * other.T).T
        
        def transpose(self):
            return self
    
    A_linop.append( aslinearoperator(ImplicitConversion(d)) )
    
    return A, A_linop

def _test_matricies(A, A_linop, x):
    Ax = A @ x
    ATx = A.transpose() @ x
    AHx = A.transpose().conjugate() @ x
    ACx = A.conjugate() @ x
    
    for Ai in A_linop:
        assert np.allclose(Ai @ x, Ax)
        assert np.allclose(Ai.T @ x, ATx)
        assert np.allclose(Ai.H @ x, AHx)
        assert np.allclose(Ai.C @ x, ACx)

def test1():
    d = np.arange(5)
    A, A_linop = define_matricies(d)

    _test_matricies(A, A_linop, np.random.rand(5))
    _test_matricies(A, A_linop, np.random.rand(5,4))
    _test_matricies(A, A_linop, np.random.rand(5) + 1j*np.random.rand(5))
    _test_matricies(A, A_linop, np.random.rand(5,4) + 1j*np.random.rand(5,4))

    d = np.arange(5) + 1j*np.arange(5)
    A, A_linop = define_matricies(d)

    _test_matricies(A, A_linop, np.random.rand(5))
    _test_matricies(A, A_linop, np.random.rand(5,4))
    _test_matricies(A, A_linop, np.random.rand(5) + 1j*np.random.rand(5))
    _test_matricies(A, A_linop, np.random.rand(5,4) + 1j*np.random.rand(5,4))
    
if __name__ == "__main__":
    test1()