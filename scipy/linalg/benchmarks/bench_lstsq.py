# -*- coding: utf-8 -*-
""" Benchmark linalg.lstsq for various LAPACK drivers

"""

from __future__ import division, print_function, absolute_import

import time
import math

import numpy as np
from numpy.testing import assert_allclose
import scipy.linalg as la
import numpy.random as rnd

def bench_lstsq_1():
    """
    Test the speed of four least-squares solvers on not full rank matrices.
    Also check the difference in the solutions.

    The matrix has the size (m,2/3*m) the rank is 1/2*m.
    Matrix values are random in the range (-5, 5), the same is for the right 
    hand side. The complex matrix is the sum of real and imaginary matrices.
    """
    np.random.seed(1234)
    print()
    print('                       Least-Squares Solver')
    print('==============================================================')
    print('      shape      |   LAPACK   |        dtype       |   time   ')
    print('                 |   driver   |                    | (seconds)')
    print('--------------------------------------------------------------')
    fmt = ' %15s |   %6s  | %18s | %6.2f '
    m_steps = range(500, 3500, 500)
    for (i,m) in enumerate(m_steps):
        for dtype in (np.float64, np.complex128):

            n = math.ceil(2./3. * m)
            k = math.ceil(1./2. * m)

            if (dtype == np.complex128):
                A = (10 * rnd.rand(m,k) - 5) + 1j*(10 * rnd.rand(m,k) - 5)
                Temp = (10 * rnd.rand(k,n) - 5) + 1j*(10 * rnd.rand(k,n) - 5)
                b = (10 * rnd.rand(m,1) - 5) + 1j*(10 * rnd.rand(m,1) - 5) 
            else:
                A = (10 * rnd.rand(m,k) - 5)
                Temp = 10 * rnd.rand(k,n) - 5
                b = 10 * rnd.rand(m,1) - 5
                
            A = A.dot(Temp)
            
            A1 = A.copy() 
            b1 = b.copy()
            t1 = time.time()
            res1 = la.lstsq(A1, b1, cond=None, overwrite_a=False, 
                            overwrite_b=False, check_finite=False, 
                            lapack_driver='gelss')
            t1 = time.time() - t1
            x1 = res1[0]
            print(fmt % (A.shape, '*gelss', A.dtype, t1))
    
            A2 = A.copy() 
            b2 = b.copy()
            t2 = time.time()
            res2 = la.lstsq(A2, b2, cond=None, overwrite_a=False, 
                         overwrite_b=False, check_finite=False, 
                         lapack_driver='gelsy')
            t2 = time.time() - t2
            x2 = res2[0]
            print(fmt % (A.shape, '*gelsy', A.dtype, t2))
    
            A3 = A.copy()
            b3 = b.copy()
            t3 = time.time()
            res3 = la.lstsq(A3, b3, cond=None, overwrite_a=False, 
                         overwrite_b=False, check_finite=False, 
                         lapack_driver='gelsd')
            t3 = time.time() - t3
            x3 = res3[0]
            print(fmt % (A.shape, '*gelsd', A.dtype, t3))
            
            A4 = A.copy()
            b4 = b.copy()
            t4 = time.time()
            res4 = np.linalg.lstsq(A4,b4,rcond=np.finfo(A4.dtype).eps*100)
            t4 = time.time() - t4
            x4 = res4[0]
            print(fmt % (A.shape, ' Numpy', A.dtype, t4))

            # Check that the results are the same for all drivers
            assert_allclose(x4,x1,atol=1000*np.finfo(A1.dtype).eps,
                            rtol=1000*np.finfo(A1.dtype).eps)
            assert_allclose(x4,x2,atol=1000*np.finfo(A2.dtype).eps,
                            rtol=1000*np.finfo(A2.dtype).eps)
            assert_allclose(x4,x3,atol=1000*np.finfo(A3.dtype).eps,
                            rtol=1000*np.finfo(A3.dtype).eps)
                            
def bench_lstsq_2():
    """
    Test the scaling of least-squares solvers only with respect to one dimesion, 
    Another dimension is fixed. The rank is k.

    The matrix has size (m,n) the rank is k.
    Matrix values are random in the range (-5, 5), the same is for right hand
    side. The complex matrix is the sum of real and imaginary matrices.
    """

    np.random.seed(1234)
    print()
    print('                       Least-Squares Solver')
    print('==============================================================')
    print('      shape      |   LAPACK   |        dtype       |   time   ')
    print('                 |   driver   |                    | (seconds)')
    print('--------------------------------------------------------------')
    fmt = ' %15s |   %6s  | %18s | %6.2f '
    m_steps = range(10000, 120000, 10000)
    for (i,m) in enumerate(m_steps):
        for dtype in (np.float64, np.complex128):
  
            n = 100
            k = 80

            if (dtype == np.complex128):
                A = (10 * rnd.rand(m,k) - 5) + 1j*(10 * rnd.rand(m,k) - 5)
                Temp = (10 * rnd.rand(k,n) - 5) + 1j*(10 * rnd.rand(k,n) - 5)
                b = (10 * rnd.rand(m,1) - 5) + 1j*(10 * rnd.rand(m,1) - 5) 
            else:
                A = (10 * rnd.rand(m,k) - 5)
                Temp = 10 * rnd.rand(k,n) - 5
                b = 10 * rnd.rand(m,1) - 5
                
            A = A.dot(Temp)

            A1 = A.copy()
            b1 = b.copy()
            t1 = time.time()
            res1 = la.lstsq(A1, b1, cond=None, overwrite_a=False, 
                            overwrite_b=False, check_finite=False, 
                            lapack_driver='gelss')
            t1 = time.time() - t1
            x1 = res1[0]
            print(fmt % (A.shape, '*gelss', A.dtype, t1))
    
            A2 = A.copy()
            b2 = b.copy()
            t2 = time.time()
            res2 = la.lstsq(A2, b2, cond=None, overwrite_a=False, 
                         overwrite_b=False, check_finite=False, 
                         lapack_driver='gelsy')
            t2 = time.time() - t2
            x2 = res2[0]
            print(fmt % (A.shape, '*gelsy', A.dtype, t2))
    
            A3 = A.copy()
            b3 = b.copy()
            t3 = time.time()
            res3 = la.lstsq(A3, b3, cond=None, overwrite_a=False, 
                         overwrite_b=False, check_finite=False, 
                         lapack_driver='gelsd')
            t3 = time.time() - t3
            x3 = res3[0]
            print(fmt % (A.shape, '*gelsd', A.dtype, t3))
            
            A4 = A.copy()
            b4 = b.copy()
            t4 = time.time()
            res4 = np.linalg.lstsq(A4,b4,rcond=np.finfo(A4.dtype).eps * 100)
            t4 = time.time() - t4
            x4 = res4[0]
            print(fmt % (A.shape, ' Numpy', A.dtype, t4))

            # Check that the results are the same for all drivers
            assert_allclose(x4,x1,atol=1000*np.finfo(A1.dtype).eps,
                            rtol=1000*np.finfo(A1.dtype).eps)
            assert_allclose(x4,x2,atol=1000*np.finfo(A2.dtype).eps,
                            rtol=1000*np.finfo(A2.dtype).eps)
            assert_allclose(x4,x3,atol=1000*np.finfo(A3.dtype).eps,
                            rtol=1000*np.finfo(A3.dtype).eps)

if __name__ == '__main__':
    bench_lstsq_1()
    bench_lstsq_2()

