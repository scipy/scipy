#!/usr/bin/env python
#
# Created by: Pearu Peterson, March 2002
#
""" Test functions for linalg.decomp module

"""
__usage__ = """
Build linalg:
  python setup_linalg.py build
Run tests if scipy is installed:
  python -c 'import scipy;scipy.linalg.test(<level>)'
Run tests if linalg is not installed:
  python tests/test_decomp.py [<level>]
"""

import Numeric
dot = Numeric.dot

import sys
from scipy_test.testing import *
set_package_path()
from linalg import eig,eigvals,lu,svd,svdvals,cholesky,qr,schur,rsf2csf
from linalg import lu_solve,lu_factor,solve,diagsvd,hessenberg
import scipy_base
del sys.path[0]

import unittest
def random(size):
    return rand(*size)

class test_eigvals(ScipyTestCase):

    def check_simple(self):
        a = [[1,2,3],[1,2,3],[2,5,6]]
        w = eigvals(a)
        exact_w = [(9+Numeric.sqrt(93))/2,0,(9-Numeric.sqrt(93))/2]
        assert_array_almost_equal(w,exact_w)

    def check_simple_tr(self):
        a = Numeric.array([[1,2,3],[1,2,3],[2,5,6]],'d')
        a = Numeric.transpose(a).copy()
        a = Numeric.transpose(a)
        w = eigvals(a)
        exact_w = [(9+Numeric.sqrt(93))/2,0,(9-Numeric.sqrt(93))/2]
        assert_array_almost_equal(w,exact_w)

    def check_simple_complex(self):
        a = [[1,2,3],[1,2,3],[2,5,6+1j]]
        w = eigvals(a)
        exact_w = [(9+1j+Numeric.sqrt(92+6j))/2,
                   0,
                   (9+1j-Numeric.sqrt(92+6j))/2]
        assert_array_almost_equal(w,exact_w)

    def bench_random(self,level=5):
        import LinearAlgebra
        Numeric_eigvals = LinearAlgebra.eigenvalues
        print
        print '           Finding matrix eigenvalues'
        print '      =================================='
        print '      |    contiguous     '#'|   non-contiguous '
        print '----------------------------------------------'
        print ' size |  scipy  '#'| Numeric |  scipy  | Numeric'
        
        for size,repeat in [(20,150),(100,7),(200,2)]:
            repeat *= 1
            print '%5s' % size,
            sys.stdout.flush()
            
            a = random([size,size])

            print '| %6.2f ' % self.measure('eigvals(a)',repeat),
            sys.stdout.flush()
            
            print '   (secs for %s calls)' % (repeat)

class test_eig(ScipyTestCase):

    def check_simple(self):
        a = [[1,2,3],[1,2,3],[2,5,6]]
        w,v = eig(a)
        exact_w = [(9+Numeric.sqrt(93))/2,0,(9-Numeric.sqrt(93))/2]
        v0 = Numeric.array([1,1,(1+Numeric.sqrt(93)/3)/2])
        v1 = Numeric.array([3.,0,-1])
        v2 = Numeric.array([1,1,(1-Numeric.sqrt(93)/3)/2])
        v0 = v0 / Numeric.sqrt(Numeric.dot(v0,Numeric.transpose(v0)))
        v1 = v1 / Numeric.sqrt(Numeric.dot(v1,Numeric.transpose(v1)))
        v2 = v2 / Numeric.sqrt(Numeric.dot(v2,Numeric.transpose(v2)))
        assert_array_almost_equal(w,exact_w)
        assert_array_almost_equal(v0,v[:,0]*Numeric.sign(v[0,0]))
        assert_array_almost_equal(v1,v[:,1]*Numeric.sign(v[0,1]))
        assert_array_almost_equal(v2,v[:,2]*Numeric.sign(v[0,2]))
        for i in range(3):
            assert_array_almost_equal(Numeric.dot(a,v[:,i]),w[i]*v[:,i])
        w,v = eig(a,left=1,right=0)
        for i in range(3):
            assert_array_almost_equal(Numeric.matrixmultiply(\
                Numeric.transpose(a),v[:,i]),w[i]*v[:,i])

    def check_simple_complex(self):
        a = [[1,2,3],[1,2,3],[2,5,6+1j]]
        w,vl,vr = eig(a,left=1,right=1)
        for i in range(3):
            assert_array_almost_equal(Numeric.dot(a,vr[:,i]),w[i]*vr[:,i])
        for i in range(3):
            assert_array_almost_equal(Numeric.matrixmultiply(\
            Numeric.conjugate(Numeric.transpose(a)),vl[:,i]),
                                      Numeric.conjugate(w[i])*vl[:,i])

class test_lu(ScipyTestCase):

    def check_simple(self):
        a = [[1,2,3],[1,2,3],[2,5,6]]
        p,l,u = lu(a)
        assert_array_almost_equal(Numeric.dot(Numeric.dot(p,l),u),a)
        pl,u = lu(a,permute_l=1)
        assert_array_almost_equal(Numeric.dot(pl,u),a)

    def check_simple_complex(self):
        a = [[1,2,3],[1,2,3],[2,5j,6]]
        p,l,u = lu(a)
        assert_array_almost_equal(Numeric.dot(Numeric.dot(p,l),u),a)
        pl,u = lu(a,permute_l=1)
        assert_array_almost_equal(Numeric.dot(pl,u),a)

    #XXX: need more tests

class test_lu_solve(unittest.TestCase):        
    def check_lu(self):
        a = random((10,10))
        b = random((10,))
        
        x1 = solve(a,b)
        
        lu_a = lu_factor(a)
        x2 = lu_solve(lu_a,b)
        
        assert_array_equal(x1,x2)

class test_svd(ScipyTestCase):

    def check_simple(self):
        a = [[1,2,3],[1,20,3],[2,5,6]]
        u,s,v = svd(a)
        sigma = Numeric.zeros((u.shape[0],v.shape[0]),s.typecode())
        for i in range(len(s)): sigma[i,i] = s[i]
        assert_array_almost_equal(dot(dot(u,sigma),v),a)

    def check_simple_singular(self):
        a = [[1,2,3],[1,2,3],[2,5,6]]
        u,s,v = svd(a)
        sigma = Numeric.zeros((u.shape[0],v.shape[0]),s.typecode())
        for i in range(len(s)): sigma[i,i] = s[i]
        assert_array_almost_equal(dot(dot(u,sigma),v),a)

    def check_simple_underdet(self):
        a = [[1,2,3],[4,5,6]]
        u,s,v = svd(a)
        sigma = Numeric.zeros((u.shape[0],v.shape[0]),s.typecode())
        for i in range(len(s)): sigma[i,i] = s[i]
        assert_array_almost_equal(dot(dot(u,sigma),v),a)

    def check_simple_overdet(self):
        a = [[1,2],[4,5],[3,4]]
        u,s,v = svd(a)
        sigma = Numeric.zeros((u.shape[0],v.shape[0]),s.typecode())
        for i in range(len(s)): sigma[i,i] = s[i]
        assert_array_almost_equal(dot(dot(u,sigma),v),a)

    def check_random(self):
        n = 20
        m = 15
        for i in range(3):
            for a in [random([n,m]),random([m,n])]:
                u,s,v = svd(a)
                sigma = Numeric.zeros((u.shape[0],v.shape[0]),s.typecode())
                for i in range(len(s)): sigma[i,i] = s[i]
                assert_array_almost_equal(dot(dot(u,sigma),v),a)

    def check_simple_complex(self):
        a = [[1,2,3],[1,2j,3],[2,5,6]]
        u,s,v = svd(a)
        sigma = Numeric.zeros((u.shape[0],v.shape[0]),s.typecode())
        for i in range(len(s)): sigma[i,i] = s[i]
        assert_array_almost_equal(dot(dot(u,sigma),v),a)

    def check_random_complex(self):
        n = 20
        m = 15
        for i in range(3):
            for a in [random([n,m]),random([m,n])]:
                a = a + 1j*random(list(a.shape))
                u,s,v = svd(a)
                sigma = Numeric.zeros((u.shape[0],v.shape[0]),s.typecode())
                for i in range(len(s)): sigma[i,i] = s[i]
                assert_array_almost_equal(dot(dot(u,sigma),v),a)

class test_svdvals(ScipyTestCase):

    def check_simple(self):
        a = [[1,2,3],[1,2,3],[2,5,6]]
        s = svdvals(a)
        assert len(s)==3
        assert s[0]>=s[1]>=s[2]

    def check_simple_underdet(self):
        a = [[1,2,3],[4,5,6]]
        s = svdvals(a)
        assert len(s)==2
        assert s[0]>=s[1]

    def check_simple_overdet(self):
        a = [[1,2],[4,5],[3,4]]
        s = svdvals(a)
        assert len(s)==2
        assert s[0]>=s[1]

    def check_simple_complex(self):
        a = [[1,2,3],[1,20,3j],[2,5,6]]
        s = svdvals(a)
        assert len(s)==3
        assert s[0]>=s[1]>=s[2]

    def check_simple_underdet_complex(self):
        a = [[1,2,3],[4,5j,6]]
        s = svdvals(a)
        assert len(s)==2
        assert s[0]>=s[1]

    def check_simple_overdet_complex(self):
        a = [[1,2],[4,5],[3j,4]]
        s = svdvals(a)
        assert len(s)==2
        assert s[0]>=s[1]

class test_diagsvd(ScipyTestCase):

    def check_simple(self):
        assert_array_almost_equal(diagsvd([1,0,0],3,3),[[1,0,0],[0,0,0],[0,0,0]])

class test_cholesky(ScipyTestCase):

    def check_simple(self):
        a = [[8,2,3],[2,9,3],[3,3,6]]
        c = cholesky(a)
        assert_array_almost_equal(Numeric.dot(Numeric.transpose(c),c),a)
        c = Numeric.transpose(c)
        a = Numeric.dot(c,Numeric.transpose(c))
        assert_array_almost_equal(cholesky(a,lower=1),c)

    def check_simple_complex(self):
        m = Numeric.array([[3+1j,3+4j,5],[0,2+2j,2+7j],[0,0,7+4j]])
        a = Numeric.dot(Numeric.transpose(Numeric.conjugate(m)),m)
        c = cholesky(a)
        a1 = Numeric.dot(Numeric.transpose(Numeric.conjugate(c)),c)
        assert_array_almost_equal(a,a1)
        c = Numeric.transpose(c)
        a = Numeric.dot(c,Numeric.transpose(Numeric.conjugate(c)))
        assert_array_almost_equal(cholesky(a,lower=1),c)

    def check_random(self):
        n = 20
        for k in range(2):
            m = random([n,n])
            for i in range(n):
                m[i,i] = 20*(.1+m[i,i])
            a = Numeric.dot(Numeric.transpose(m),m)
            c = cholesky(a)
            a1 = Numeric.dot(Numeric.transpose(c),c)
            assert_array_almost_equal(a,a1)
            c = Numeric.transpose(c)
            a = Numeric.dot(c,Numeric.transpose(c))
            assert_array_almost_equal(cholesky(a,lower=1),c)

    def check_random_complex(self):
        n = 20
        for k in range(2):
            m = random([n,n])+1j*random([n,n])
            for i in range(n):
                m[i,i] = 20*(.1+abs(m[i,i]))
            a = Numeric.dot(Numeric.transpose(Numeric.conjugate(m)),m)
            c = cholesky(a)
            a1 = Numeric.dot(Numeric.transpose(Numeric.conjugate(c)),c)
            assert_array_almost_equal(a,a1)
            c = Numeric.transpose(c)
            a = Numeric.dot(c,Numeric.transpose(Numeric.conjugate(c)))
            assert_array_almost_equal(cholesky(a,lower=1),c)


class test_qr(ScipyTestCase):

    def check_simple(self):
        a = [[8,2,3],[2,9,3],[5,3,6]]
        q,r = qr(a)
        assert_array_almost_equal(Numeric.dot(q,r),a)

    def check_simple_complex(self):
        a = [[3,3+4j,5],[5,2,2+7j],[3,2,7]]
        q,r = qr(a)
        assert_array_almost_equal(Numeric.dot(q,r),a)

    def check_random(self):
        n = 20
        for k in range(2):
            a = random([n,n])
            q,r = qr(a)
            assert_array_almost_equal(Numeric.dot(q,r),a)

    def check_random_complex(self):
        n = 20
        for k in range(2):
            a = random([n,n])+1j*random([n,n])
            q,r = qr(a)
            assert_array_almost_equal(Numeric.dot(q,r),a)

dot = Numeric.dot
transp = Numeric.transpose
conj = Numeric.conjugate
any = Numeric.sometrue
ravel = Numeric.ravel
iscomplex = scipy_base.iscomplex

class test_schur(ScipyTestCase):

    def check_simple(self):
        a = [[8,12,3],[2,9,3],[10,3,6]]
        t,z = schur(a)
        assert_array_almost_equal(dot(dot(z,t),transp(conj(z))),a)
        tc,zc = schur(a,'complex')
        assert(any(ravel(iscomplex(zc))) and any(ravel(iscomplex(tc))))
        assert_array_almost_equal(dot(dot(zc,tc),transp(conj(zc))),a)
        tc2,zc2 = rsf2csf(tc,zc)
        assert_array_almost_equal(dot(dot(zc2,tc2),transp(conj(zc2))),a)

class test_hessenberg(ScipyTestCase):

    def check_simple(self):
        a = [[-149, -50,-154],
             [ 537, 180, 546],
             [ -27,  -9, -25]]
        h1 = [[-149.0000,42.2037,-156.3165],
              [-537.6783,152.5511,-554.9272],
              [0,0.0728, 2.4489]]
        h,q = hessenberg(a,calc_q=1)
        assert_array_almost_equal(dot(transp(q),dot(a,q)),h)
        assert_array_almost_equal(h,h1,decimal=4)

    def check_simple_complex(self):
        a = [[-149, -50,-154],
             [ 537, 180j, 546],
             [ -27j,  -9, -25]]
        h,q = hessenberg(a,calc_q=1)
        h1 = dot(transp(conj(q)),dot(a,q))
        assert_array_almost_equal(h1,h)

    def check_simple2(self):
        a = [[1,2,3,4,5,6,7],
             [0,2,3,4,6,7,2],
             [0,2,2,3,0,3,2],
             [0,0,2,8,0,0,2],
             [0,3,1,2,0,1,2],
             [0,1,2,3,0,1,0],
             [0,0,0,0,0,1,2]]
        h,q = hessenberg(a,calc_q=1)
        assert_array_almost_equal(dot(transp(q),dot(a,q)),h)

    def check_random(self):
        n = 20
        for k in range(2):
            a = random([n,n])
            h,q = hessenberg(a,calc_q=1)
            assert_array_almost_equal(dot(transp(q),dot(a,q)),h)

    def check_random_complex(self):
        n = 20
        for k in range(2):
            a = random([n,n])+1j*random([n,n])
            h,q = hessenberg(a,calc_q=1)
            h1 = dot(transp(conj(q)),dot(a,q))
            assert_array_almost_equal(h1,h)

if __name__ == "__main__":
    ScipyTest('linalg.decomp').run()
