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

from scipy_test.testing import rand
def random(size):
    return rand(*size)

import sys
from scipy_test.testing import set_package_path
set_package_path()
from linalg import eig,eigvals,lu,svd,svdvals,cholesky,qr,schur,rsf2csf
del sys.path[0]

from scipy_test.testing import assert_array_almost_equal
from scipy_test.testing import assert_almost_equal
from scipy_test.testing import ScipyTestCase
import unittest


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

    def bench_random(self):
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

    #XXX: need more tests

class test_lu_solve(unittest.TestCase):        
    def check_lu(self):
        a = scipy.stats.random((10,10))
        b = scipy.stats.random(10)
        
        x1 = scipy.linalg.solve(a,b)
        
        lu_a = scipy.linalg.lu_factor(a)
        x2 = scipy.linalg.lu_solve(lu_a,b)
        
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

class test_cholesky(ScipyTestCase):

    def check_simple(self):
        a = [[8,2,3],[2,9,3],[3,3,6]]
        c = cholesky(a)
        assert_array_almost_equal(Numeric.dot(Numeric.transpose(c),c),a)
        c = Numeric.transpose(c)
        a = Numeric.dot(c,Numeric.transpose(c))
        assert_array_almost_equal(cholesky(a,lower=1),c)

    def check_simple_complex(self):
        c = [[3,3+4j,5],[0,2,2+7j],[0,0,7]]
        a = Numeric.dot(Numeric.transpose(Numeric.conjugate(c)),c)
        assert_array_almost_equal(cholesky(a),c)
        c = Numeric.transpose(c)
        a = Numeric.dot(c,Numeric.transpose(Numeric.conjugate(c)))
        assert_array_almost_equal(cholesky(a,lower=1),c)

    def check_random(self):
        n = 20
        for k in range(2):
            c = random([n,n])
            for i in range(n):
                c[i,i] = 20*(.1+c[i,i])
                for j in range(i): c[i,j] = 0
            a = Numeric.dot(Numeric.transpose(c),c)
            assert_array_almost_equal(cholesky(a),c)
            c = Numeric.transpose(c)
            a = Numeric.dot(c,Numeric.transpose(c))
            assert_array_almost_equal(cholesky(a,lower=1),c)

    def check_random_complex(self):
        n = 20
        for k in range(2):
            c = random([n,n])+1j*random([n,n])
            for i in range(n):
                c[i,i] = 20*(.1+abs(c[i,i]))
                for j in range(i): c[i,j] = 0
            a = Numeric.dot(Numeric.transpose(Numeric.conjugate(c)),c)
            assert_array_almost_equal(cholesky(a),c)
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

#####################################
def test_suite(level=1):
    suites = []
    if level > 0:
        suites.append( unittest.makeSuite(test_eigvals,'check_') )
        suites.append( unittest.makeSuite(test_eig,'check_') )
        suites.append( unittest.makeSuite(test_lu,'check_') )
        suites.append( unittest.makeSuite(test_svd,'check_') )
        suites.append( unittest.makeSuite(test_svdvals,'check_') )
        suites.append( unittest.makeSuite(test_cholesky,'check_') )
        suites.append( unittest.makeSuite(test_qr,'check_') )
    if level > 5:
        suites.append( unittest.makeSuite(test_eigvals,'bench_') )
        suites.append( unittest.makeSuite(test_eig,'bench_') )

    total_suite = unittest.TestSuite(suites)
    return total_suite

def test(level=10):
    all_tests = test_suite(level)
    runner = unittest.TextTestRunner()
    runner.run(all_tests)
    return runner

if __name__ == "__main__":
    if len(sys.argv)>1:
        level = eval(sys.argv[1])
    else:
        level = 1
    test(level)
