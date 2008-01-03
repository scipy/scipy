#!/usr/bin/env python
__usage__ = """
First ensure that scipy core modules are installed.

Build interface to arpack
  python setup.py build
Run tests locally:
  python tests/test_arpack.py [-l<int>] [-v<int>]

"""

from scipy.testing import *

from scipy.sandbox.arpack import *

import numpy
from scipy.linalg import eig,eigh,norm

class TestEigenNonsymmetric(TestCase):

    def get_a1(self,typ):
        mat=numpy.array([[-2., -8.,  1.,  2., -5.],
                         [ 6.,  6.,  0.,  2.,  1.],
                         [ 0.,  4., -2., 11.,  0.],
                         [ 1.,  6.,  1.,  0., -4.],
                         [ 2., -6.,  4.,  9., -3]],typ)

        w=numpy.array([-2.21691+8.59661*1j,-2.21691-8.59661*1j,\
                       4.45961+3.80078*1j, 4.45961-3.80078*1j,\
                       -5.48541+0j],typ.upper())
        return mat,w


    def large_magnitude(self,typ,k):
        a,aw = self.get_a1(typ)
        w,v = eigen(a,k,which='LM')
        for i in range(k):
            assert_array_almost_equal(sb.dot(a,v[:,i]),w[i]*v[:,i],decimal=5)
        exact=numpy.abs(aw)
        num=numpy.abs(w)
        exact.sort()
        num.sort()
        assert_array_almost_equal(num[-k:],exact[-k:],decimal=5)

    def small_magnitude(self,typ,k):
        a,aw = self.get_a1(typ)
        w,v = eigen(a,k,which='SM')
        for i in range(k):
            assert_array_almost_equal(sb.dot(a,v[:,i]),w[i]*v[:,i],decimal=5)
        exact=numpy.abs(aw)
        num=numpy.abs(w)
        exact.sort()
        num.sort()
        assert_array_almost_equal(num[:k],exact[:k],decimal=5)


    def large_real(self,typ,k):
        a,aw = self.get_a1(typ)
        w,v = eigen(a,k,which='LR')
        for i in range(k):
            assert_array_almost_equal(sb.dot(a,v[:,i]),w[i]*v[:,i],decimal=5)
        exact=numpy.real(aw)
        num=numpy.real(w)
        exact.sort()
        num.sort()
        assert_array_almost_equal(num[-k:],exact[-k:],decimal=5)

    def small_real(self,typ,k):
        a,aw = self.get_a1(typ)
        w,v = eigen(a,k,which='SR')
        for i in range(k):
            assert_array_almost_equal(sb.dot(a,v[:,i]),w[i]*v[:,i],decimal=5)
        exact=numpy.real(aw)
        num=numpy.real(w)
        exact.sort()
        num.sort()
        assert_array_almost_equal(num[:k],exact[:k],decimal=5)


    def large_imag(self,typ,k):
        a,aw = self.get_a1(typ)
        w,v = eigen(a,k,which='LI')
        for i in range(k):
            assert_array_almost_equal(sb.dot(a,v[:,i]),w[i]*v[:,i],decimal=5)
        print w
        print aw
        exact=numpy.imag(aw)
        num=numpy.imag(w)
        exact.sort()
        num.sort()
        assert_array_almost_equal(num[-k:],exact[-k:],decimal=5)

    def small_imag(self,typ,k):
        a,aw = self.get_a1(typ)
        w,v = eigen(a,k,which='SI')
        for i in range(k):
            assert_array_almost_equal(sb.dot(a,v[:,i]),w[i]*v[:,i],decimal=5)
        exact=numpy.imag(aw)
        num=numpy.imag(w)
        exact.sort()
        num.sort()
        print num
        assert_array_almost_equal(num[:k],exact[:k],decimal=5)


    def test_type(self):
        k=2
        for typ in 'fd':
            self.large_magnitude(typ,k)
            self.small_magnitude(typ,k)
            self.large_real(typ,k)
            self.small_real(typ,k)
# Maybe my understanding of small imaginary and large imaginary
# isn't too keen.  I don't understand why these return
# different answers than in the complex case (the latter seems correct)
#            self.large_imag(typ,k)
#            self.small_imag(typ,k)



class TestEigenComplexNonsymmetric(TestCase):

    def get_a1(self,typ):
        mat=numpy.array([[-2., -8.,  1.,  2., -5.],
                         [ 6.,  6.,  0.,  2.,  1.],
                         [ 0.,  4., -2., 11.,  0.],
                         [ 1.,  6.,  1.,  0., -4.],
                         [ 2., -6.,  4.,  9., -3]],typ)

        w=numpy.array([-2.21691+8.59661*1j,-2.21691-8.59661*1j,\
                       4.45961+3.80078*1j, 4.45961-3.80078*1j,\
                       -5.48541+0j],typ.upper())
        return mat,w


    def large_magnitude(self,typ,k):
        a,aw = self.get_a1(typ)
        w,v = eigen(a,k,which='LM')
        for i in range(k):
            assert_array_almost_equal(sb.dot(a,v[:,i]),w[i]*v[:,i],decimal=5)
        exact=numpy.abs(aw)
        num=numpy.abs(w)
        exact.sort()
        num.sort()
        assert_array_almost_equal(num,exact[-k:],decimal=5)

    def small_magnitude(self,typ,k):
        a,aw = self.get_a1(typ)
        w,v = eigen(a,k,which='SM')
        for i in range(k):
            assert_array_almost_equal(sb.dot(a,v[:,i]),w[i]*v[:,i],decimal=5)
        exact=numpy.abs(aw)
        num=numpy.abs(w)
        exact.sort()
        num.sort()
        assert_array_almost_equal(num,exact[:k],decimal=5)


    def large_real(self,typ,k):
        a,aw = self.get_a1(typ)
        w,v = eigen(a,k,which='LR')
        for i in range(k):
            assert_array_almost_equal(sb.dot(a,v[:,i]),w[i]*v[:,i],decimal=5)
        exact=numpy.real(aw)
        num=numpy.real(w)
        exact.sort()
        num.sort()
        assert_array_almost_equal(num,exact[-k:],decimal=5)

    def small_real(self,typ,k):
        a,aw = self.get_a1(typ)
        w,v = eigen(a,k,which='SR')
        for i in range(k):
            assert_array_almost_equal(sb.dot(a,v[:,i]),w[i]*v[:,i],decimal=5)
        exact=numpy.real(aw)
        num=numpy.real(w)
        exact.sort()
        num.sort()
        assert_array_almost_equal(num,exact[:k],decimal=5)


    def large_imag(self,typ,k):
        a,aw = self.get_a1(typ)
        w,v = eigen(a,k,which='LI')
        for i in range(k):
            assert_array_almost_equal(sb.dot(a,v[:,i]),w[i]*v[:,i],decimal=5)
        exact=numpy.imag(aw)
        num=numpy.imag(w)
        exact.sort()
        num.sort()
        assert_array_almost_equal(num,exact[-k:],decimal=5)

    def small_imag(self,typ,k):
        a,aw = self.get_a1(typ)
        w,v = eigen(a,k,which='SI')
        for i in range(k):
            assert_array_almost_equal(sb.dot(a,v[:,i]),w[i]*v[:,i],decimal=5)
        exact=numpy.imag(aw)
        num=numpy.imag(w)
        exact.sort()
        num.sort()
        assert_array_almost_equal(num,exact[:k],decimal=5)


    def test_type(self):
        k=2
        for typ in 'FD':
            self.large_magnitude(typ,k)
            self.small_magnitude(typ,k)
            self.large_real(typ,k)
            self.small_real(typ,k)
            self.large_imag(typ,k)
            self.small_imag(typ,k)




class TestEigenSymmetric(TestCase):

    def get_a1(self,typ):
        mat_a1=numpy.array([[ 2.,  0.,  0., -1.,  0., -1.],
                            [ 0.,  2.,  0., -1.,  0., -1.],
                            [ 0.,  0.,  2., -1.,  0., -1.],
                            [-1., -1., -1.,  4.,  0., -1.],
                            [ 0.,  0.,  0.,  0.,  1., -1.],
                            [-1., -1., -1., -1., -1.,  5.]],
                           typ)
        w = [0,1,2,2,5,6] # eigenvalues of a1
        return mat_a1,w

    def large_eigenvalues(self,typ,k):
        a,aw = self.get_a1(typ)
        w,v = eigen_symmetric(a,k,which='LM',tol=1e-7)
        assert_array_almost_equal(w,aw[-k:])

    def small_eigenvalues(self,typ,k):
        a,aw = self.get_a1(typ)
        w,v = eigen_symmetric(a,k,which='SM')
        assert_array_almost_equal(w,aw[:k])

    def end_eigenvalues(self,typ,k):
        a,aw = self.get_a1(typ)
        w,v = eigen_symmetric(a,k,which='BE')
        exact=[aw[0],aw[-1]]
        assert_array_almost_equal(w,exact)

    def large_eigenvectors(self,typ,k):
        a,aw = self.get_a1(typ)
        w,v = eigen_symmetric(a,k,which='LM')
        ew,ev = eigh(a)
        ind=ew.argsort()
        assert_array_almost_equal(w,numpy.take(ew,ind[-k:]))
        for i in range(k):
            assert_array_almost_equal(sb.dot(a,v[:,i]),w[i]*v[:,i])

    def small_eigenvectors(self,typ,k):
        a,aw = self.get_a1(typ)
        w,v = eigen_symmetric(a,k,which='SM',tol=1e-7)
        ew,ev = eigh(a)
        ind=ew.argsort()
        assert_array_almost_equal(w,numpy.take(ew,ind[:k]))
        for i in range(k):
            assert_array_almost_equal(sb.dot(a,v[:,i]),w[i]*v[:,i])

    def end_eigenvectors(self,typ,k):
        a,aw = self.get_a1(typ)
        w,v = eigen_symmetric(a,k,which='BE')
        ew,ev = eigh(a)
        ind=ew.argsort()
        exact=numpy.concatenate(([ind[:k/2],ind[-k/2:]]))
        assert_array_almost_equal(w,numpy.take(ew,exact))
        for i in range(k):
            assert_array_almost_equal(sb.dot(a,v[:,i]),w[i]*v[:,i])

    def test_eigenvectors(self):
        k=2
        for typ in 'fd':
            self.large_eigenvectors(typ,k)
            self.small_eigenvectors(typ,k)
            self.end_eigenvectors(typ,k)

    def test_type(self):
        k=2
        for typ in 'fd':
            self.large_eigenvalues(typ,k)
            self.small_eigenvalues(typ,k)
            self.end_eigenvalues(typ,k)


class TestEigenComplexSymmetric(TestCase):

    def get_a1(self,typ):
        mat_a1=numpy.array([[ 2.,  0.,  0., -1.,  0., -1.],
                            [ 0.,  2.,  0., -1.,  0., -1.],
                            [ 0.,  0.,  2., -1.,  0., -1.],
                            [-1., -1., -1.,  4.,  0., -1.],
                            [ 0.,  0.,  0.,  0.,  1., -1.],
                            [-1., -1., -1., -1., -1.,  5.]],
                           typ)
        w = numpy.array([0+0j,1+0j,2+0j,2+0j,5+0j,6+0j]) # eigenvalues of a1
        return mat_a1,w

    def large_magnitude(self,typ,k):
        a,aw = self.get_a1(typ)
        w,v = eigen(a,k,which='LM')
        for i in range(k):
            assert_array_almost_equal(sb.dot(a,v[:,i]),w[i]*v[:,i],decimal=5)
        aw.real.sort()
        w.real.sort()
        assert_array_almost_equal(w,aw[-k:])


    def small_magnitude(self,typ,k):
        a,aw = self.get_a1(typ)
        w,v = eigen(a,k,which='SM')
        for i in range(k):
            assert_array_almost_equal(sb.dot(a,v[:,i]),w[i]*v[:,i])
        aw.real.sort()
        w.real.sort()
        assert_array_almost_equal(w,aw[:k])

    def large_real(self,typ,k):
        a,aw = self.get_a1(typ)
        w,v = eigen(a,k,which='LR')
        for i in range(k):
            assert_array_almost_equal(sb.dot(a,v[:,i]),w[i]*v[:,i],decimal=5)
        aw.real.sort()
        w.real.sort()
        assert_array_almost_equal(w,aw[-k:],decimal=5)


    def small_real(self,typ,k):
        a,aw = self.get_a1(typ)
        w,v = eigen(a,k,which='SR')
        for i in range(k):
            assert_array_almost_equal(sb.dot(a,v[:,i]),w[i]*v[:,i])
        aw.real.sort()
        w.real.sort()
        assert_array_almost_equal(w,aw[:k])

    def test_complex_symmetric(self):
        k=2
        for typ in 'FD':
            self.large_magnitude(typ,k)
            self.small_magnitude(typ,k)
            self.large_real(typ,k)
            self.small_real(typ,k)



if __name__ == "__main__":
    unittest.main()
