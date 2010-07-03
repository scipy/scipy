#!/usr/bin/env python
__usage__ = """
To run tests locally:
  python tests/test_arpack.py [-l<int>] [-v<int>]

"""
import sys, platform

import numpy as np

from numpy.testing import *

from numpy import array, finfo, argsort, dot, round, conj, random
from scipy.sparse.linalg.eigen.arpack import eigen_symmetric, eigen, svd

from scipy.linalg import svd as dsvd

def assert_almost_equal_cc(actual,desired,decimal=7,err_msg='',verbose=True):
    # almost equal or complex conjugates almost equal
    try:
        assert_almost_equal(actual,desired,decimal,err_msg,verbose)
    except:
        assert_almost_equal(actual,conj(desired),decimal,err_msg,verbose)


def assert_array_almost_equal_cc(actual,desired,decimal=7,
                                 err_msg='',verbose=True):
    # almost equal or complex conjugates almost equal
    try:
        assert_array_almost_equal(actual,desired,decimal,err_msg,verbose)
    except:
        assert_array_almost_equal(actual,conj(desired),decimal,err_msg,verbose)


# check if we're on 64-bit OS X, there these tests fail.
if sys.platform == 'darwin' and platform.architecture()[0] == '64bit':
    _osx64bit = True
else:
    _osx64bit = False

# precision for tests
_ndigits = {'f':4, 'd':12, 'F':4, 'D':12}

class TestArpack(TestCase):

    def setUp(self):
        self.symmetric=[]
        self.nonsymmetric=[]

        S1={}
        S1['mat']=\
        array([[ 2.,  0.,  0., -1.,  0., -1.],
               [ 0.,  2.,  0., -1.,  0., -1.],
               [ 0.,  0.,  2., -1.,  0., -1.],
               [-1., -1., -1.,  4.,  0., -1.],
               [ 0.,  0.,  0.,  0.,  1., -1.],
               [-1., -1., -1., -1., -1.,  5.]])

        S1['eval']=array([0,1,2,2,5,6])
        self.symmetric.append(S1)

        N1={}
        N1['mat']=\
            array([[-2., -8.,  1.,  2., -5.],
                   [ 6.,  6.,  0.,  2.,  1.],
                   [ 0.,  4., -2., 11.,  0.],
                   [ 1.,  6.,  1.,  0., -4.],
                   [ 2., -6.,  4.,  9., -3]])

        N1['eval']=\
            array([ -5.4854094033782888+0.0j,
                     -2.2169058544873783+8.5966096591588261j,
                     -2.2169058544873783-8.5966096591588261j,
                     4.4596105561765107+3.8007839204319454j,
                     4.4596105561765107-3.8007839204319454j],'D')



        self.nonsymmetric.append(N1)


class TestEigenSymmetric(TestArpack):

    def get_exact_eval(self,d,typ,k,which):
        eval=d['eval'].astype(typ)
        ind=argsort(eval)
        eval=eval[ind]
        if which=='LM':
            return eval[-k:]
        if which=='SM':
            return eval[:k]
        if which=='BE':
            # one ev from each end - if k is odd, extra ev on high end
            l=k/2
            h=k/2+k%2
            low=range(len(eval))[:l]
            high=range(len(eval))[-h:]
            return eval[low+high]

    def eval_evec(self,d,typ,k,which,**kwds):
        a=d['mat'].astype(typ)
        exact_eval=self.get_exact_eval(d,typ,k,which)
        eval,evec=eigen_symmetric(a,k,which=which,**kwds)
        # check eigenvalues
        assert_array_almost_equal(eval,exact_eval,decimal=_ndigits[typ])
        # check eigenvectors A*evec=eval*evec
        for i in range(k):
            assert_array_almost_equal(dot(a,evec[:,i]),
                                      eval[i]*evec[:,i],
                                      decimal=_ndigits[typ])

    @dec.knownfailureif(_osx64bit, "Currently fails on 64-bit OS X 10.6")
    def test_symmetric_modes(self):
        k=2
        for typ in 'fd':
            for which in ['LM','SM','BE']:
                self.eval_evec(self.symmetric[0],typ,k,which)

    @dec.knownfailureif(_osx64bit, "Currently fails on 64-bit OS X 10.6")
    def test_starting_vector(self):
        k=2
        for typ in 'fd':
            A=self.symmetric[0]['mat']
            n=A.shape[0]
            v0 = random.rand(n).astype(typ)
            self.eval_evec(self.symmetric[0],typ,k,which='LM',v0=v0)


class TestEigenComplexSymmetric(TestArpack):

    def sort_choose(self,eval,typ,k,which):
        # sort and choose the eigenvalues and eigenvectors
        # both for the exact answer and that returned from ARPACK
        reval=round(eval,decimals=_ndigits[typ])
        ind=argsort(reval)
        if which=='LM' or which=='LR':
            return ind[-k:]
        if which=='SM' or which=='SR':
            return ind[:k]

    def eval_evec(self,d,typ,k,which):
        a=d['mat'].astype(typ)
        # get exact eigenvalues
        exact_eval=d['eval'].astype(typ)
        ind=self.sort_choose(exact_eval,typ,k,which)
        exact_eval=exact_eval[ind]
        # compute eigenvalues
        eval,evec=eigen(a,k,which=which)
        ind=self.sort_choose(eval,typ,k,which)
        eval=eval[ind]
        evec=evec[:,ind]

        # check eigenvalues
        assert_array_almost_equal(eval,exact_eval,decimal=_ndigits[typ])
        # check eigenvectors A*evec=eval*evec
        for i in range(k):
            assert_array_almost_equal(dot(a,evec[:,i]),
                                      eval[i]*evec[:,i],
                                      decimal=_ndigits[typ])

    @dec.knownfailureif(_osx64bit, "Currently fails on 64-bit OS X 10.6")
    def test_complex_symmetric_modes(self):
        k=2
        for typ in 'FD':
            for which in ['LM','SM','LR','SR']:
                self.eval_evec(self.symmetric[0],typ,k,which)



class TestEigenNonSymmetric(TestArpack):


    def sort_choose(self,eval,typ,k,which):
        reval=round(eval,decimals=_ndigits[typ])
        if which in ['LR','SR']:
            ind=argsort(reval.real)
        elif which in ['LI','SI']:
            # for LI,SI ARPACK returns largest,smallest abs(imaginary) why?
            ind=argsort(abs(reval.imag))
        else:
            ind=argsort(abs(reval))

        if which in ['LR','LM','LI']:
            return ind[-k:]
        if which in ['SR','SM','SI']:
            return ind[:k]


    def eval_evec(self,d,typ,k,which,**kwds):
        a=d['mat'].astype(typ)
        # get exact eigenvalues
        exact_eval=d['eval'].astype(typ.upper())
        ind=self.sort_choose(exact_eval,typ,k,which)
        exact_eval=exact_eval[ind]
        # compute eigenvalues
        eval,evec=eigen(a,k,which=which,**kwds)
        ind=self.sort_choose(eval,typ,k,which)
        eval=eval[ind]
        evec=evec[:,ind]
        # check eigenvalues
        # check eigenvectors A*evec=eval*evec
        for i in range(k):
            assert_almost_equal_cc(eval[i],exact_eval[i],decimal=_ndigits[typ])
            assert_array_almost_equal_cc(dot(a,evec[:,i]),
                                      eval[i]*evec[:,i],
                                      decimal=_ndigits[typ])


    @dec.knownfailureif(_osx64bit, "Currently fails on 64-bit OS X 10.6")
    def test_nonsymmetric_modes(self):
        k=2
        for typ in 'fd':
            for which in ['LI','LR','LM','SM','SR','SI']:
                for m in self.nonsymmetric:
                    self.eval_evec(m,typ,k,which)



    @dec.knownfailureif(_osx64bit, "Currently fails on 64-bit OS X 10.6")
    def test_starting_vector(self):
        k=2
        for typ in 'fd':
            A=self.symmetric[0]['mat']
            n=A.shape[0]
            v0 = random.rand(n).astype(typ)
            self.eval_evec(self.symmetric[0],typ,k,which='LM',v0=v0)




class TestEigenComplexNonSymmetric(TestArpack):

    def sort_choose(self,eval,typ,k,which):
        eps=finfo(typ).eps
        reval=round(eval,decimals=_ndigits[typ])
        if which in ['LR','SR']:
            ind=argsort(reval)
        elif which in ['LI','SI']:
            ind=argsort(reval.imag)
        else:
            ind=argsort(abs(reval))

        if which in ['LR','LI','LM']:
            return ind[-k:]
        if which in ['SR','SI','SM']:
            return ind[:k]

    def eval_evec(self,d,typ,k,which):
        a=d['mat'].astype(typ)
        # get exact eigenvalues
        exact_eval=d['eval'].astype(typ.upper())
        ind=self.sort_choose(exact_eval,typ,k,which)
        exact_eval=exact_eval[ind]
        if verbose >= 3:
            print "exact"
            print exact_eval


        # compute eigenvalues
        eval,evec=eigen(a,k,which=which)
        ind=self.sort_choose(eval,typ,k,which)
        eval=eval[ind]
        evec=evec[:,ind]
        if verbose >= 3:
            print eval
        # check eigenvalues
        # check eigenvectors A*evec=eval*evec
        for i in range(k):
            assert_almost_equal_cc(eval[i],exact_eval[i],decimal=_ndigits[typ])
            assert_array_almost_equal_cc(dot(a,evec[:,i]),
                                      eval[i]*evec[:,i],
                                      decimal=_ndigits[typ])

    @dec.knownfailureif(_osx64bit, "Currently fails on 64-bit OS X 10.6")
    def test_complex_nonsymmetric_modes(self):
        k=2
        for typ in 'FD':
            for which in ['LI','LR','LM','SI','SR','SM']:
                for m in self.nonsymmetric:
                    self.eval_evec(m,typ,k,which)

def sorted_svd(m, k):
    """Compute svd of a dense matrix m, and return singular vectors/values
    sorted."""
    u, s, vh = dsvd(m)
    ii = np.argsort(s)[-k:]

    return u[:, ii], s[ii], vh[ii]

def svd_estimate(u, s, vh):
    return np.dot(u, np.dot(np.diag(s), vh))

class TestSparseSvd(TestCase):
    def test_simple_real(self):
        x = np.array([[1, 2, 3],
                      [3, 4, 3],
                      [1, 0, 2],
                      [0, 0, 1]], np.float)

        for m in [x.T, x]:
            for k in range(1, 3):
                u, s, vh = sorted_svd(m, k)
                su, ss, svh = svd(m, k)

                m_hat = svd_estimate(u, s, vh)
                sm_hat = svd_estimate(su, ss, svh)

                assert_array_almost_equal_nulp(m_hat, sm_hat, nulp=1000)

    @dec.knownfailureif(True, "Complex sparse SVD not implemented (depends on "
                              "Hermitian support in eigen_symmetric")
    def test_simple_complex(self):
        x = np.array([[1, 2, 3],
                      [3, 4, 3],
                      [1+1j, 0, 2],
                      [0, 0, 1]], np.complex)

        for m in [x, x.T.conjugate()]:
            for k in range(1, 3):
                u, s, vh = sorted_svd(m, k)
                su, ss, svh = svd(m, k)

                m_hat = svd_estimate(u, s, vh)
                sm_hat = svd_estimate(su, ss, svh)

                assert_array_almost_equal_nulp(m_hat, sm_hat, nulp=1000)

if __name__ == "__main__":
    run_module_suite()
