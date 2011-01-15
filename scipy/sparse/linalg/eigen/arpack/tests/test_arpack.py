#!/usr/bin/env python
__usage__ = """
To run tests locally:
  python tests/test_arpack.py [-l<int>] [-v<int>]

"""
import sys, platform

import numpy as np

from numpy.testing import assert_almost_equal, assert_array_almost_equal, \
        assert_array_almost_equal_nulp, TestCase, run_module_suite, dec, \
        assert_raises, verbose, assert_equal

from numpy import array, finfo, argsort, dot, round, conj, random
from scipy.sparse import csc_matrix, isspmatrix
from scipy.sparse.linalg import LinearOperator
from scipy.sparse.linalg.eigen.arpack import eigs, eigsh, svds, \
     ArpackNoConvergence

from scipy.linalg import svd

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

        S1['v0']= array([0.39574246391875789,
                         0.00086496039750016962,
                         -0.9227205789982591,
                         -0.9165671495278005,
                         0.1175963848841306,
                         -0.29962625203712179])

        S1['eval']=array([0,1,2,2,5,6])
        self.symmetric.append(S1)

        N1={}
        N1['mat']=\
            array([[-2., -8.,  1.,  2., -5.],
                   [ 6.,  6.,  0.,  2.,  1.],
                   [ 0.,  4., -2., 11.,  0.],
                   [ 1.,  6.,  1.,  0., -4.],
                   [ 2., -6.,  4.,  9., -3]])

        N1['v0'] = array([0.39574246391875789,
                          0.00086496039750016962,
                          -0.9227205789982591,
                          -0.9165671495278005,
                          0.1175963848841306])

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
            l=k//2
            h=k//2+k%2
            low=range(len(eval))[:l]
            high=range(len(eval))[-h:]
            return eval[low+high]

    def eval_evec(self,d,typ,k,which,v0=None):
        a=d['mat'].astype(typ)
        if v0 == None:
            v0 = d['v0']
        exact_eval=self.get_exact_eval(d,typ,k,which)
        eval,evec=eigsh(a,k,which=which,v0=v0)
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


    def test_no_convergence(self):
        np.random.seed(1234)
        m = np.random.rand(30, 30)
        m = m + m.T
        try:
            w, v = eigsh(m, 4, which='LM', v0=m[:,0], maxiter=5)
            raise AssertionError("Spurious no-error exit")
        except ArpackNoConvergence, err:
            k = len(err.eigenvalues)
            if k <= 0:
                raise AssertionError("Spurious no-eigenvalues-found case")
            w, v = err.eigenvalues, err.eigenvectors
            for ww, vv in zip(w, v.T):
                assert_array_almost_equal(dot(m, vv), ww*vv,
                                          decimal=_ndigits['d'])


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

    def eval_evec(self,d,typ,k,which,v0=None):
        a=d['mat'].astype(typ)
        if v0 == None:
            v0 = d['v0']
        # get exact eigenvalues
        exact_eval=d['eval'].astype(typ)
        ind=self.sort_choose(exact_eval,typ,k,which)
        exact_eval=exact_eval[ind]
        # compute eigenvalues
        eval,evec=eigs(a,k,which=which,v0=v0)
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


    def test_no_convergence(self):
        np.random.seed(1234)
        m = np.random.rand(30, 30) + 1j*np.random.rand(30, 30)
        try:
            w, v = eigs(m, 3, which='LM', v0=m[:,0], maxiter=30)
            raise AssertionError("Spurious no-error exit")
        except ArpackNoConvergence, err:
            k = len(err.eigenvalues)
            if k <= 0:
                raise AssertionError("Spurious no-eigenvalues-found case")
            w, v = err.eigenvalues, err.eigenvectors
            for ww, vv in zip(w, v.T):
                assert_array_almost_equal(dot(m, vv), ww*vv,
                                          decimal=_ndigits['D'])

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


    def eval_evec(self,d,typ,k,which,v0=None):
        a=d['mat'].astype(typ)
        if v0 == None:
            v0 = d['v0']
        # get exact eigenvalues
        exact_eval=d['eval'].astype(typ.upper())
        ind=self.sort_choose(exact_eval,typ,k,which)
        exact_eval=exact_eval[ind]
        # compute eigenvalues
        eval,evec=eigs(a,k,which=which,v0=v0)
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

    def test_no_convergence(self):
        np.random.seed(1234)
        m = np.random.rand(30, 30)
        try:
            w, v = eigs(m, 3, which='LM', v0=m[:,0], maxiter=30)
            raise AssertionError("Spurious no-error exit")
        except ArpackNoConvergence, err:
            k = len(err.eigenvalues)
            if k <= 0:
                raise AssertionError("Spurious no-eigenvalues-found case")
            w, v = err.eigenvalues, err.eigenvectors
            for ww, vv in zip(w, v.T):
                assert_array_almost_equal(dot(m, vv), ww*vv,
                                          decimal=_ndigits['d'])

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

    def eval_evec(self,d,typ,k,which,v0=None):
        a=d['mat'].astype(typ)
        if v0 == None:
            v0 = d['v0']
        # get exact eigenvalues
        exact_eval=d['eval'].astype(typ.upper())
        ind=self.sort_choose(exact_eval,typ,k,which)
        exact_eval=exact_eval[ind]
        if verbose >= 3:
            print "exact"
            print exact_eval


        # compute eigenvalues
        eval,evec=eigs(a,k,which=which,v0=v0)
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


    def test_no_convergence(self):
        np.random.seed(1234)
        m = np.random.rand(30, 30) + 1j*np.random.rand(30, 30)
        try:
            w, v = eigs(m, 3, which='LM', v0=m[:,0], maxiter=30)
            raise AssertionError("Spurious no-error exit")
        except ArpackNoConvergence, err:
            k = len(err.eigenvalues)
            if k <= 0:
                raise AssertionError("Spurious no-eigenvalues-found case")
            w, v = err.eigenvalues, err.eigenvectors
            for ww, vv in zip(w, v.T):
                assert_array_almost_equal(dot(m, vv), ww*vv,
                                          decimal=_ndigits['D'])

def test_eigen_bad_shapes():
    # A is not square.
    A = csc_matrix(np.zeros((2,3)))
    assert_raises(ValueError, eigs, A)

def test_eigs_operator():
    # Check inferring LinearOperator dtype
    fft_op = LinearOperator((6, 6), np.fft.fft)
    w, v = eigs(fft_op, k=3)
    assert_equal(w.dtype, np.complex_)

def sorted_svd(m, k):
    """Compute svd of a dense matrix m, and return singular vectors/values
    sorted."""
    if isspmatrix(m):
        m = m.todense()
    u, s, vh = svd(m)
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
        y = np.array([[1, 2, 3, 8],
                      [3, 4, 3, 5],
                      [1, 0, 2, 3],
                      [0, 0, 1, 0]], np.float)
        z = csc_matrix(x)

        for m in [x.T, x, y, z, z.T]:
            for k in range(1, min(m.shape)):
                u, s, vh = sorted_svd(m, k)
                su, ss, svh = svds(m, k)

                m_hat = svd_estimate(u, s, vh)
                sm_hat = svd_estimate(su, ss, svh)

                assert_array_almost_equal_nulp(m_hat, sm_hat, nulp=1000)

    def test_simple_complex(self):
        x = np.array([[1, 2, 3],
                      [3, 4, 3],
                      [1+1j, 0, 2],
                      [0, 0, 1]], np.complex)
        y = np.array([[1, 2, 3, 8+5j],
                      [3-2j, 4, 3, 5],
                      [1, 0, 2, 3],
                      [0, 0, 1, 0]], np.complex)
        z = csc_matrix(x)

        for m in [x, x.T.conjugate(), x.T, y, y.conjugate(), z, z.T]:
            for k in range(1, min(m.shape)-1):
                u, s, vh = sorted_svd(m, k)
                su, ss, svh = svds(m, k)

                m_hat = svd_estimate(u, s, vh)
                sm_hat = svd_estimate(su, ss, svh)

                assert_array_almost_equal_nulp(m_hat, sm_hat, nulp=1000)

if __name__ == "__main__":
    run_module_suite()
