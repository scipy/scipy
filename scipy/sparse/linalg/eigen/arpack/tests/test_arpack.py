__usage__ = """
To run tests locally:
  python tests/test_arpack.py [-l<int>] [-v<int>]

"""

import numpy as np

from numpy.testing import assert_almost_equal, assert_array_almost_equal, \
        assert_array_almost_equal_nulp, TestCase, run_module_suite, dec, \
        assert_raises, verbose, assert_equal

from numpy import array, finfo, argsort, dot, round, conj, random
from scipy.linalg import eig, eigh
from scipy.sparse import csc_matrix, csr_matrix, lil_matrix, isspmatrix
from scipy.sparse.linalg import LinearOperator, aslinearoperator
from scipy.sparse.linalg.eigen.arpack import eigs, eigsh, svds, \
     ArpackNoConvergence

from scipy.linalg import svd

def _aslinearoperator_with_dtype(m):
    m = aslinearoperator(m)
    if not hasattr(m, 'dtype'):
        x = np.zeros(m.shape[1])
        m.dtype = (m*x).dtype
    return m

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

# precision for tests
_ndigits = {'f':4, 'd':12, 'F':4, 'D':12}

# types of matrices to use
_mattypes = [csr_matrix,
             _aslinearoperator_with_dtype, 
             lambda x:x]

class TestArpack(TestCase):
    def argsort_which(self,eval,typ,k,which,sigma=None):
        if sigma is None:
            reval=round(eval,decimals=_ndigits[typ])
        else:
            reval=round(1./(eval-sigma),decimals=_ndigits[typ])
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

    def eval_evec(self,d,typ,k,which,sigma=None,v0=None,conv=None):
        general = ('bmat' in d)
        a=d['mat'].astype(typ)
        if conv is not None:
            a=conv(a)
        if v0 == None:
            v0 = d['v0']
        #get exact eigenvalues
        exact_eval=d['eval'].astype(typ.upper())
        ind=self.argsort_which(exact_eval,typ,k,which,sigma)
        exact_eval=exact_eval[ind]
        # compute eigenvalues
        if general:
            b=d['bmat'].astype(typ)
            if conv is not None:
                b = conv(b)
            eval,evec=self.eigs(a,k,b,which=which,v0=v0,sigma=sigma)
        else:
            eval,evec=self.eigs(a,k,which=which,v0=v0,sigma=sigma)
        ind=self.argsort_which(eval,typ,k,which,sigma)
        eval=eval[ind]

class TestSymmetric(TestArpack):
    def setUp(self):
        self.eigs = eigsh
        self.modes = ['LM','SM','BE']

        # standard symmetric problem
        SS = {}
        SS['mat']=array([[ 2.,  0.,  0., -1.,  0., -1.],
                         [ 0.,  2.,  0., -1.,  0., -1.],
                         [ 0.,  0.,  2., -1.,  0., -1.],
                         [-1., -1., -1.,  4.,  0., -1.],
                         [ 0.,  0.,  0.,  0.,  1., -1.],
                         [-1., -1., -1., -1., -1.,  5.]])
        SS['eval']=array([0,1,2,2,5,6])        
        SS['v0'] = array([0.39574246391875789,
                          0.00086496039750016962,
                          -0.9227205789982591,
                          -0.9165671495278005,
                          0.1175963848841306,
                          -0.29962625203712179])
        self.standard = [SS]
        
        # general symmetric problem
        SG = {}
        SG['mat'] = SS['mat']
        SG['bmat'] = array([[ 3., -3.,  0., -1.,  0., -1.],
                            [-3., 10.,  0.,  7.,  3.,  5.],
                            [ 0.,  0., 10.,  1.,  3.,  0.],
                            [-1.,  7.,  1., 11.,  0.,  1.],
                            [ 0.,  3.,  3.,  0.,  8., -1.],
                            [-1.,  5.,  0.,  1., -1., 10.]])
        SG['eval'] = array([0.0,
                            1.4360859298239867e-01,
                            2.0496461767247975e-01,
                            5.0507316521393908e-01,   
                            6.0317191682460070e-01,   
                            1.9075146642093827e+01])
        SG['v0']= SS['v0']
        self.general = [SG]
    
    def test_standard_symmetric_modes(self):
        k=2
        for d in self.standard:
            for typ in 'fd':
                for which in self.modes:
                    for sigma in (None,0.5):
                        for conv in _mattypes:
                            self.eval_evec(d,typ,k,which,sigma=sigma,conv=conv)
    
    def test_general_symmetric_modes(self):
        k=2
        for d in self.general:
            for typ in 'fd':
                for which in self.modes:
                    for sigma in (None,0.5):
                        for conv in _mattypes:
                            self.eval_evec(d,typ,k,which,sigma=sigma,conv=conv)

    def test_standard_symmetric_starting_vector(self):
        k=2
        for d in self.standard:
            for typ in 'fd':
                A=d['mat']
                n=A.shape[0]
                v0 = random.rand(n).astype(typ)
                self.eval_evec(d,typ,k,which='LM',v0=v0)

    def test_general_symmetric_starting_vector(self):
        k=2
        for d in self.general:
            for typ in 'fd':
                A=d['mat']
                n=A.shape[0]
                v0 = random.rand(n).astype(typ)
                self.eval_evec(d,typ,k,which='LM',v0=v0)

    def test_standard_symmetric_no_convergence(self):
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

class TestNonSymmetric(TestArpack):
    def setUp(self):
        self.eigs = eigs
        self.modes = ['LM','SM','LR', 'SR', 'LI', 'SI']

        #standard nonsymmetric test case A*x = w*x
        NS={}
        NS['mat']= array([[-2., -8.,  1.,  2., -5.],
                          [ 6.,  6.,  0.,  2.,  1.],
                          [ 0.,  4., -2., 11.,  0.],
                          [ 1.,  6.,  1.,  0., -4.],
                          [ 2., -6.,  4.,  9., -3]])
        
        NS['v0'] = array([0.39574246391875789,
                          0.00086496039750016962,
                          -0.9227205789982591,
                          -0.9165671495278005,
                          0.1175963848841306])
        NS['eval']=array([ -5.4854094033782888+0.0j,
                            -2.2169058544873783+8.5966096591588261j,
                            -2.2169058544873783-8.5966096591588261j,
                            4.4596105561765107+3.8007839204319454j,
                            4.4596105561765107-3.8007839204319454j],'D')
        self.standard = [NS]

        #general nonsymmetric test case A*x = w*B*x
        NG={}
        NG['mat'] =  NS['mat']

        NG['bmat'] = array([[11.,  0.,  0., -1., -2.],
                            [ 0.,  9.,  2.,  0.,  0.],
                            [ 0.,  2.,  1.,  0., -1.],
                            [-1.,  0.,  0., 10.,  1.],
                            [-2.,  0., -1.,  1., 13.]])

        NG['eval'] = array([-4.1690059397428412+0.j,
                             0.3168657188570727+0.4432894295551382j,
                             0.3168657188570727-0.4432894295551382j,
                             -0.5102635838808636+1.3352613149875483j,
                             -0.5102635838808637-1.3352613149875483j])
        
        NG['v0'] = array([0.39574246391875789,
                          0.00086496039750016962,
                          -0.9227205789982591,
                          -0.9165671495278005,
                          0.1175963848841306])
        self.general = [NG]

        #standard complex nonsymmetric case A*x = w*x
        NCS={}
        NCS['mat'] = NS['mat'] + 1j*(NS['mat']+1)
        NCS['eval'] = eig(NCS['mat'],left=False,right=False)
        NCS['v0'] = NG['v0']
        self.standard_complex = [NCS]

        #general complex nonsymmetric case A*x = w*M*x
        NCG={}
        NCG['mat'] = NCS['mat']
        NCG['bmat'] = NG['bmat']
        NCG['eval'] = eig(NCS['mat'],NCG['bmat'],left=False,right=False)
        NCG['v0'] = NG['v0']
        self.general_complex = [NCG]
        
    
    def test_standard_nonsymmetric_modes(self):
        k=2
        sigma = None
        for d in self.standard:
            for typ in 'fd':
                for which in self.modes:
                    for conv in _mattypes:
                        self.eval_evec(d,typ,k,which,sigma=sigma,conv=conv)
    
    def test_general_nonsymmetric_modes(self):
        k=2
        sigma = None
        for d in self.general:
            for typ in 'fd':
                for which in self.modes:
                    for conv in _mattypes:
                        self.eval_evec(d,typ,k,which,sigma=sigma,conv=conv)

    def test_complex_standard_nonsymmetric_modes(self):
        k=2
        for d in self.standard_complex:
            for typ in 'FD':
                for sigma in (None,0.1+0.1j):
                    for which in self.modes:
                        for conv in _mattypes:
                            self.eval_evec(d,typ,k,which,sigma=sigma,conv=conv)

    def test_complex_general_nonsymmetric_modes(self):
        k=2
        for d in self.general_complex:
            for typ in 'FD':
                for sigma in (None,0.1+0.1j):
                    for which in self.modes:
                        for conv in _mattypes:
                            self.eval_evec(d,typ,k,which,sigma=sigma,conv=conv)
        

    def test_standard_nonsymmetric_starting_vector(self):
        k=2
        sigma = None
        for d in self.standard:
            for typ in 'fd':
                A=d['mat']
                n=A.shape[0]
                v0 = random.rand(n).astype(typ)
                self.eval_evec(d,typ,k,which='LM',sigma=sigma,v0=v0)

    def test_general_nonsymmetric_starting_vector(self):
        k=2
        sigma = None
        for d in self.general:
            for typ in 'fd':
                A=d['mat']
                n=A.shape[0]
                v0 = random.rand(n).astype(typ)
                self.eval_evec(d,typ,k,which='LM',sigma=sigma,v0=v0)

    def test_standard_nonsymmetric_no_convergence(self):
        np.random.seed(1234)
        m = np.random.rand(30, 30)
        try:
            w, v = eigs(m, 4, which='LM', v0=m[:,0], maxiter=5)
            raise AssertionError("Spurious no-error exit")
        except ArpackNoConvergence, err:
            k = len(err.eigenvalues)
            if k <= 0:
                raise AssertionError("Spurious no-eigenvalues-found case")
            w, v = err.eigenvalues, err.eigenvectors
            for ww, vv in zip(w, v.T):
                assert_array_almost_equal(dot(m, vv), ww*vv,
                                          decimal=_ndigits['d'])

def test_eigen_bad_shapes():
    # A is not square.
    A = csc_matrix(np.zeros((2,3)))
    assert_raises(ValueError, eigs, A)


def test_eigen_bad_kwargs():
    # Test eigen on wrong keyword argument
    A = csc_matrix(np.zeros((2,2)))
    assert_raises(ValueError, eigs, A, which='XX')

def test_ticket_1459_arpack_crash():
    for dtype in [np.float32, np.float64]:
        # XXX: this test does not seem to catch the issue for float32,
        #      but we made the same fix there, just to be sure

        N = 6
        k = 2

        np.random.seed(2301)
        A = np.random.random((N, N)).astype(dtype)
        v0 = np.array([-0.71063568258907849895, -0.83185111795729227424,
                       -0.34365925382227402451, 0.46122533684552280420,
                       -0.58001341115969040629, -0.78844877570084292984e-01],
                      dtype=dtype)

        # Should not crash:
        evals, evecs = eigs(A, k, v0=v0)

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
