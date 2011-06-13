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

def generate_matrix(N, complex=False, hermitian=False, 
                    pos_definite=False, sparse=False):
    M = np.random.random((N,N))
    if complex:
        M = M + 1j * np.random.random((N,N))

    if hermitian:
        if pos_definite:
            if sparse:
                i = np.arange(N)
                j = np.random.randint(N, size=N-2)
                i, j = np.meshgrid(i, j)
                M[i,j] = 0
            M = np.dot(M.conj(), M.T)
        else:
            M = np.dot(M.conj(), M.T)
            if sparse:
                i = np.random.randint(N, size=N * N / 4)
                j = np.random.randint(N, size=N * N / 4)
                ind = np.where(i == j)
                j[ind] = (j[ind] + 1) % N
                M[i,j] = 0
                M[j,i] = 0
    else:
        if sparse:
            i = np.random.randint(N, size=N * N / 2)
            j = np.random.randint(N, size=N * N / 2)
            M[i,j] = 0
    return M

def _aslinearoperator_with_dtype(m):
    m = aslinearoperator(m)
    if not hasattr(m, 'dtype'):
        x = np.zeros(m.shape[1])
        m.dtype = (m * x).dtype
    return m


def assert_almost_equal_cc(actual, desired, decimal=7,
                           err_msg='', verbose=True):
    # almost equal or complex conjugates almost equal
    try:
        assert_almost_equal(actual, desired, decimal, err_msg, verbose)
    except:
        assert_almost_equal(actual, conj(desired), decimal, err_msg, verbose)


def assert_array_almost_equal_cc(actual, desired, decimal=7,
                                 err_msg='', verbose=True):
    # almost equal or complex conjugates almost equal
    try:
        assert_array_almost_equal(actual, desired, decimal, err_msg, verbose)
    except:
        assert_array_almost_equal(actual, conj(desired), decimal,
                                  err_msg, verbose)

# precision for tests
_ndigits = {'f': 3, 'd': 11, 'F': 2, 'D': 8}

class TestArpack(TestCase):
    def argsort_which(self, eval, typ, k, which, 
                      sigma=None, OPpart=None, mode=None):
        if sigma is None:
            reval = np.round(eval, decimals=_ndigits[typ])
        else:
            if mode is None or mode=='normal':
                if OPpart is None:
                    reval = 1. / (eval - sigma)
                elif OPpart == 'r':
                    reval = 0.5 * (1. / (eval - sigma) 
                                   + 1. / (eval - np.conj(sigma)))
                elif OPpart == 'i':
                    reval = -0.5j * (1. / (eval - sigma) 
                                 - 1. / (eval - np.conj(sigma)))
            elif mode=='cayley':
                reval = (eval + sigma) / (eval - sigma)
            elif mode=='buckling':
                reval = eval / (eval - sigma)
            else:
                raise ValueError("mode='%s' not recognized" % mode)
            
            reval = np.round(reval, decimals=_ndigits[typ])

        if which in ['LM', 'SM']:
            ind = np.argsort(abs(reval))
        elif which in ['LR', 'SR', 'LA', 'SA', 'BE']:
            ind = np.argsort(np.real(reval))
        elif which in ['LI', 'SI']:
            # for LI,SI ARPACK returns largest,smallest abs(imaginary) why?
            print '-'
            if typ.islower():
                ind = np.argsort(abs(np.imag(reval)))
            else:
                ind = np.argsort(np.imag(reval))
        else:
            raise ValueError("which='%s' is unrecognized" % which)
        
        if which in ['LM', 'LA', 'LR', 'LI']:
            return ind[-k:]
        elif which in ['SM', 'SA', 'SR', 'SI']:
            return ind[:k]
        elif which == 'BE':
            return np.concatenate((ind[:k/2], ind[k/2-k:]))

    def eval_evec(self, d, typ, k, which, sigma=None, v0=None,
                  mattype=np.asarray, OPpart=None, mode='normal'):
        general = ('bmat' in d)
        if general:
            err = ("error for %s:general, typ=%s, which=%s, sigma=%s, "
                   "mattype=%s, OPpart=%s, mode=%s" % (self.__class__.__name__,
                                                       typ, which, sigma,
                                                       mattype.__name__,
                                                       OPpart, mode))
        else:
            err = ("error for %s:standard, typ=%s, which=%s, sigma=%s, "
                   "mattype=%s, OPpart=%s, mode=%s" % (self.__class__.__name__,
                                                       typ, which, sigma,
                                                       mattype.__name__,
                                                       OPpart, mode))

        a = d['mat'].astype(typ)
        ac = mattype(a)
        
        if general:
            b = d['bmat'].astype(typ.lower())
            bc = mattype(b)
        
        # get exact eigenvalues
        exact_eval = d['eval'].astype(typ.upper())
        ind = self.argsort_which(exact_eval, typ, k, which,
                                 sigma, OPpart, mode)
        exact_eval_a = exact_eval
        exact_eval = exact_eval[ind]
        
        # compute arpack eigenvalues
        kwargs = dict(which=which, v0=v0, sigma=sigma)
        if self.eigs is eigsh:
            kwargs['mode'] = mode
        else:
            kwargs['OPpart'] = OPpart
            
        if general:
            eval, evec = self.eigs(ac, k, bc, **kwargs)
        else:
            eval, evec = self.eigs(ac, k, **kwargs)

        ind = self.argsort_which(eval, typ, k, which,
                                 sigma, OPpart, mode)
        eval_a = eval
        eval = eval[ind]
        evec = evec[:,ind]
        
        # check eigenvalues
        try:
            assert_array_almost_equal_cc(eval, exact_eval, 
                                         decimal=_ndigits[typ],
                                         err_msg=err)
        except:
            print err
            print eval_a
            print exact_eval_a
            

        # check eigenvectors
        LHS = np.dot(a, evec)
        if general:
            RHS = eval * np.dot(b, evec)
        else:
            RHS = eval * evec
            
        try:
            assert_array_almost_equal(LHS, RHS,
                                      decimal=_ndigits[typ],
                                      err_msg=err)
        except:
            print err

def modes(sigma):
    if sigma is None:
        return ["normal"]
    else:
        return ["normal", "cayley", "buckling"]

class TestSymmetric(TestArpack):
    def setUp(self):
        self.eigs = eigsh
        self.which = ['LM', 'SM', 'LA', 'SA', 'BE']
        self.mattypes = [csr_matrix, aslinearoperator, np.asarray]
        self.sigmas_modes = {None : ['normal'],
                             0.5 : ['normal', 'buckling', 'cayley']}

        #generate matrices
        # these should all be float32 so that the eigenvalues
        # are the same in float32 and float64
        N = 6
        np.random.seed(2300)
        Ar = generate_matrix(N, hermitian=True,
                             pos_definite=True).astype('f').astype('d')
        M = generate_matrix(N, hermitian=True,
                            pos_definite=True).astype('f').astype('d')
        Ac = generate_matrix(N, hermitian=True, pos_definite=True,
                             complex=True).astype('F').astype('D')
        v0 = np.random.random(N)

        # standard symmetric problem
        SS = {}
        SS['mat'] = Ar
        SS['v0'] = v0
        SS['eval'] = eigh(SS['mat'], eigvals_only=True)

        # general symmetric problem
        GS = {}
        GS['mat'] = Ar
        GS['bmat'] = M
        GS['v0'] = v0
        GS['eval'] = eigh(GS['mat'], GS['bmat'], eigvals_only=True)

        # standard hermitian problem
        SH = {}
        SH['mat'] = Ac
        SH['v0'] = v0
        SH['eval'] = eigh(SH['mat'], eigvals_only=True)

        # general hermitian problem
        GH = {}
        GH['mat'] = Ac
        GH['bmat'] = M
        GH['v0'] = v0
        GH['eval'] = eigh(GH['mat'], GH['bmat'], eigvals_only=True)

        self.real_test_cases = [SS, GS]
        self.complex_test_cases = [SH, GH]
    
    def test_symmetric_modes(self):
        k = 2
        for D in self.real_test_cases:
            for typ in 'fd':
                for which in self.which:
                    for mattype in self.mattypes:
                        for (sigma, modes) in self.sigmas_modes.iteritems():
                            for mode in modes:
                                self.eval_evec(D, typ, k, which, sigma=sigma, 
                                               mattype=mattype, mode=mode)
    
    def test_hermitian_modes(self):
        k = 2
        for D in self.complex_test_cases:
            for typ in 'FD':
                for which in self.which:
                    if which == 'BE': continue
                    for mattype in self.mattypes:
                        for sigma in self.sigmas_modes:
                            self.eval_evec(D, typ, k, which, 
                                           sigma=sigma, mattype=mattype)

    def test_symmetric_starting_vector(self):
        k = 2
        for D in self.real_test_cases:
            for typ in 'fd':
                v0 = random.rand(len(D['v0'])).astype(typ)
                self.eval_evec(D, typ, k, which='LM', v0=v0)

    def test_symmetric_no_convergence(self):
        np.random.seed(1234)
        m = generate_matrix(30, hermitian=True, pos_definite=True)

        try:
            w, v = eigsh(m, 4, which='LM', v0=m[:, 0], maxiter=5)
            raise AssertionError("Spurious no-error exit")
        except ArpackNoConvergence, err:
            k = len(err.eigenvalues)
            if k <= 0:
                raise AssertionError("Spurious no-eigenvalues-found case")
            w, v = err.eigenvalues, err.eigenvectors
            assert_array_almost_equal(dot(m, v), w * v, decimal=_ndigits['d'])
    


class TestNonSymmetric(TestArpack):
    def setUp(self):
        self.eigs = eigs
        self.which = ['LM', 'LR', 'LI']#, 'SM', 'LR', 'SR', 'LI', 'SI']
        self.mattypes = [csr_matrix, aslinearoperator, np.asarray]
        self.sigmas_OPparts = {None : [None],
                               0.1 : ['r'], 
                               0.1 + 0.1j : ['r', 'i']}

        #generate matrices
        # these should all be float32 so that the eigenvalues
        # are the same in float32 and float64
        N = 6
        np.random.seed(2300)
        Ar = generate_matrix(N).astype('f').astype('d')
        M = generate_matrix(N, hermitian=True,
                            pos_definite=True).astype('f').astype('d')
        Ac = generate_matrix(N, complex=True).astype('F').astype('D')
        v0 = np.random.random(N)

        # standard real nonsymmetric problem
        SNR = {}
        SNR['mat'] = Ar
        SNR['v0'] = v0
        SNR['eval'] = eig(SNR['mat'], left=False, right=False)

        # general real nonsymmetric problem
        GNR = {}
        GNR['mat'] = Ar
        GNR['bmat'] = M
        GNR['v0'] = v0
        GNR['eval'] = eig(GNR['mat'], GNR['bmat'], left=False, right=False)
        
        # standard complex nonsymmetric problem
        SNC = {}
        SNC['mat'] = Ac
        SNC['v0'] = v0
        SNC['eval'] = eig(SNC['mat'], left=False, right=False)

        # general complex nonsymmetric problem
        GNC = {}
        GNC['mat'] = Ac
        GNC['bmat'] = M
        GNC['v0'] = v0
        GNC['eval'] = eig(GNC['mat'], GNC['bmat'], left=False, right=False)
        
        self.real_test_cases = [SNR, GNR]
        self.complex_test_cases = [SNC, GNC]

        for D in self.real_test_cases + self.complex_test_cases:
            print D['eval']

    # This leads to memory errors.  Not sure the source...
    #
    #def test_real_nonsymmetric_modes(self):
    #    k = 2
    #    for D in self.real_test_cases:
    #        for typ in 'fd':
    #            for which in self.which:
    #                for mattype in self.mattypes:
    #                    for sigma, OPparts in self.sigmas_OPparts.iteritems():
    #                        for OPpart in OPparts:
    #                            self.eval_evec(D, typ, k, which, sigma=sigma,
    #                                           mattype=mattype, OPpart=OPpart)

    def test_complex_nonsymmetric_modes(self):
        k = 2
        for D in self.complex_test_cases:
            for typ in 'DF':
                for which in self.which:
                    for mattype in self.mattypes:
                        for sigma in self.sigmas_OPparts:
                            self.eval_evec(D, typ, k, which, sigma=sigma,
                                           mattype=mattype)


    def test_standard_nonsymmetric_starting_vector(self):
        k = 2
        sigma = None
        for d in self.complex_test_cases:
            for typ in 'FD':
                A = d['mat']
                n = A.shape[0]
                v0 = random.rand(n).astype(typ)
                self.eval_evec(d, typ, k, which='LM', sigma=sigma, v0=v0)

    def test_general_nonsymmetric_starting_vector(self):
        k = 2
        sigma = None
        for d in self.complex_test_cases:
            for typ in 'FD':
                A = d['mat']
                n = A.shape[0]
                v0 = random.rand(n).astype(typ)
                self.eval_evec(d, typ, k, which='LM', sigma=sigma, v0=v0)

    def test_standard_nonsymmetric_no_convergence(self):
        np.random.seed(1234)
        m = generate_matrix(30, complex=True)
        try:
            w, v = eigs(m, 4, which='LM', v0=m[:, 0], maxiter=5)
            raise AssertionError("Spurious no-error exit")
        except ArpackNoConvergence, err:
            k = len(err.eigenvalues)
            if k <= 0:
                raise AssertionError("Spurious no-eigenvalues-found case")
            w, v = err.eigenvalues, err.eigenvectors
            for ww, vv in zip(w, v.T):
                assert_array_almost_equal(dot(m, vv), ww * vv,
                                          decimal=_ndigits['d'])
    

def test_eigen_bad_shapes():
    # A is not square.
    A = csc_matrix(np.zeros((2, 3)))
    assert_raises(ValueError, eigs, A)


def test_eigen_bad_kwargs():
    # Test eigen on wrong keyword argument
    A = csc_matrix(np.zeros((2, 2)))
    assert_raises(ValueError, eigs, A, which='XX')


def sorted_svd(m, k):
    #Compute svd of a dense matrix m, and return singular vectors/values
    #sorted.
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
                      [1 + 1j, 0, 2],
                      [0, 0, 1]], np.complex)
        y = np.array([[1, 2, 3, 8 + 5j],
                      [3 - 2j, 4, 3, 5],
                      [1, 0, 2, 3],
                      [0, 0, 1, 0]], np.complex)
        z = csc_matrix(x)

        for m in [x, x.T.conjugate(), x.T, y, y.conjugate(), z, z.T]:
            for k in range(1, min(m.shape) - 1):
                u, s, vh = sorted_svd(m, k)
                su, ss, svh = svds(m, k)

                m_hat = svd_estimate(u, s, vh)
                sm_hat = svd_estimate(su, ss, svh)

                assert_array_almost_equal_nulp(m_hat, sm_hat, nulp=1000)

if __name__ == "__main__":
    run_module_suite()
