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
  python -c 'import scipy;scipy.linalg.test()'
Run tests if linalg is not installed:
  python tests/test_decomp.py
"""

import numpy as np
from numpy.testing import TestCase, assert_equal, assert_array_almost_equal, \
        assert_array_equal, assert_raises, assert_, run_module_suite, dec

from scipy.linalg import eig, eigvals, lu, svd, svdvals, cholesky, qr, \
     schur, rsf2csf, lu_solve, lu_factor, solve, diagsvd, hessenberg, rq, \
     eig_banded, eigvals_banded, eigh, eigvalsh, LinAlgError
from scipy.linalg.flapack import dgbtrf, dgbtrs, zgbtrf, zgbtrs, \
     dsbev, dsbevd, dsbevx, zhbevd, zhbevx

from numpy import array, transpose, sometrue, diag, ones, linalg, \
     argsort, zeros, arange, float32, complex64, dot, conj, identity, \
     ravel, sqrt, iscomplex, shape, sort, conjugate, bmat, sign, \
     asarray, matrix, isfinite, all, ndarray, outer, eye, dtype, empty,\
     triu, tril

from numpy.random import rand, normal, seed

from scipy.linalg._testutils import assert_no_overwrite

# digit precision to use in asserts for different types
DIGITS = {'d':11, 'D':11, 'f':4, 'F':4}

# XXX: This function should be available through numpy.testing
def assert_dtype_equal(act, des):
    if isinstance(act, ndarray):
        act = act.dtype
    else:
        act = dtype(act)

    if isinstance(des, ndarray):
        des = des.dtype
    else:
        des = dtype(des)

    assert_(act == des, 'dtype mismatch: "%s" (should be "%s") ' % (act, des))

# XXX: This function should not be defined here, but somewhere in
#      scipy.linalg namespace
def symrand(dim_or_eigv):
    """Return a random symmetric (Hermitian) matrix.

    If 'dim_or_eigv' is an integer N, return a NxN matrix, with eigenvalues
        uniformly distributed on (-1,1).

    If 'dim_or_eigv' is  1-D real array 'a', return a matrix whose
                      eigenvalues are 'a'.
    """
    if isinstance(dim_or_eigv, int):
        dim = dim_or_eigv
        d = (rand(dim)*2)-1
    elif (isinstance(dim_or_eigv, ndarray) and
          len(dim_or_eigv.shape) == 1):
        dim = dim_or_eigv.shape[0]
        d = dim_or_eigv
    else:
        raise TypeError("input type not supported.")

    v = random_rot(dim)
    h = dot(dot(v.T.conj(), diag(d)), v)
    # to avoid roundoff errors, symmetrize the matrix (again)
    h = 0.5*(h.T+h)
    return h

# XXX: This function should not be defined here, but somewhere in
#      scipy.linalg namespace
def random_rot(dim):
    """Return a random rotation matrix, drawn from the Haar distribution
    (the only uniform distribution on SO(n)).
    The algorithm is described in the paper
    Stewart, G.W., 'The efficient generation of random orthogonal
    matrices with an application to condition estimators', SIAM Journal
    on Numerical Analysis, 17(3), pp. 403-409, 1980.
    For more information see
    http://en.wikipedia.org/wiki/Orthogonal_matrix#Randomization"""
    H = eye(dim)
    D = ones((dim, ))
    for n in range(1, dim):
        x = normal(size=(dim-n+1, ))
        D[n-1] = sign(x[0])
        x[0] -= D[n-1]*sqrt((x*x).sum())
        # Householder transformation

        Hx = eye(dim-n+1) - 2.*outer(x, x)/(x*x).sum()
        mat = eye(dim)
        mat[n-1:,n-1:] = Hx
        H = dot(H, mat)
    # Fix the last sign such that the determinant is 1
    D[-1] = -D.prod()
    H = (D*H.T).T
    return H

def random(size):
    return rand(*size)

class TestEigVals(TestCase):

    def test_simple(self):
        a = [[1,2,3],[1,2,3],[2,5,6]]
        w = eigvals(a)
        exact_w = [(9+sqrt(93))/2,0,(9-sqrt(93))/2]
        assert_array_almost_equal(w,exact_w)

    def test_simple_tr(self):
        a = array([[1,2,3],[1,2,3],[2,5,6]],'d')
        a = transpose(a).copy()
        a = transpose(a)
        w = eigvals(a)
        exact_w = [(9+sqrt(93))/2,0,(9-sqrt(93))/2]
        assert_array_almost_equal(w,exact_w)

    def test_simple_complex(self):
        a = [[1,2,3],[1,2,3],[2,5,6+1j]]
        w = eigvals(a)
        exact_w = [(9+1j+sqrt(92+6j))/2,
                   0,
                   (9+1j-sqrt(92+6j))/2]
        assert_array_almost_equal(w,exact_w)


class TestEig(object):

    def test_simple(self):
        a = [[1,2,3],[1,2,3],[2,5,6]]
        w,v = eig(a)
        exact_w = [(9+sqrt(93))/2,0,(9-sqrt(93))/2]
        v0 = array([1,1,(1+sqrt(93)/3)/2])
        v1 = array([3.,0,-1])
        v2 = array([1,1,(1-sqrt(93)/3)/2])
        v0 = v0 / sqrt(dot(v0,transpose(v0)))
        v1 = v1 / sqrt(dot(v1,transpose(v1)))
        v2 = v2 / sqrt(dot(v2,transpose(v2)))
        assert_array_almost_equal(w,exact_w)
        assert_array_almost_equal(v0,v[:,0]*sign(v[0,0]))
        assert_array_almost_equal(v1,v[:,1]*sign(v[0,1]))
        assert_array_almost_equal(v2,v[:,2]*sign(v[0,2]))
        for i in range(3):
            assert_array_almost_equal(dot(a,v[:,i]),w[i]*v[:,i])
        w,v = eig(a,left=1,right=0)
        for i in range(3):
            assert_array_almost_equal(dot(transpose(a),v[:,i]),w[i]*v[:,i])

    def test_simple_complex_eig(self):
        a = [[1,2],[-2,1]]
        w,vl,vr = eig(a,left=1,right=1)
        assert_array_almost_equal(w, array([1+2j, 1-2j]))
        for i in range(2):
            assert_array_almost_equal(dot(a,vr[:,i]),w[i]*vr[:,i])
        for i in range(2):
            assert_array_almost_equal(dot(conjugate(transpose(a)),vl[:,i]),
                                      conjugate(w[i])*vl[:,i])

    def test_simple_complex(self):
        a = [[1,2,3],[1,2,3],[2,5,6+1j]]
        w,vl,vr = eig(a,left=1,right=1)
        for i in range(3):
            assert_array_almost_equal(dot(a,vr[:,i]),w[i]*vr[:,i])
        for i in range(3):
            assert_array_almost_equal(dot(conjugate(transpose(a)),vl[:,i]),
                                      conjugate(w[i])*vl[:,i])

    def _check_gen_eig(self, A, B):
        A, B = asarray(A), asarray(B)
        msg = "\n%r\n%r" % (A, B)
        w, vr = eig(A,B)
        wt = eigvals(A,B)
        val1 = dot(A, vr)
        val2 = dot(B, vr) * w
        res = val1 - val2
        for i in range(res.shape[1]):
            if all(isfinite(res[:, i])):
                assert_array_almost_equal(res[:, i], 0, err_msg=msg)

        assert_array_almost_equal(sort(w[isfinite(w)]), sort(wt[isfinite(wt)]),
                                  err_msg=msg)

    def test_singular(self):
        """Test singular pair"""
        # Example taken from
        # http://www.cs.umu.se/research/nla/singular_pairs/guptri/matlab.html
        A = array(( [22,34,31,31,17], [45,45,42,19,29], [39,47,49,26,34],
            [27,31,26,21,15], [38,44,44,24,30]))
        B = array(( [13,26,25,17,24], [31,46,40,26,37], [26,40,19,25,25],
            [16,25,27,14,23], [24,35,18,21,22]))

        olderr = np.seterr(all='ignore')
        try:
            self._check_gen_eig(A, B)
        finally:
            np.seterr(**olderr)

    def test_falker(self):
        """Test matrices giving some Nan generalized eigen values."""
        M = diag(array(([1,0,3])))
        K = array(([2,-1,-1],[-1,2,-1],[-1,-1,2]))
        D = array(([1,-1,0],[-1,1,0],[0,0,0]))
        Z = zeros((3,3))
        I = identity(3)
        A = bmat([[I,Z],[Z,-K]])
        B = bmat([[Z,I],[M,D]])

        olderr = np.seterr(all='ignore')
        try:
            self._check_gen_eig(A, B)
        finally:
            np.seterr(**olderr)

    def test_bad_geneig(self):
        # Ticket #709 (strange return values from DGGEV)

        def matrices(omega):
            c1 = -9 + omega**2
            c2 = 2*omega
            A = [[1, 0,  0,  0],
                 [0, 1,  0,  0],
                 [0, 0,  c1, 0],
                 [0, 0,  0, c1]]
            B = [[0, 0,  1,   0],
                 [0, 0,  0,   1],
                 [1, 0,  0, -c2],
                 [0, 1, c2,   0]]
            return A, B

        # With a buggy LAPACK, this can fail for different omega on different
        # machines -- so we need to test several values
        olderr = np.seterr(all='ignore')
        try:
            for k in xrange(100):
                A, B = matrices(omega=k*5./100)
                self._check_gen_eig(A, B)
        finally:
            np.seterr(**olderr)

    def test_not_square_error(self):
        """Check that passing a non-square array raises a ValueError."""
        A = np.arange(6).reshape(3,2)
        assert_raises(ValueError, eig, A)

    def test_shape_mismatch(self):
        """Check that passing arrays of with different shapes raises a ValueError."""
        A = identity(2)
        B = np.arange(9.0).reshape(3,3)
        assert_raises(ValueError, eig, A, B)
        assert_raises(ValueError, eig, B, A)

class TestEigBanded(TestCase):

    def __init__(self, *args):
        TestCase.__init__(self, *args)

        self.create_bandmat()

    def create_bandmat(self):
        """Create the full matrix `self.fullmat` and
           the corresponding band matrix `self.bandmat`."""
        N  = 10
        self.KL = 2   # number of subdiagonals (below the diagonal)
        self.KU = 2   # number of superdiagonals (above the diagonal)

        # symmetric band matrix
        self.sym_mat = ( diag(1.0*ones(N))
                     +  diag(-1.0*ones(N-1), -1) + diag(-1.0*ones(N-1), 1)
                     + diag(-2.0*ones(N-2), -2) + diag(-2.0*ones(N-2), 2) )

        # hermitian band matrix
        self.herm_mat = ( diag(-1.0*ones(N))
                     + 1j*diag(1.0*ones(N-1), -1) - 1j*diag(1.0*ones(N-1), 1)
                     + diag(-2.0*ones(N-2), -2) + diag(-2.0*ones(N-2), 2) )

        # general real band matrix
        self.real_mat = ( diag(1.0*ones(N))
                     +  diag(-1.0*ones(N-1), -1) + diag(-3.0*ones(N-1), 1)
                     + diag(2.0*ones(N-2), -2) + diag(-2.0*ones(N-2), 2) )

        # general complex band matrix
        self.comp_mat = ( 1j*diag(1.0*ones(N))
                     +  diag(-1.0*ones(N-1), -1) + 1j*diag(-3.0*ones(N-1), 1)
                     + diag(2.0*ones(N-2), -2) + diag(-2.0*ones(N-2), 2) )


        # Eigenvalues and -vectors from linalg.eig
        ew, ev = linalg.eig(self.sym_mat)
        ew = ew.real
        args = argsort(ew)
        self.w_sym_lin = ew[args]
        self.evec_sym_lin = ev[:,args]

        ew, ev = linalg.eig(self.herm_mat)
        ew = ew.real
        args = argsort(ew)
        self.w_herm_lin = ew[args]
        self.evec_herm_lin = ev[:,args]


        # Extract upper bands from symmetric and hermitian band matrices
        # (for use in dsbevd, dsbevx, zhbevd, zhbevx
        #  and their single precision versions)
        LDAB = self.KU + 1
        self.bandmat_sym  = zeros((LDAB, N), dtype=float)
        self.bandmat_herm = zeros((LDAB, N), dtype=complex)
        for i in xrange(LDAB):
            self.bandmat_sym[LDAB-i-1,i:N]  = diag(self.sym_mat, i)
            self.bandmat_herm[LDAB-i-1,i:N] = diag(self.herm_mat, i)


        # Extract bands from general real and complex band matrix
        # (for use in dgbtrf, dgbtrs and their single precision versions)
        LDAB = 2*self.KL + self.KU + 1
        self.bandmat_real = zeros((LDAB, N), dtype=float)
        self.bandmat_real[2*self.KL,:] = diag(self.real_mat)     # diagonal
        for i in xrange(self.KL):
            # superdiagonals
            self.bandmat_real[2*self.KL-1-i,i+1:N]   = diag(self.real_mat, i+1)
            # subdiagonals
            self.bandmat_real[2*self.KL+1+i,0:N-1-i] = diag(self.real_mat,-i-1)

        self.bandmat_comp = zeros((LDAB, N), dtype=complex)
        self.bandmat_comp[2*self.KL,:] = diag(self.comp_mat)     # diagonal
        for i in xrange(self.KL):
            # superdiagonals
            self.bandmat_comp[2*self.KL-1-i,i+1:N]   = diag(self.comp_mat, i+1)
            # subdiagonals
            self.bandmat_comp[2*self.KL+1+i,0:N-1-i] = diag(self.comp_mat,-i-1)

        # absolute value for linear equation system A*x = b
        self.b = 1.0*arange(N)
        self.bc = self.b *(1 + 1j)


    #####################################################################


    def test_dsbev(self):
        """Compare dsbev eigenvalues and eigenvectors with
           the result of linalg.eig."""
        w, evec, info  = dsbev(self.bandmat_sym, compute_v=1)
        evec_ = evec[:,argsort(w)]
        assert_array_almost_equal(sort(w), self.w_sym_lin)
        assert_array_almost_equal(abs(evec_), abs(self.evec_sym_lin))



    def test_dsbevd(self):
        """Compare dsbevd eigenvalues and eigenvectors with
           the result of linalg.eig."""
        w, evec, info = dsbevd(self.bandmat_sym, compute_v=1)
        evec_ = evec[:,argsort(w)]
        assert_array_almost_equal(sort(w), self.w_sym_lin)
        assert_array_almost_equal(abs(evec_), abs(self.evec_sym_lin))



    def test_dsbevx(self):
        """Compare dsbevx eigenvalues and eigenvectors
           with the result of linalg.eig."""
        N,N = shape(self.sym_mat)
        ## Achtung: Argumente 0.0,0.0,range?
        w, evec, num, ifail, info = dsbevx(self.bandmat_sym, 0.0, 0.0, 1, N,
                                       compute_v=1, range=2)
        evec_ = evec[:,argsort(w)]
        assert_array_almost_equal(sort(w), self.w_sym_lin)
        assert_array_almost_equal(abs(evec_), abs(self.evec_sym_lin))


    def test_zhbevd(self):
        """Compare zhbevd eigenvalues and eigenvectors
           with the result of linalg.eig."""
        w, evec, info = zhbevd(self.bandmat_herm, compute_v=1)
        evec_ = evec[:,argsort(w)]
        assert_array_almost_equal(sort(w), self.w_herm_lin)
        assert_array_almost_equal(abs(evec_), abs(self.evec_herm_lin))



    def test_zhbevx(self):
        """Compare zhbevx eigenvalues and eigenvectors
           with the result of linalg.eig."""
        N,N = shape(self.herm_mat)
        ## Achtung: Argumente 0.0,0.0,range?
        w, evec, num, ifail, info = zhbevx(self.bandmat_herm, 0.0, 0.0, 1, N,
                                       compute_v=1, range=2)
        evec_ = evec[:,argsort(w)]
        assert_array_almost_equal(sort(w), self.w_herm_lin)
        assert_array_almost_equal(abs(evec_), abs(self.evec_herm_lin))



    def test_eigvals_banded(self):
        """Compare eigenvalues of eigvals_banded with those of linalg.eig."""
        w_sym = eigvals_banded(self.bandmat_sym)
        w_sym = w_sym.real
        assert_array_almost_equal(sort(w_sym), self.w_sym_lin)

        w_herm = eigvals_banded(self.bandmat_herm)
        w_herm = w_herm.real
        assert_array_almost_equal(sort(w_herm), self.w_herm_lin)

        # extracting eigenvalues with respect to an index range
        ind1 = 2
        ind2 = 6
        w_sym_ind = eigvals_banded(self.bandmat_sym,
                                    select='i', select_range=(ind1, ind2) )
        assert_array_almost_equal(sort(w_sym_ind),
                                  self.w_sym_lin[ind1:ind2+1])
        w_herm_ind = eigvals_banded(self.bandmat_herm,
                                    select='i', select_range=(ind1, ind2) )
        assert_array_almost_equal(sort(w_herm_ind),
                                  self.w_herm_lin[ind1:ind2+1])

        # extracting eigenvalues with respect to a value range
        v_lower = self.w_sym_lin[ind1] - 1.0e-5
        v_upper = self.w_sym_lin[ind2] + 1.0e-5
        w_sym_val = eigvals_banded(self.bandmat_sym,
                                select='v', select_range=(v_lower, v_upper) )
        assert_array_almost_equal(sort(w_sym_val),
                                  self.w_sym_lin[ind1:ind2+1])

        v_lower = self.w_herm_lin[ind1] - 1.0e-5
        v_upper = self.w_herm_lin[ind2] + 1.0e-5
        w_herm_val = eigvals_banded(self.bandmat_herm,
                                select='v', select_range=(v_lower, v_upper) )
        assert_array_almost_equal(sort(w_herm_val),
                                  self.w_herm_lin[ind1:ind2+1])



    def test_eig_banded(self):
        """Compare eigenvalues and eigenvectors of eig_banded
           with those of linalg.eig. """
        w_sym, evec_sym = eig_banded(self.bandmat_sym)
        evec_sym_ = evec_sym[:,argsort(w_sym.real)]
        assert_array_almost_equal(sort(w_sym), self.w_sym_lin)
        assert_array_almost_equal(abs(evec_sym_), abs(self.evec_sym_lin))

        w_herm, evec_herm = eig_banded(self.bandmat_herm)
        evec_herm_ = evec_herm[:,argsort(w_herm.real)]
        assert_array_almost_equal(sort(w_herm), self.w_herm_lin)
        assert_array_almost_equal(abs(evec_herm_), abs(self.evec_herm_lin))

        # extracting eigenvalues with respect to an index range
        ind1 = 2
        ind2 = 6
        w_sym_ind, evec_sym_ind = eig_banded(self.bandmat_sym,
                                    select='i', select_range=(ind1, ind2) )
        assert_array_almost_equal(sort(w_sym_ind),
                                  self.w_sym_lin[ind1:ind2+1])
        assert_array_almost_equal(abs(evec_sym_ind),
                                  abs(self.evec_sym_lin[:,ind1:ind2+1]) )

        w_herm_ind, evec_herm_ind = eig_banded(self.bandmat_herm,
                                    select='i', select_range=(ind1, ind2) )
        assert_array_almost_equal(sort(w_herm_ind),
                                  self.w_herm_lin[ind1:ind2+1])
        assert_array_almost_equal(abs(evec_herm_ind),
                                  abs(self.evec_herm_lin[:,ind1:ind2+1]) )

        # extracting eigenvalues with respect to a value range
        v_lower = self.w_sym_lin[ind1] - 1.0e-5
        v_upper = self.w_sym_lin[ind2] + 1.0e-5
        w_sym_val, evec_sym_val = eig_banded(self.bandmat_sym,
                                select='v', select_range=(v_lower, v_upper) )
        assert_array_almost_equal(sort(w_sym_val),
                                  self.w_sym_lin[ind1:ind2+1])
        assert_array_almost_equal(abs(evec_sym_val),
                                  abs(self.evec_sym_lin[:,ind1:ind2+1]) )

        v_lower = self.w_herm_lin[ind1] - 1.0e-5
        v_upper = self.w_herm_lin[ind2] + 1.0e-5
        w_herm_val, evec_herm_val = eig_banded(self.bandmat_herm,
                                select='v', select_range=(v_lower, v_upper) )
        assert_array_almost_equal(sort(w_herm_val),
                                  self.w_herm_lin[ind1:ind2+1])
        assert_array_almost_equal(abs(evec_herm_val),
                                  abs(self.evec_herm_lin[:,ind1:ind2+1]) )


    def test_dgbtrf(self):
        """Compare dgbtrf  LU factorisation with the LU factorisation result
           of linalg.lu."""
        M,N = shape(self.real_mat)
        lu_symm_band, ipiv, info = dgbtrf(self.bandmat_real, self.KL, self.KU)

        # extract matrix u from lu_symm_band
        u = diag(lu_symm_band[2*self.KL,:])
        for i in xrange(self.KL + self.KU):
            u += diag(lu_symm_band[2*self.KL-1-i,i+1:N], i+1)

        p_lin, l_lin, u_lin = lu(self.real_mat, permute_l=0)
        assert_array_almost_equal(u, u_lin)


    def test_zgbtrf(self):
        """Compare zgbtrf  LU factorisation with the LU factorisation result
           of linalg.lu."""
        M,N = shape(self.comp_mat)
        lu_symm_band, ipiv, info = zgbtrf(self.bandmat_comp, self.KL, self.KU)

        # extract matrix u from lu_symm_band
        u = diag(lu_symm_band[2*self.KL,:])
        for i in xrange(self.KL + self.KU):
            u += diag(lu_symm_band[2*self.KL-1-i,i+1:N], i+1)

        p_lin, l_lin, u_lin =lu(self.comp_mat, permute_l=0)
        assert_array_almost_equal(u, u_lin)



    def test_dgbtrs(self):
        """Compare dgbtrs  solutions for linear equation system  A*x = b
           with solutions of linalg.solve."""

        lu_symm_band, ipiv, info = dgbtrf(self.bandmat_real, self.KL, self.KU)
        y, info = dgbtrs(lu_symm_band, self.KL, self.KU, self.b, ipiv)

        y_lin = linalg.solve(self.real_mat, self.b)
        assert_array_almost_equal(y, y_lin)

    def test_zgbtrs(self):
        """Compare zgbtrs  solutions for linear equation system  A*x = b
           with solutions of linalg.solve."""

        lu_symm_band, ipiv, info = zgbtrf(self.bandmat_comp, self.KL, self.KU)
        y, info = zgbtrs(lu_symm_band, self.KL, self.KU, self.bc, ipiv)

        y_lin = linalg.solve(self.comp_mat, self.bc)
        assert_array_almost_equal(y, y_lin)

def test_eigh():
    DIM = 6
    v = {'dim': (DIM, ),
         'dtype': ('f','d','F','D'),
         'overwrite': (True, False),
         'lower': (True, False),
         'turbo': (True, False),
         'eigvals': (None, (2, DIM-2))}

    for dim in v['dim']:
        for typ in v['dtype']:
            for overwrite in v['overwrite']:
                for turbo in v['turbo']:
                    for eigvals in v['eigvals']:
                        for lower in v['lower']:
                            yield (eigenhproblem_standard,
                                   'ordinary',
                                   dim, typ, overwrite, lower,
                                   turbo, eigvals)
                            yield (eigenhproblem_general,
                                   'general ',
                                   dim, typ, overwrite, lower,
                                   turbo, eigvals)

def _complex_symrand(dim, dtype):
    a1, a2 = symrand(dim), symrand(dim)
    # add antisymmetric matrix as imag part
    a = a1 +1j*(triu(a2)-tril(a2))
    return a.astype(dtype)

def eigenhproblem_standard(desc, dim, dtype,
                           overwrite, lower, turbo,
                           eigvals):
    """Solve a standard eigenvalue problem."""
    if iscomplex(empty(1, dtype=dtype)):
        a = _complex_symrand(dim, dtype)
    else:
        a = symrand(dim).astype(dtype)

    if overwrite:
        a_c = a.copy()
    else:
        a_c = a
    w, z = eigh(a, overwrite_a=overwrite, lower=lower, eigvals=eigvals)
    assert_dtype_equal(z.dtype, dtype)
    w = w.astype(dtype)
    diag_ = diag(dot(z.T.conj(), dot(a_c, z))).real
    assert_array_almost_equal(diag_, w, DIGITS[dtype])

def eigenhproblem_general(desc, dim, dtype,
                          overwrite, lower, turbo,
                          eigvals):
    """Solve a generalized eigenvalue problem."""
    if iscomplex(empty(1, dtype=dtype)):
        a = _complex_symrand(dim, dtype)
        b = _complex_symrand(dim, dtype)+diag([2.1]*dim).astype(dtype)
    else:
        a = symrand(dim).astype(dtype)
        b = symrand(dim).astype(dtype)+diag([2.1]*dim).astype(dtype)

    if overwrite:
        a_c, b_c = a.copy(), b.copy()
    else:
        a_c, b_c = a, b

    w, z = eigh(a, b, overwrite_a=overwrite, lower=lower,
                overwrite_b=overwrite, turbo=turbo, eigvals=eigvals)
    assert_dtype_equal(z.dtype, dtype)
    w = w.astype(dtype)
    diag1_ = diag(dot(z.T.conj(), dot(a_c, z))).real
    assert_array_almost_equal(diag1_, w, DIGITS[dtype])
    diag2_ = diag(dot(z.T.conj(), dot(b_c, z))).real
    assert_array_almost_equal(diag2_, ones(diag2_.shape[0]), DIGITS[dtype])

def test_eigh_integer():
    a = array([[1,2],[2,7]])
    b = array([[3,1],[1,5]])
    w,z = eigh(a)
    w,z = eigh(a,b)

class TestLU(TestCase):

    def __init__(self, *args, **kw):
        TestCase.__init__(self, *args, **kw)

        self.a = array([[1,2,3],[1,2,3],[2,5,6]])
        self.ca = array([[1,2,3],[1,2,3],[2,5j,6]])
        # Those matrices are more robust to detect problems in permutation
        # matrices than the ones above
        self.b = array([[1,2,3],[4,5,6],[7,8,9]])
        self.cb = array([[1j,2j,3j],[4j,5j,6j],[7j,8j,9j]])

        # Reectangular matrices
        self.hrect = array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 12, 12]])
        self.chrect = 1.j * array([[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 12, 12]])

        self.vrect = array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 12, 12]])
        self.cvrect = 1.j * array([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 12, 12]])

        # Medium sizes matrices
        self.med = rand(30, 40)
        self.cmed = rand(30, 40) + 1.j * rand(30, 40)

    def _test_common(self, data):
        p,l,u = lu(data)
        assert_array_almost_equal(dot(dot(p,l),u),data)
        pl,u = lu(data,permute_l=1)
        assert_array_almost_equal(dot(pl,u),data)

    # Simple tests
    def test_simple(self):
        self._test_common(self.a)

    def test_simple_complex(self):
        self._test_common(self.ca)

    def test_simple2(self):
        self._test_common(self.b)

    def test_simple2_complex(self):
        self._test_common(self.cb)

    # rectangular matrices tests
    def test_hrectangular(self):
        self._test_common(self.hrect)

    def test_vrectangular(self):
        self._test_common(self.vrect)

    def test_hrectangular_complex(self):
        self._test_common(self.chrect)

    def test_vrectangular_complex(self):
        self._test_common(self.cvrect)

    # Bigger matrices
    def test_medium1(self):
        """Check lu decomposition on medium size, rectangular matrix."""
        self._test_common(self.med)

    def test_medium1_complex(self):
        """Check lu decomposition on medium size, rectangular matrix."""
        self._test_common(self.cmed)

class TestLUSingle(TestLU):
    """LU testers for single precision, real and double"""
    def __init__(self, *args, **kw):
        TestLU.__init__(self, *args, **kw)

        self.a = self.a.astype(float32)
        self.ca = self.ca.astype(complex64)
        self.b = self.b.astype(float32)
        self.cb = self.cb.astype(complex64)

        self.hrect = self.hrect.astype(float32)
        self.chrect = self.hrect.astype(complex64)

        self.vrect = self.vrect.astype(float32)
        self.cvrect = self.vrect.astype(complex64)

        self.med = self.vrect.astype(float32)
        self.cmed = self.vrect.astype(complex64)

class TestLUSolve(TestCase):
    def setUp(self):
        seed(1234)

    def test_lu(self):
        a = random((10,10))
        b = random((10,))

        x1 = solve(a,b)

        lu_a = lu_factor(a)
        x2 = lu_solve(lu_a,b)

        assert_array_equal(x1,x2)

class TestSVD(TestCase):
    def setUp(self):
        seed(1234)

    def test_simple(self):
        a = [[1,2,3],[1,20,3],[2,5,6]]
        u,s,vh = svd(a)
        assert_array_almost_equal(dot(transpose(u),u),identity(3))
        assert_array_almost_equal(dot(transpose(vh),vh),identity(3))
        sigma = zeros((u.shape[0],vh.shape[0]),s.dtype.char)
        for i in range(len(s)): sigma[i,i] = s[i]
        assert_array_almost_equal(dot(dot(u,sigma),vh),a)

    def test_simple_singular(self):
        a = [[1,2,3],[1,2,3],[2,5,6]]
        u,s,vh = svd(a)
        assert_array_almost_equal(dot(transpose(u),u),identity(3))
        assert_array_almost_equal(dot(transpose(vh),vh),identity(3))
        sigma = zeros((u.shape[0],vh.shape[0]),s.dtype.char)
        for i in range(len(s)): sigma[i,i] = s[i]
        assert_array_almost_equal(dot(dot(u,sigma),vh),a)

    def test_simple_underdet(self):
        a = [[1,2,3],[4,5,6]]
        u,s,vh = svd(a)
        assert_array_almost_equal(dot(transpose(u),u),identity(2))
        assert_array_almost_equal(dot(transpose(vh),vh),identity(3))
        sigma = zeros((u.shape[0],vh.shape[0]),s.dtype.char)
        for i in range(len(s)): sigma[i,i] = s[i]
        assert_array_almost_equal(dot(dot(u,sigma),vh),a)

    def test_simple_overdet(self):
        a = [[1,2],[4,5],[3,4]]
        u,s,vh = svd(a)
        assert_array_almost_equal(dot(transpose(u),u),identity(3))
        assert_array_almost_equal(dot(transpose(vh),vh),identity(2))
        sigma = zeros((u.shape[0],vh.shape[0]),s.dtype.char)
        for i in range(len(s)): sigma[i,i] = s[i]
        assert_array_almost_equal(dot(dot(u,sigma),vh),a)

    def test_random(self):
        n = 20
        m = 15
        for i in range(3):
            for a in [random([n,m]),random([m,n])]:
                u,s,vh = svd(a)
                assert_array_almost_equal(dot(transpose(u),u),identity(len(u)))
                assert_array_almost_equal(dot(transpose(vh),vh),identity(len(vh)))
                sigma = zeros((u.shape[0],vh.shape[0]),s.dtype.char)
                for i in range(len(s)): sigma[i,i] = s[i]
                assert_array_almost_equal(dot(dot(u,sigma),vh),a)

    def test_simple_complex(self):
        a = [[1,2,3],[1,2j,3],[2,5,6]]
        u,s,vh = svd(a)
        assert_array_almost_equal(dot(conj(transpose(u)),u),identity(3))
        assert_array_almost_equal(dot(conj(transpose(vh)),vh),identity(3))
        sigma = zeros((u.shape[0],vh.shape[0]),s.dtype.char)
        for i in range(len(s)): sigma[i,i] = s[i]
        assert_array_almost_equal(dot(dot(u,sigma),vh),a)

    def test_random_complex(self):
        n = 20
        m = 15
        for i in range(3):
            for a in [random([n,m]),random([m,n])]:
                a = a + 1j*random(list(a.shape))
                u,s,vh = svd(a)
                assert_array_almost_equal(dot(conj(transpose(u)),u),identity(len(u)))
                # This fails when [m,n]
                #assert_array_almost_equal(dot(conj(transpose(vh)),vh),identity(len(vh),dtype=vh.dtype.char))
                sigma = zeros((u.shape[0],vh.shape[0]),s.dtype.char)
                for i in range(len(s)): sigma[i,i] = s[i]
                assert_array_almost_equal(dot(dot(u,sigma),vh),a)

class TestSVDVals(TestCase):

    def test_simple(self):
        a = [[1,2,3],[1,2,3],[2,5,6]]
        s = svdvals(a)
        assert_(len(s) == 3)
        assert_(s[0] >= s[1] >= s[2])

    def test_simple_underdet(self):
        a = [[1,2,3],[4,5,6]]
        s = svdvals(a)
        assert_(len(s) == 2)
        assert_(s[0] >= s[1])

    def test_simple_overdet(self):
        a = [[1,2],[4,5],[3,4]]
        s = svdvals(a)
        assert_(len(s) == 2)
        assert_(s[0] >= s[1])

    def test_simple_complex(self):
        a = [[1,2,3],[1,20,3j],[2,5,6]]
        s = svdvals(a)
        assert_(len(s) == 3)
        assert_(s[0] >= s[1] >= s[2])

    def test_simple_underdet_complex(self):
        a = [[1,2,3],[4,5j,6]]
        s = svdvals(a)
        assert_(len(s) == 2)
        assert_(s[0] >= s[1])

    def test_simple_overdet_complex(self):
        a = [[1,2],[4,5],[3j,4]]
        s = svdvals(a)
        assert_(len(s) == 2)
        assert_(s[0] >= s[1])

class TestDiagSVD(TestCase):

    def test_simple(self):
        assert_array_almost_equal(diagsvd([1,0,0],3,3),[[1,0,0],[0,0,0],[0,0,0]])


class TestQR(TestCase):

    def setUp(self):
        seed(1234)

    def test_simple(self):
        a = [[8,2,3],[2,9,3],[5,3,6]]
        q,r = qr(a)
        assert_array_almost_equal(dot(transpose(q),q),identity(3))
        assert_array_almost_equal(dot(q,r),a)

    def test_simple_trap(self):
        a = [[8,2,3],[2,9,3]]
        q,r = qr(a)
        assert_array_almost_equal(dot(transpose(q),q),identity(2))
        assert_array_almost_equal(dot(q,r),a)

    def test_simple_tall(self):
        # full version
        a = [[8,2],[2,9],[5,3]]
        q,r = qr(a)
        assert_array_almost_equal(dot(transpose(q),q),identity(3))
        assert_array_almost_equal(dot(q,r),a)

    def test_simple_tall_e(self):
        # economy version
        a = [[8,2],[2,9],[5,3]]
        q,r = qr(a, mode='economic')
        assert_array_almost_equal(dot(transpose(q),q),identity(2))
        assert_array_almost_equal(dot(q,r),a)
        assert_equal(q.shape, (3,2))
        assert_equal(r.shape, (2,2))

    def test_simple_fat(self):
        # full version
        a = [[8,2,5],[2,9,3]]
        q,r = qr(a)
        assert_array_almost_equal(dot(transpose(q),q),identity(2))
        assert_array_almost_equal(dot(q,r),a)
        assert_equal(q.shape, (2,2))
        assert_equal(r.shape, (2,3))

    def test_simple_fat_e(self):
        # economy version
        a = [[8,2,3],[2,9,5]]
        q,r = qr(a, mode='economic')
        assert_array_almost_equal(dot(transpose(q),q),identity(2))
        assert_array_almost_equal(dot(q,r),a)
        assert_equal(q.shape, (2,2))
        assert_equal(r.shape, (2,3))

    def test_simple_complex(self):
        a = [[3,3+4j,5],[5,2,2+7j],[3,2,7]]
        q,r = qr(a)
        assert_array_almost_equal(dot(conj(transpose(q)),q),identity(3))
        assert_array_almost_equal(dot(q,r),a)

    def test_random(self):
        n = 20
        for k in range(2):
            a = random([n,n])
            q,r = qr(a)
            assert_array_almost_equal(dot(transpose(q),q),identity(n))
            assert_array_almost_equal(dot(q,r),a)

    def test_random_tall(self):
        # full version
        m = 200
        n = 100
        for k in range(2):
            a = random([m,n])
            q,r = qr(a)
            assert_array_almost_equal(dot(transpose(q),q),identity(m))
            assert_array_almost_equal(dot(q,r),a)

    def test_random_tall_e(self):
        # economy version
        m = 200
        n = 100
        for k in range(2):
            a = random([m,n])
            q,r = qr(a, mode='economic')
            assert_array_almost_equal(dot(transpose(q),q),identity(n))
            assert_array_almost_equal(dot(q,r),a)
            assert_equal(q.shape, (m,n))
            assert_equal(r.shape, (n,n))

    def test_random_trap(self):
        m = 100
        n = 200
        for k in range(2):
            a = random([m,n])
            q,r = qr(a)
            assert_array_almost_equal(dot(transpose(q),q),identity(m))
            assert_array_almost_equal(dot(q,r),a)

    def test_random_complex(self):
        n = 20
        for k in range(2):
            a = random([n,n])+1j*random([n,n])
            q,r = qr(a)
            assert_array_almost_equal(dot(conj(transpose(q)),q),identity(n))
            assert_array_almost_equal(dot(q,r),a)

class TestRQ(TestCase):

    def setUp(self):
        seed(1234)

    def test_simple(self):
        a = [[8,2,3],[2,9,3],[5,3,6]]
        r,q = rq(a)
        assert_array_almost_equal(dot(q, transpose(q)),identity(3))
        assert_array_almost_equal(dot(r,q),a)

    def test_r(self):
        a = [[8,2,3],[2,9,3],[5,3,6]]
        r,q = rq(a)
        r2 = rq(a, mode='r')
        assert_array_almost_equal(r, r2)

    def test_random(self):
        n = 20
        for k in range(2):
            a = random([n,n])
            r,q = rq(a)
            assert_array_almost_equal(dot(q, transpose(q)),identity(n))
            assert_array_almost_equal(dot(r,q),a)

    def test_simple_trap(self):
        a = [[8,2,3],[2,9,3]]
        r,q = rq(a)
        assert_array_almost_equal(dot(transpose(q),q),identity(3))
        assert_array_almost_equal(dot(r,q),a)

    def test_simple_tall(self):
        a = [[8,2],[2,9],[5,3]]
        r,q = rq(a)
        assert_array_almost_equal(dot(transpose(q),q),identity(2))
        assert_array_almost_equal(dot(r,q),a)

    def test_simple_fat(self):
        a = [[8,2,5],[2,9,3]]
        r,q = rq(a)
        assert_array_almost_equal(dot(transpose(q),q),identity(3))
        assert_array_almost_equal(dot(r,q),a)

    def test_simple_complex(self):
        a = [[3,3+4j,5],[5,2,2+7j],[3,2,7]]
        r,q = rq(a)
        assert_array_almost_equal(dot(q, conj(transpose(q))),identity(3))
        assert_array_almost_equal(dot(r,q),a)

    def test_random_tall(self):
        m = 200
        n = 100
        for k in range(2):
            a = random([m,n])
            r,q = rq(a)
            assert_array_almost_equal(dot(q, transpose(q)),identity(n))
            assert_array_almost_equal(dot(r,q),a)

    def test_random_trap(self):
        m = 100
        n = 200
        for k in range(2):
            a = random([m,n])
            r,q = rq(a)
            assert_array_almost_equal(dot(q, transpose(q)),identity(n))
            assert_array_almost_equal(dot(r,q),a)

    def test_random_trap_economic(self):
        m = 100
        n = 200
        for k in range(2):
            a = random([m,n])
            r,q = rq(a, mode='economic')
            assert_array_almost_equal(dot(q,transpose(q)),identity(m))
            assert_array_almost_equal(dot(r,q),a)
            assert_equal(q.shape, (m, n))
            assert_equal(r.shape, (m, m))

    def test_random_complex(self):
        n = 20
        for k in range(2):
            a = random([n,n])+1j*random([n,n])
            r,q = rq(a)
            assert_array_almost_equal(dot(q, conj(transpose(q))),identity(n))
            assert_array_almost_equal(dot(r,q),a)

    def test_random_complex_economic(self):
        m = 100
        n = 200
        for k in range(2):
            a = random([m,n])+1j*random([m,n])
            r,q = rq(a, mode='economic')
            assert_array_almost_equal(dot(q,conj(transpose(q))),identity(m))
            assert_array_almost_equal(dot(r,q),a)
            assert_equal(q.shape, (m, n))
            assert_equal(r.shape, (m, m))

transp = transpose
any = sometrue

class TestSchur(TestCase):

    def test_simple(self):
        a = [[8,12,3],[2,9,3],[10,3,6]]
        t,z = schur(a)
        assert_array_almost_equal(dot(dot(z,t),transp(conj(z))),a)
        tc,zc = schur(a,'complex')
        assert_(any(ravel(iscomplex(zc))) and any(ravel(iscomplex(tc))))
        assert_array_almost_equal(dot(dot(zc,tc),transp(conj(zc))),a)
        tc2,zc2 = rsf2csf(tc,zc)
        assert_array_almost_equal(dot(dot(zc2,tc2),transp(conj(zc2))),a)
        
    def test_sort(self):
        a = [[4.,3.,1.,-1.],[-4.5,-3.5,-1.,1.],[9.,6.,-4.,4.5],[6.,4.,-3.,3.5]]
        s,u,sdim = schur(a,sort='lhp')
        assert_array_almost_equal([[0.1134,0.5436,0.8316,0.], 
                                   [-0.1134,-0.8245,0.5544,0.], 
                                   [-0.8213,0.1308,0.0265,-0.5547], 
                                   [-0.5475,0.0872,0.0177,0.8321]],
                                  u,3)
        assert_array_almost_equal([[-1.4142,0.1456,-11.5816,-7.7174], 
                                   [0.,-0.5000,9.4472,-0.7184], 
                                   [0.,0.,1.4142,-0.1456], 
                                   [0.,0.,0.,0.5]],
                                  s,3)
        assert_equal(2,sdim)
                                  
        s,u,sdim = schur(a,sort='rhp')
        assert_array_almost_equal([[0.4862,-0.4930,0.1434,-0.7071],
                                   [-0.4862,0.4930,-0.1434,-0.7071],
                                   [0.6042,0.3944,-0.6924,0.],
                                   [0.4028,0.5986,0.6924,0.]],
                                  u,3)
        assert_array_almost_equal([[1.4142,-0.9270,4.5368,-14.4130],
                                   [0.,0.5,6.5809,-3.1870],
                                   [0.,0.,-1.4142,0.9270],
                                   [0.,0.,0.,-0.5]],
                                  s,3)
        assert_equal(2,sdim)
                                  
        s,u,sdim = schur(a,sort='iuc')
        assert_array_almost_equal([[0.5547,0.,-0.5721,-0.6042],
                                   [-0.8321,0.,-0.3814,-0.4028],
                                   [0.,0.7071,-0.5134,0.4862],
                                   [0.,0.7071,0.5134,-0.4862]],
                                  u,3)
        assert_array_almost_equal([[-0.5000,0.0000,-6.5809,-4.0974],
                                   [0.,0.5000,-3.3191,-14.4130],
                                   [0.,0.,1.4142,2.1573],
                                   [0.,0.,0.,-1.4142]],
                                  s,3)
        assert_equal(2,sdim)
                                  
        s,u,sdim = schur(a,sort='ouc')
        assert_array_almost_equal([[0.4862,-0.5134,0.7071,0.],
                                   [-0.4862,0.5134,0.7071,0.],
                                   [0.6042,0.5721,0.,-0.5547],
                                   [0.4028,0.3814,0.,0.8321]],
                                  u,3)
        assert_array_almost_equal([[1.4142,-2.1573,14.4130,4.0974],
                                   [0.,-1.4142,3.3191,6.5809],
                                   [0.,0.,-0.5000,0.],
                                   [0.,0.,0.,0.5000]],
                                  s,3)
        assert_equal(2,sdim)
                                  
        rhp_function = lambda x: x >= 0.0
        s,u,sdim = schur(a,sort=rhp_function)
        assert_array_almost_equal([[0.4862,-0.4930,0.1434,-0.7071],
                                   [-0.4862,0.4930,-0.1434,-0.7071],
                                   [0.6042,0.3944,-0.6924,0.],
                                   [0.4028,0.5986,0.6924,0.]],
                                  u,3)
        assert_array_almost_equal([[1.4142,-0.9270,4.5368,-14.4130],
                                   [0.,0.5,6.5809,-3.1870],
                                   [0.,0.,-1.4142,0.9270],
                                   [0.,0.,0.,-0.5]],
                                  s,3)
        assert_equal(2,sdim)
                                  
    def test_sort_errors(self):
        a = [[4.,3.,1.,-1.],[-4.5,-3.5,-1.,1.],[9.,6.,-4.,4.5],[6.,4.,-3.,3.5]]
        assert_raises(ValueError, schur, a, sort='unsupported')
        assert_raises(ValueError, schur, a, sort=1)

class TestHessenberg(TestCase):

    def test_simple(self):
        a = [[-149, -50,-154],
             [ 537, 180, 546],
             [ -27,  -9, -25]]
        h1 = [[-149.0000,42.2037,-156.3165],
              [-537.6783,152.5511,-554.9272],
              [0,0.0728, 2.4489]]
        h,q = hessenberg(a,calc_q=1)
        assert_array_almost_equal(dot(transp(q),dot(a,q)),h)
        assert_array_almost_equal(h,h1,decimal=4)

    def test_simple_complex(self):
        a = [[-149, -50,-154],
             [ 537, 180j, 546],
             [ -27j,  -9, -25]]
        h,q = hessenberg(a,calc_q=1)
        h1 = dot(transp(conj(q)),dot(a,q))
        assert_array_almost_equal(h1,h)

    def test_simple2(self):
        a = [[1,2,3,4,5,6,7],
             [0,2,3,4,6,7,2],
             [0,2,2,3,0,3,2],
             [0,0,2,8,0,0,2],
             [0,3,1,2,0,1,2],
             [0,1,2,3,0,1,0],
             [0,0,0,0,0,1,2]]
        h,q = hessenberg(a,calc_q=1)
        assert_array_almost_equal(dot(transp(q),dot(a,q)),h)

    def test_random(self):
        n = 20
        for k in range(2):
            a = random([n,n])
            h,q = hessenberg(a,calc_q=1)
            assert_array_almost_equal(dot(transp(q),dot(a,q)),h)

    def test_random_complex(self):
        n = 20
        for k in range(2):
            a = random([n,n])+1j*random([n,n])
            h,q = hessenberg(a,calc_q=1)
            h1 = dot(transp(conj(q)),dot(a,q))
            assert_array_almost_equal(h1,h)



class TestDatacopied(TestCase):

    def test_datacopied(self):
        from scipy.linalg.decomp import _datacopied

        M = matrix([[0,1],[2,3]])
        A = asarray(M)
        L = M.tolist()
        M2 = M.copy()

        class Fake1:
            def __array__(self):
                return A

        class Fake2:
            __array_interface__ = A.__array_interface__

        F1 = Fake1()
        F2 = Fake2()

        AF1 = asarray(F1)
        AF2 = asarray(F2)

        for item, status in [(M, False), (A, False), (L, True),
                             (M2, False), (F1, False), (F2, False)]:
            arr = asarray(item)
            assert_equal(_datacopied(arr, item), status,
                         err_msg=repr(item))


def test_aligned_mem_float():
    """Check linalg works with non-aligned memory"""
    # Allocate 402 bytes of memory (allocated on boundary)
    a = arange(402, dtype=np.uint8)

    # Create an array with boundary offset 4
    z = np.frombuffer(a.data, offset=2, count=100, dtype=float32)
    z.shape = 10, 10

    eig(z, overwrite_a=True)
    eig(z.T, overwrite_a=True)


def test_aligned_mem():
    """Check linalg works with non-aligned memory"""
    # Allocate 804 bytes of memory (allocated on boundary)
    a = arange(804, dtype=np.uint8)

    # Create an array with boundary offset 4
    z = np.frombuffer(a.data, offset=4, count=100, dtype=float)
    z.shape = 10, 10

    eig(z, overwrite_a=True)
    eig(z.T, overwrite_a=True)

def test_aligned_mem_complex():
    """Check that complex objects don't need to be completely aligned"""
    # Allocate 1608 bytes of memory (allocated on boundary)
    a = zeros(1608, dtype=np.uint8)

    # Create an array with boundary offset 8
    z = np.frombuffer(a.data, offset=8, count=100, dtype=complex)
    z.shape = 10, 10

    eig(z, overwrite_a=True)
    # This does not need special handling
    eig(z.T, overwrite_a=True)

def check_lapack_misaligned(func, args, kwargs):
    args = list(args)
    for i in range(len(args)):
        a = args[:]
        if isinstance(a[i],np.ndarray):
            # Try misaligning a[i]
            aa = np.zeros(a[i].size*a[i].dtype.itemsize+8, dtype=np.uint8)
            aa = np.frombuffer(aa.data, offset=4, count=a[i].size, dtype=a[i].dtype)
            aa.shape = a[i].shape
            aa[...] = a[i]
            a[i] = aa
            func(*a,**kwargs)
            if len(a[i].shape)>1:
                a[i] = a[i].T
                func(*a,**kwargs)


@dec.knownfailureif(True, "Ticket #1152, triggers a segfault in rare cases.")
def test_lapack_misaligned():
    M = np.eye(10,dtype=float)
    R = np.arange(100)
    R.shape = 10,10
    S = np.arange(20000,dtype=np.uint8)
    S = np.frombuffer(S.data, offset=4, count=100, dtype=np.float)
    S.shape = 10, 10
    b = np.ones(10)
    v = np.ones(3,dtype=float)
    LU, piv = lu_factor(S)
    for (func, args, kwargs) in [
            (eig,(S,),dict(overwrite_a=True)), # crash
            (eigvals,(S,),dict(overwrite_a=True)), # no crash
            (lu,(S,),dict(overwrite_a=True)), # no crash
            (lu_factor,(S,),dict(overwrite_a=True)), # no crash
            (lu_solve,((LU,piv),b),dict(overwrite_b=True)),
            (solve,(S,b),dict(overwrite_a=True,overwrite_b=True)),
            (svd,(M,),dict(overwrite_a=True)), # no crash
            (svd,(R,),dict(overwrite_a=True)), # no crash
            (svd,(S,),dict(overwrite_a=True)), # crash
            (svdvals,(S,),dict()), # no crash
            (svdvals,(S,),dict(overwrite_a=True)), #crash
            (cholesky,(M,),dict(overwrite_a=True)), # no crash
            (qr,(S,),dict(overwrite_a=True)), # crash
            (rq,(S,),dict(overwrite_a=True)), # crash
            (hessenberg,(S,),dict(overwrite_a=True)), # crash
            (schur,(S,),dict(overwrite_a=True)), # crash
            ]:
        yield check_lapack_misaligned, func, args, kwargs
# not properly tested
# cholesky, rsf2csf, lu_solve, solve, eig_banded, eigvals_banded, eigh, diagsvd


class TestOverwrite(object):
    def test_eig(self):
        assert_no_overwrite(eig, [(3,3)])
        assert_no_overwrite(eig, [(3,3), (3,3)])
    def test_eigh(self):
        assert_no_overwrite(eigh, [(3,3)])
        assert_no_overwrite(eigh, [(3,3), (3,3)])
    def test_eig_banded(self):
        assert_no_overwrite(eig_banded, [(3,2)])
    def test_eigvals(self):
        assert_no_overwrite(eigvals, [(3,3)])
    def test_eigvalsh(self):
        assert_no_overwrite(eigvalsh, [(3,3)])
    def test_eigvals_banded(self):
        assert_no_overwrite(eigvals_banded, [(3,2)])
    def test_hessenberg(self):
        assert_no_overwrite(hessenberg, [(3,3)])
    def test_lu_factor(self):
        assert_no_overwrite(lu_factor, [(3,3)])
    def test_lu_solve(self):
        x = np.array([[1,2,3], [4,5,6], [7,8,8]])
        xlu = lu_factor(x)
        assert_no_overwrite(lambda b: lu_solve(xlu, b), [(3,)])
    def test_lu(self):
        assert_no_overwrite(lu, [(3,3)])
    def test_qr(self):
        assert_no_overwrite(qr, [(3,3)])
    def test_rq(self):
        assert_no_overwrite(rq, [(3,3)])
    def test_schur(self):
        assert_no_overwrite(schur, [(3,3)])
    def test_schur_complex(self):
        assert_no_overwrite(lambda a: schur(a, 'complex'), [(3,3)],
                            dtypes=[np.float32, np.float64])
    def test_svd(self):
        assert_no_overwrite(svd, [(3,3)])
    def test_svdvals(self):
        assert_no_overwrite(svdvals, [(3,3)])

if __name__ == "__main__":
    run_module_suite()
