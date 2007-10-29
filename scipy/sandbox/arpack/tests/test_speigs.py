#!/usr/bin/env python

import sys
from numpy.testing import *
set_package_path()
from arpack.speigs import *
restore_path()

import numpy as N

class TestEigs(NumpyTestCase):
    def test(self):
        maxn=15                # Dimension of square matrix to be solved
        # Use a PDP^-1 factorisation to construct matrix with known
        # eiegevalues/vectors. Used random eiegenvectors initially.
        P = N.mat(N.random.random((maxn,)*2))
        P /= map(N.linalg.norm, P.T)            # Normalise the eigenvectors
        D = N.mat(N.zeros((maxn,)*2))
        D[range(maxn), range(maxn)] = (N.arange(maxn, dtype=float)+1)/N.sqrt(maxn)
        A = P*D*N.linalg.inv(P)
        vals = N.array(D.diagonal())[0]
        vecs = P
        uv_sortind = vals.argsort()
        vals = vals[uv_sortind]
        vecs = vecs[:,uv_sortind]

        from scipy.linalg.iterative import get_matvec
        matvec = get_matvec(A)
        #= lambda x: N.asarray(A*x)[0]
        nev=4
        eigvs = ARPACK_eigs(matvec, A.shape[0], nev=nev)
        calc_vals = eigvs[0]
        # Ensure the calculate eigenvectors have the same sign as the refence values
        calc_vecs = eigvs[1] / [N.sign(x[0]) for x in eigvs[1].T]
        assert_array_almost_equal(calc_vals, vals[0:nev], decimal=7)
        assert_array_almost_equal(calc_vecs,  N.array(vecs)[:,0:nev], decimal=7)


# class TestGeneigs(NumpyTestCase):
#     def test(self):
#         import pickle
#         import scipy.linsolve
#         A,B = pickle.load(file('mats.pickle'))
#         sigma = 27.
#         sigma_solve = scipy.linsolve.splu(A - sigma*B).solve
#         w = ARPACK_gen_eigs(B.matvec, sigma_solve, B.shape[0], sigma, 10)[0]
#         assert_array_almost_equal(w,
#         [27.346442255386375,  49.100299170945405,  56.508474856551544, 56.835800191692492,
#          65.944215785041365, 66.194792400328367, 78.003788872725238, 79.550811647295944,
#          94.646308846854879, 95.30841709116271], decimal=11)

if __name__ == "__main__":
    NumpyTest().run()
