"""
The algorithm of LOBPCG is described in detail in:

A. V. Knyazev, Toward the Optimal Preconditioned Eigensolver: Locally Optimal Block Preconditioned Conjugate Gradient Method. SIAM Journal on Scientific Computing 23 (2001), no. 2, pp. 517-541. http://dx.doi.org/10.1137/S1064827500366124

A. V. Knyazev, I. Lashuk, M. E. Argentati, and E. Ovchinnikov, Block Locally Optimal Preconditioned Eigenvalue Xolvers (BLOPEX) in hypre and PETSc  (2007). http://arxiv.org/abs/0705.2626


Depends upon symeig (http://mdp-toolkit.sourceforge.net/symeig.html) for the
moment, as the symmetric eigenvalue solvers were not available in scipy.

Usage: XXXXX


"""
postpone_import = 1
