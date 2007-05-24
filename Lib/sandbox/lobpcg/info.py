"""
The algorithm of LOBPCG is described in detail in:

A. V. Knyazev, Toward the Optimal Preconditioned Eigensolver: Locally Optimal
Block Preconditioned Conjugate Gradient Method. SIAM Journal on Scientific
Computing 23 (2001), no. 2,
pp. 517-541. http://dx.doi.org/10.1137/S1064827500366124

A. V. Knyazev, I. Lashuk, M. E. Argentati, and E. Ovchinnikov, Block Locally
Optimal Preconditioned Eigenvalue Xolvers (BLOPEX) in hypre and PETSc
(2007). http://arxiv.org/abs/0705.2626

Call the function lobpcg - see help for lobpcg.lobpcg. See also lobpcg.as2d,
which can be used in the preconditioner (example below)

Example:

    # Solve A x = lambda B x with constraints and preconditioning.

    n = 100

    vals = [nm.arange( n, dtype = nm.float64 ) + 1]

    # Matrix A.
    operatorA = spdiags( vals, 0, n, n )

    # Matrix B
    operatorB = nm.eye( n, n )

    # Constraints.
    Y = nm.eye( n, 3 )

    # Initial guess for eigenvectors, should have linearly independent
    # columns. Column dimension = number of requested eigenvalues.
    X = sc.rand( n, 3 )

    # Preconditioner - inverse of A.
    ivals = [1./vals[0]]
    def precond( x ):
        invA = spdiags( ivals, 0, n, n )
        y = invA  * x
        if sp.issparse( y ):
            y = y.toarray()

        return as2d( y )

    # Alternative way of providing the same preconditioner.
    #precond = spdiags( ivals, 0, n, n )

    tt = time.clock()
    eigs, vecs = lobpcg( X, operatorA, operatorB, blockVectorY = Y,
                         operatorT = precond,
                         residualTolerance = 1e-4, maxIterations = 40,
                         largest = False, verbosityLevel = 1 )
    print 'solution time:', time.clock() - tt
    print eigs

Usage notes:

Notation: n - matrix size, m - number of required eigenvalues (smallest or
largest)

1) The LOBPCG code internally solves eigenproblems of the size 3m on every
   iteration by calling the"standard" eigensolver, so if m is not small enough
   compared to n, it does not make sense to call the LOBPCG code, but rather
   one should use the "standard" eigensolver, e.g. symeig function in this
   case. If one calls the LOBPCG algorithm for 5m>n, it will most likely break
   internally, so the code tries to call symeig instead.

   It is not that n should be large for the LOBPCG to work, but rather the
   ratio n/m should be large. It you call the LOBPCG code with m=1 and n=10, it
   should work, though n is small. The method is intended for extremely large
   n/m, see e.g., reference [28] in http://arxiv.org/abs/0705.2626

2) The convergence speed depends basically on two factors:

  a) how well relatively separated the seeking eigenvalues are from the rest of
  the eigenvalues. One can try to vary m to make this better.

  b) how "well conditioned" the problem is. This can be changed by using proper
  preconditioning. For example, a rod vibration test problem (under tests
  directory) is ill-conditioned
  for large n, so convergence will be slow, unless efficient preconditioning is
  used. For this specific problem, a good simple preconditioner function would
  be a linear solve for A, which is easy to code since A is tridiagonal.

"""
postpone_import = 1
