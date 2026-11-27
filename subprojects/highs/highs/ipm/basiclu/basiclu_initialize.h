lu_int basiclu_initialize
(
    lu_int m,
    lu_int istore[],
    double xstore[]
);

/*
Purpose:

    Initialize istore, xstore to a BASICLU instance. Set parameters to defaults
    and reset counters. The initialization fixes the dimension of matrices
    which can be processed by this instance.

    This routine must be called once before passing istore, xstore to any other
    basiclu_ routine.

Return:

    BASICLU_OK

        m, istore, xstore were valid arguments. Only in this case are istore,
        xstore initialized.

    BASICLU_ERROR_argument_missing

        istore or xstore is NULL.

    BASICLU_ERROR_invalid_argument

        m is less than or equal to zero.

Arguments:

    lu_int m

        The dimension of matrices which can be processed. m > 0.

    lu_int istore[]
    double xstore[]

        Fixed size arrays. These must be allocated by the user as follows:

          length of istore: BASICLU_SIZE_ISTORE_1 + BASICLU_SIZE_ISTORE_M * m
          length of xstore: BASICLU_SIZE_XSTORE_1 + BASICLU_SIZE_XSTORE_M * m

Info:

    After initialization, the following entries of xstore are maintained
    throughout by all basiclu_ routines:

    xstore[BASICLU_DIM] Matrix dimension (constant).

    xstore[BASICLU_NUPDATE] Number of updates since last factorization. This is
                            the sum of Forrest-Tomlin updates and permutation
                            updates.

    xstore[BASICLU_NFORREST] Number of Forrest-Tomlin updates since last
                             factorization. The upper limit on Forrest-Tomlin
                             updates before refactorization is m, but that is
                             far too much for performance reasons and numerical
                             stability.

    xstore[BASICLU_NFACTORIZE] Number of factorizations since initialization.

    xstore[BASICLU_NUPDATE_TOTAL] Number of updates since initialization.

    xstore[BASICLU_NFORREST_TOTAL] Number of Forrest-Tomlin updates since
                                   initialization.

    xstore[BASICLU_NSYMPERM_TOTAL] Number of symmetric permutation updates since
                                   initialization. A permutation update is
                                   "symmetric" if the row and column
                                   permutation can be updated symmetrically.

    xstore[BASICLU_LNZ] Number of nonzeros in L excluding diagonal elements
                        (not changed by updates).

    xstore[BASICLU_UNZ] Number of nonzeros in U excluding diagonal elements
                        (changed by updates).

    xstore[BASICLU_RNZ] Number of nonzeros in update ETA vectors excluding
                        diagonal elements (zero after factorization, increased
                        by Forrest-Tomlin updates).

    xstore[BASICLU_MIN_PIVOT]
    xstore[BASICLU_MAX_PIVOT] After factorization these are the smallest and
                              largest pivot element. xstore[BASICLU_MIN_PIVOT]
                              is replaced when a smaller pivot occurs in an
                              update. xstore[BASICLU_MAX_PIVOT] is replaced when
                              a larger pivot occurs in an update.

    xstore[BASICLU_UPDATE_COST] Deterministic measure of solve/update cost
                                compared to cost of last factorization. This
                                value is zero after factorization and
                                monotonically increases with solves/updates.
                                When xstore[BASICLU_UPDATE_COST] > 1.0, then
                                a refactorization is good for performance.

    xstore[BASICLU_TIME_FACTORIZE] Wall clock time for last factorization.

    xstore[BASICLU_TIME_SOLVE] Wall clock time for all calls to
                               basiclu_solve_sparse and basiclu_solve_for_update
                               since last factorization.

    xstore[BASICLU_TIME_UPDATE] Wall clock time for all calls to basiclu_update
                                since last factorization.

    xstore[BASICLU_TIME_FACTORIZE_TOTAL]
    xstore[BASICLU_TIME_SOLVE_TOTAL]
    xstore[BASICLU_TIME_UPDATE_TOTAL] Analogous to above, but summing up all
                                      calls since initialization.

    xstore[BASICLU_LFLOPS]
    xstore[BASICLU_UFLOPS]
    xstore[BASICLU_RFLOPS] Number of flops for operations with L, U and update
                           ETA vectors in calls to basiclu_solve_sparse and
                           basiclu_solve_for_update since last factorization.
*/
