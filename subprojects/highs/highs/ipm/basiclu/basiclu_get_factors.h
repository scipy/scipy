lu_int basiclu_get_factors
(
    lu_int istore[],
    double xstore[],
    lu_int Li[],
    double Lx[],
    lu_int Ui[],
    double Ux[],
    lu_int Wi[],
    double Wx[],
    lu_int rowperm[],
    lu_int colperm[],
    lu_int Lcolptr[],
    lu_int Lrowidx[],
    double Lvalue[],
    lu_int Ucolptr[],
    lu_int Urowidx[],
    double Uvalue[]
);

/*
Purpose:

    Extract the row and column permutation and the LU factors. This routine can
    be used only after basiclu_factorize() has completed and before a call to
    basiclu_update(). At that point the factorized form of matrix B is

        B[rowperm,colperm] = L*U,

    where L is unit lower triangular and U is upper triangular. If the
    factorization was singular (rank < m), then columns colperm[rank..m-1]
    of B have been replaced by unit columns with entry 1 in position
    rowperm[rank..m-1].

    basiclu_get_factors() is intended when the user needs direct access to the
    matrix factors. It is not required to solve linear systems with the factors
    (see basiclu_solve_dense() and basiclu_solve_sparse() instead).

Return:

    BASICLU_ERROR_invalid_store if istore, xstore do not hold a BASICLU
    instance. In this case xstore[BASICLU_STATUS] is not set.

    Otherwise return the status code. See xstore[BASICLU_STATUS] below.

Arguments:

    lu_int istore[]
    double xstore[]
    lu_int Li[]
    double Lx[]
    lu_int Ui[]
    double Ux[]
    lu_int Wi[]
    double Wx[]

        The BASICLU instance after basiclu_factorize() has completed.

    lu_int rowperm[m]

        Returns the row permutation. If the row permutation is not required,
        then NULL can be passed (this is not an error).

    lu_int colperm[m]

        Returns the column permutation. If the column permutation is not
        required, then NULL can be passed (this is not an error).

    lu_int Lcolptr[m+1]
    lu_int Lrowidx[m+Lnz]
    double Lvalue[m+Lnz], where Lnz = xstore[BASICLU_LNZ]

        If all three arguments are not NULL, then they are filled with L in
        compressed column form. The indices in each column are sorted with the
        unit diagonal element at the front.

        If any of the three arguments is NULL, then L is not returned
        (this is not an error).

    lu_int Ucolptr[m+1]
    lu_int Urowidx[m+Unz]
    double Uvalue[m+Unz], where Unz = xstore[BASICLU_UNZ]

        If all three arguments are not NULL, then they are filled with U in
        compressed column form. The indices in each column are sorted with the
        diagonal element at the end.

        If any of the three arguments is NULL, then U is not returned
        (this is not an error).

Info:

    xstore[BASICLU_STATUS]: status code.

        BASICLU_OK

            The requested quantities have been returned successfully.

        BASICLU_ERROR_argument_missing

            One or more of the mandatory pointer/array arguments are NULL.

        BASICLU_ERROR_invalid_call

            The BASICLU instance does not hold a fresh factorization (either
            basiclu_factorize() has not completed or basiclu_update() has been
            called in the meanwhile).
*/
