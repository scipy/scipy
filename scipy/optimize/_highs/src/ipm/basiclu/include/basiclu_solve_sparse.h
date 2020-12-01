lu_int basiclu_solve_sparse
(
    lu_int istore[],
    double xstore[],
    lu_int Li[],
    double Lx[],
    lu_int Ui[],
    double Ux[],
    lu_int Wi[],
    double Wx [],
    lu_int nzrhs,
    const lu_int irhs[],
    const double xrhs[],
    lu_int *p_nzlhs,
    lu_int ilhs[],
    double lhs[],
    char trans
);

/*
Purpose:

    Given the factorization computed by basiclu_factorize() or basiclu_update()
    and the sparse right-hand side, rhs, solve a linear system for the solution
    lhs.

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

        Factorization computed by basiclu_factorize() or basiclu_update().

    lu_int nzrhs
    const lu_int irhs[nzrhs]
    const double xrhs[nzrhs]

        The right-hand side vector in compressed format. irhs[0..nzrhs-1] are
        the indices of nonzeros and xrhs[0..nzrhs-1] the corresponding values.
        irhs must not contain duplicates.

    lu_int *p_nzlhs
    lu_int ilhs[m]
    lu_int lhs[m]

        *p_nzlhs is uninitialized on entry. On return *p_nzlhs holds
        the number of nonzeros in the solution.
        The contents of ilhs is uninitialized on entry. On return
        ilhs[0..*p_nzlhs-1] holds the indices of nonzeros in the solution.
        The contents lhs must be initialized to zero on entry. On return
        the solution is scattered into lhs.

    char trans

        Defines which system to solve. 't' or 'T' for the transposed system,
        any other character for the forward system.

Parameters:

    xstore[BASICLU_SPARSE_THRESHOLD]

        Defines which method is used for solving a triangular system. A
        triangular solve can be done either by the two phase method of Gilbert
        and Peierls ("sparse solve") or by a sequential pass through the vector
        ("sequential solve").

        Solving B*x=b requires two triangular solves. The first triangular solve
        is done sparse. The second triangular solve is done sparse if its
        right-hand side has not more than m * xstore[BASICLU_SPARSE_THRESHOLD]
        nonzeros. Otherwise the sequential solve is used.

        Default: 0.05

    xstore[BASICLU_DROP_TOLERANCE]

        Nonzeros which magnitude is less than or equal to the drop tolerance
        are removed after each triangular solve. Default: 1e-20

Info:

    xstore[BASICLU_STATUS]: status code.

        BASICLU_OK

            The linear system has been successfully solved.

        BASICLU_ERROR_argument_missing

            One or more of the pointer/array arguments are NULL.

        BASICLU_ERROR_invalid_call

            The factorization is invalid.

        BASICLU_ERROR_invalid_argument

            The right-hand side is invalid (nzrhs < 0 or nzrhs > m or one or
            more indices out of range).
*/
