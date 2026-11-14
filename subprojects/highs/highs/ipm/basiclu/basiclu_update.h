lu_int basiclu_update
(
    lu_int istore[],
    double xstore[],
    lu_int Li[],
    double Lx[],
    lu_int Ui[],
    double Ux[],
    lu_int Wi[],
    double Wx[],
    double xtbl
);

/*
Purpose:

    Update the factorization to replace one column of the factorized matrix.
    A call to basiclu_update() must be preceded by calls to
    basiclu_solve_for_update() to provide the column to be inserted and the
    index of the column to be replaced.

    The column to be inserted is defined as the right-hand side in the last call
    to basiclu_solve_for_update() in which the forward system was solved.

    The index of the column to be replaced is defined by the unit vector in the
    last call to basiclu_solve_for_update() in which the transposed system was
    solved.

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

    double xtbl

        This is an optional argument to monitor numerical stability. xtbl can be
        either of

        (a) element j0 of the solution to the forward system computed by
            basiclu_solve_for_update(), where j0 is the column to be replaced;

        (b) the dot product of the incoming column and the solution to the
            transposed system computed by basiclu_solve_for_update().

        In either case xstore[BASICLU_PIVOT_ERROR] (see below) has a defined
        value. If monitoring stability is not desired, xtbl can be any value.

Parameters:

    xstore[BASICLU_MEMORYL]: length of Li and Lx
    xstore[BASICLU_MEMORYU]: length of Ui and Ux
    xstore[BASICLU_MEMORYW]: length of Wi and Wx

    xstore[BASICLU_DROP_TOLERANCE]

        Nonzeros which magnitude is less than or equal to the drop tolerance
        are removed from the row eta matrix. Default: 1e-20

Info:

    xstore[BASICLU_STATUS]: status code.

        BASICLU_OK

            The update has successfully completed.

        BASICLU_ERROR_argument_missing

            One or more of the pointer/array arguments are NULL.

        BASICLU_ERROR_invalid_call

            The factorization is invalid or the update was not prepared by two
            calls to basiclu_solve_for_update().

        BASICLU_REALLOCATE

            Insufficient memory in Wi,Wx. The number of additional elements
            required is given by

                xstore[BASICLU_ADD_MEMORYW] > 0

            The user must reallocate Wi,Wx. It is recommended to reallocate for
            the requested number of additional elements plus some extra space
            for further updates (e.g. 0.5 times the current array length). The
            new array length must be provided in

                xstore[BASICLU_MEMORYW]: length of Wi and Wx

            basiclu_update will start from scratch in the next call.

        BASICLU_ERROR_singular_update

            The updated factorization would be (numerically) singular. No update
            has been computed and the old factorization is still valid.

    xstore[BASICLU_PIVOT_ERROR]

        When xtbl was given (see above), then xstore[BASICLU_PIVOT_ERROR] is a
        measure for numerical stability. It is the difference between two
        computations of the new pivot element relative to the new pivot element.
        A value larger than 1e-10 indicates numerical instability and suggests
        refactorization (and possibly tightening the pivot tolerance).

    xstore[BASICLU_MAX_ETA]

        The maximum entry (in absolute value) in the eta vectors from the
        Forrest-Tomlin update. A large value, say > 1e6, indicates that pivoting
        on diagonal element was unstable and refactorization might be necessary.
*/
