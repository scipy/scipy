lu_int basiclu_solve_for_update
(
    lu_int istore[],
    double xstore[],
    lu_int Li[],
    double Lx[],
    lu_int Ui[],
    double Ux[],
    lu_int Wi[],
    double Wx[],
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

    Given the factorization computed by basiclu_factorize() or basiclu_update(),
    solve a linear system in preparation to update the factorization.

    When the forward system is solved, then the right-hand side is the column
    to be inserted into the factorized matrix. When the transposed system is
    solved, then the right-hand side is a unit vector with entry 1 in position
    of the column to be replaced in the factorized matrix.

    For BASICLU to prepare the update, it is sufficient to compute only a
    partial solution. If the left-hand side is not requested by the user (see
    below), then only one triangular solve is done. If the left-hand side is
    requested, then a second triangular solve is required.

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

        The right-hand side vector in compressed format.

        When the forward system is solved, irhs[0..nzrhs-1] are the indices of
        nonzeros and xrhs[0..nzrhs-1] the corresponding values. irhs must not
        contain duplicates.
        
        When the transposed system is solved, the right-hand side is a unit
        vector with entry 1 in position irhs[0]. nzrhs, xrhs and elements of
        irhs other than irhs[0] are not accessed. xrhs can be NULL.

    lu_int *p_nzlhs
    lu_int ilhs[m]
    lu_int lhs[m]

        If any of p_nzlhs, ilhs or lhs is NULL, then the solution to the linear
        system is not requested. In this case only the update is prepared.

        Otherwise:

        *p_nzlhs is uninitialized on entry. On return *p_nzlhs holds
        the number of nonzeros in the solution.
        The contents of ilhs is uninitialized on entry. On return
        ilhs[0..*p_nzlhs-1] holds the indices of nonzeros in the solution.
        The contents of lhs must be initialized to zero on entry. On return
        the solution is  scattered into lhs.

    char trans

        Defines which system to solve. 't' or 'T' for the transposed system,
        any other character for the forward system.

Parameters:

    xstore[BASICLU_MEMORYL]: length of Li and Lx
    xstore[BASICLU_MEMORYU]: length of Ui and Ux
    xstore[BASICLU_MEMORYW]: length of Wi and Wx

    xstore[BASICLU_SPARSE_THRESHOLD]

        Defines which method is used for solving a triangular system. A
        triangular solve can be done either by the two phase method of Gilbert
        and Peierls ("sparse solve") or by a sequential pass through the vector
        ("sequential solve").

        When the solution to the linear system is requested, then two triangular
        systems are solved. The first triangular solve is done sparse. The
        second triangular solve is done sparse if its right-hand side has not
        more than m * xstore[BASICLU_SPARSE_THRESHOLD] nonzeros. Otherwise the
        sequential solve is used.

        When the solution to the linear system is not requested, then this
        parameter has no effect.

        Default: 0.05

    xstore[BASICLU_DROP_TOLERANCE]

        Nonzeros which magnitude is less than or equal to the drop tolerance
        are removed after each triangular solve. Default: 1e-20

Info:

    xstore[BASICLU_STATUS]: status code.

        BASICLU_OK

            The updated has been successfully prepared and, if requested, the
            solution to the linear system has been computed.

        BASICLU_ERROR_argument_missing

            One or more of the mandatory pointer/array arguments are NULL.

        BASICLU_ERROR_invalid_call

            The factorization is invalid.

        BASICLU_ERROR_maximum_updates

            There have already been m Forrest-Tomlin updates, see
            xstore[BASICLU_NFORREST]. The factorization cannot be updated any
            more and must be recomputed by basiclu_factorize().
            The solution to the linear system has not been computed.

        BASICLU_ERROR_invalid_argument

            The right-hand side is invalid (forward system: nzrhs < 0 or
            nzrhs > m or one or more indices out of range; backward system:
            irhs[0] out of range).

        BASICLU_REALLOCATE

            The solve was aborted because of insufficient memory in Li,Lx or
            Ui,Ux to store data for basiclu_update(). The number of additional
            elements required is given by

                xstore[BASICLU_ADD_MEMORYL] >= 0
                xstore[BASICLU_ADD_MEMORYU] >= 0

            The user must reallocate the arrays for which additional memory is
            required. It is recommended to reallocate for the requested number
            of additional elements plus some extra space for further updates
            (e.g. 0.5 times the current array length). The new array lengths
            must be provided in

                xstore[BASICLU_MEMORYL]: length of Li and Lx
                xstore[BASICLU_MEMORYU]: length of Ui and Ux

            basiclu_solve_for_update() will start from scratch in the next call.
*/
