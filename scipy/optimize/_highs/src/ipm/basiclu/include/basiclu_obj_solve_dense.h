lu_int basiclu_obj_solve_dense
(
    struct basiclu_object *obj,
    const double rhs[],
    double lhs[],
    char trans
);

/*
Purpose:

    Call basiclu_solve_dense() on a BASICLU object.

Return:

    BASICLU_ERROR_invalid_object

        obj is NULL or initialized to a null object.

    Other return codes are passed through from basiclu_solve_dense().

Arguments:

    struct basiclu_object *obj

        Pointer to an initialized BASICLU object.

    The other arguments are passed through to basiclu_solve_dense().
*/
