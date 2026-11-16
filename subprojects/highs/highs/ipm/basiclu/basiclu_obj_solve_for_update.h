lu_int basiclu_obj_solve_for_update
(
    struct basiclu_object *obj,
    lu_int nzrhs,
    const lu_int irhs[],
    const double xrhs[],
    char trans,
    lu_int want_solution
);

/*
Purpose:

    Call basiclu_solve_for_update() on a BASICLU object. On success, if the
    solution was requested, it is provided in obj->lhs and the nonzero pattern
    is stored in obj->ilhs[0..obj->nzlhs-1].

Return:

    BASICLU_ERROR_invalid_object

        obj is NULL or initialized to a null object.

    BASICLU_ERROR_out_of_memory

        reallocation failed because of insufficient memory.

    Other return codes are passed through from basiclu_solve_for_update().

Arguments:

    struct basiclu_object *obj

        Pointer to an initialized BASICLU object.

    lu_int want_solution

        Nonzero to compute the solution to the linear system,
        zero to only prepare the update.

    The other arguments are passed through to basiclu_solve_for_update().
*/
