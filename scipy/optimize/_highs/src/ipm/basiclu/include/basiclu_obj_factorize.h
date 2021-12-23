lu_int basiclu_obj_factorize
(
    struct basiclu_object *obj,
    const lu_int *Bbegin,
    const lu_int *Bend,
    const lu_int *Bi,
    const double *Bx
);

/*
Purpose:

    Call basiclu_factorize() on a BASICLU object.

Return:

    BASICLU_ERROR_invalid_object

        obj is NULL or initialized to a null object.

    BASICLU_ERROR_out_of_memory

        reallocation failed because of insufficient memory.

    Other return codes are passed through from basiclu_factorize().

Arguments:

    struct basiclu_object *obj

        Pointer to an initialized BASICLU object.

    The other arguments are passed through to basiclu_factorize().
*/
