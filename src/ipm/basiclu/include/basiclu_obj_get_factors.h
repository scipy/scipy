lu_int basiclu_obj_get_factors
(
    struct basiclu_object *obj,
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

    Call basiclu_get_factors() on a BASICLU object.

Return:

    BASICLU_ERROR_invalid_object

        obj is NULL or initialized to a null object.

    Other return codes are passed through from basiclu_get_factors().

Arguments:

    struct basiclu_object *obj

        Pointer to an initialized BASICLU object.

    The other arguments are passed through to basiclu_get_factors().
*/
