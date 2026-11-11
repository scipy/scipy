lu_int basiclu_obj_initialize
(
    struct basiclu_object *obj,
    lu_int m
);

/*
Purpose:

    Initialize a BASICLU object. When m is positive, then *obj is initialized to
    process matrices of dimension m. When m is zero, then *obj is initialized to
    a "null" object, which cannot be used for factorization, but can be passed
    to basiclu_obj_free().

    This routine must be called once before passing obj to any other
    basiclu_obj_ routine. When obj is initialized to a null object, then the
    routine can be called again to reinitialize obj.

Return:

    BASICLU_OK

        *obj successfully initialized.

    BASICLU_ERROR_argument_missing

        obj is NULL.

    BASICLU_ERROR_invalid_argument

        m is negative.

    BASICLU_ERROR_out_of_memory

        insufficient memory to initialize object.

Arguments:

    struct basiclu_object *obj

        Pointer to the object to be initialized.

    lu_int m

        The dimension of matrices which can be processed, or 0.
*/
