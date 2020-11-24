void basiclu_obj_free
(
    struct basiclu_object *obj
);

/*
Purpose:

    Free memory allocated from a BASICLU object. The object must have been
    initialized before by basiclu_obj_initialize(). Subsequent calls to
    basiclu_obj_free() will do nothing.

Arguments:

    struct basiclu_object *obj

        Pointer to the object which memory is to be freed. When obj is NULL,
        then the routine does nothing.
*/
