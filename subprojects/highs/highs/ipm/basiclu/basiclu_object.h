struct basiclu_object
{
    lu_int *istore;
    double *xstore;
    lu_int *Li, *Ui, *Wi;
    double *Lx, *Ux, *Wx;
    double *lhs;
    lu_int *ilhs;
    lu_int nzlhs;
    double realloc_factor;
};

/*
A variable of type struct basiclu_object must be defined in user code. Its
members are set and maintained by basiclu_obj_* routines. User code should only
access the following members:

    xstore (read/write)

        set parameters and get info values

    lhs, ilhs, nzlhs (read only)

        holds solution after solve_sparse() and solve_for_update()

    realloc_factor (read/write)

        Arrays are reallocated for max(realloc_factor, 1.0) times the
        required size. Default: 1.5
*/
