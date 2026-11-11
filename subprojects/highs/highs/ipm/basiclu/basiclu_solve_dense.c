/*
 * basiclu_solve_dense.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 */

#include "ipm/basiclu/lu_internal.h"

lu_int basiclu_solve_dense
(
    lu_int istore[],
    double xstore[],
    lu_int Li[],
    double Lx[],
    lu_int Ui[],
    double Ux[],
    lu_int Wi[],
    double Wx[],
    const double rhs[],
    double lhs[],
    char trans
)
{
    struct lu this;
    lu_int status;

    status = lu_load(&this, istore, xstore, Li, Lx, Ui, Ux, Wi, Wx);
    if (status != BASICLU_OK)
        return status;

    if (! (Li && Lx && Ui && Ux && Wi && Wx && rhs && lhs))
    {
        status = BASICLU_ERROR_argument_missing;
    }
    else if (this.nupdate < 0)
    {
        status = BASICLU_ERROR_invalid_call;
    }
    else
    {
        lu_solve_dense(&this, rhs, lhs, trans);
    }

    return lu_save(&this, istore, xstore, status);
}
