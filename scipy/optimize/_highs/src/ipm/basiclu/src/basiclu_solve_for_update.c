/*
 * basiclu_solve_for_update.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 */

#include "lu_internal.h"

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
)
{
    struct lu this;
    lu_int status, n, ok;

    status = lu_load(&this, istore, xstore, Li, Lx, Ui, Ux, Wi, Wx);
    if (status != BASICLU_OK)
        return status;

    if (! (Li && Lx && Ui && Ux && Wi && Wx && irhs))
    {
        status = BASICLU_ERROR_argument_missing;
    }
    else if (trans != 't' && trans != 'T' && !xrhs)
    {
        status = BASICLU_ERROR_argument_missing;
    }
    else if (this.nupdate < 0)
    {
        status = BASICLU_ERROR_invalid_call;
    }
    else if (this.nforrest == this.m)
    {
        status = BASICLU_ERROR_maximum_updates;
    }
    else
    {
        /* check RHS indices */
        if (trans == 't' || trans == 'T')
        {
            ok = irhs[0] >= 0 && irhs[0] < this.m;
        }
        else
        {
            ok = nzrhs >= 0 && nzrhs <= this.m;
            for (n = 0; n < nzrhs && ok; n++)
            {
                ok = ok && irhs[n] >= 0 && irhs[n] < this.m;
            }
        }
        if (!ok)
            status = BASICLU_ERROR_invalid_argument;
    }

    if (status == BASICLU_OK)
    {
        /* may request reallocation */
        status = lu_solve_for_update(&this, nzrhs, irhs, xrhs, p_nzlhs, ilhs,
                                     lhs, trans);
    }

    return lu_save(&this, istore, xstore, status);
}
