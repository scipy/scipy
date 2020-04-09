/*
 * basiclu_update.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 */

#include "lu_internal.h"

lu_int basiclu_update
(
    lu_int istore[],
    double xstore[],
    lu_int Li[],
    double Lx[],
    lu_int Ui[],
    double Ux[],
    lu_int Wi[],
    double Wx[],
    double xtbl
)
{
    struct lu this;
    lu_int status;

    status = lu_load(&this, istore, xstore, Li, Lx, Ui, Ux, Wi, Wx);
    if (status != BASICLU_OK)
        return status;

    if (! (Li && Lx && Ui && Ux && Wi && Wx))
    {
        status = BASICLU_ERROR_argument_missing;
    }
    else if (this.nupdate < 0 || this.ftran_for_update < 0 ||
             this.btran_for_update < 0)
    {
        status = BASICLU_ERROR_invalid_call;
    }
    else
    {
        status = lu_update(&this, xtbl);
    }
    return lu_save(&this, istore, xstore, status);
}
