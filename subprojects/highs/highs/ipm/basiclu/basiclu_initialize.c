/*
 * basiclu_initialize.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 */

#include "ipm/basiclu/lu_internal.h"

lu_int basiclu_initialize
(
    lu_int m,
    lu_int istore[],
    double xstore[]
)
{
    if (!istore || !xstore)
        return BASICLU_ERROR_argument_missing;
    if (m <= 0)
        return BASICLU_ERROR_invalid_argument;

    lu_initialize(m, istore, xstore);
    return BASICLU_OK;
}
