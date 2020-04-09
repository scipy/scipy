/* 
 * lu_initialize.c
 *
 * Copyright (C) 2016-2018  ERGO-Code
 *
 */

#include "lu_internal.h"


/*
 * lu_initialize()
 *
 * Make @istore, @xstore a BASICLU instance. Set parameters to defaults and
 * initialize global counters. Reset instance for a fresh factorization.
 *
 */
void lu_initialize(lu_int m, lu_int *istore, double *xstore)
{
    struct lu this;

    /* set constant entries */
    istore[0]                               = BASICLU_HASH;
    xstore[0]                               = BASICLU_HASH;
    xstore[BASICLU_DIM]                     = m;

    /* set default parameters */
    xstore[BASICLU_MEMORYL]                 = 0;
    xstore[BASICLU_MEMORYU]                 = 0;
    xstore[BASICLU_MEMORYW]                 = 0;
    xstore[BASICLU_DROP_TOLERANCE]          = 1e-20;
    xstore[BASICLU_ABS_PIVOT_TOLERANCE]     = 1e-14;
    xstore[BASICLU_REL_PIVOT_TOLERANCE]     = 0.1;
    xstore[BASICLU_BIAS_NONZEROS]           = 1;
    xstore[BASICLU_MAXN_SEARCH_PIVOT]       = 3;
    xstore[BASICLU_PAD]                     = 4;
    xstore[BASICLU_STRETCH]                 = 0.3;
    xstore[BASICLU_COMPRESSION_THRESHOLD]   = 0.5;
    xstore[BASICLU_SPARSE_THRESHOLD]        = 0.05;
    xstore[BASICLU_REMOVE_COLUMNS]          = 0;
    xstore[BASICLU_SEARCH_ROWS]             = 1;

    /* initialize global counters */
    xstore[BASICLU_NFACTORIZE]              = 0;
    xstore[BASICLU_NUPDATE_TOTAL]           = 0;
    xstore[BASICLU_NFORREST_TOTAL]          = 0;
    xstore[BASICLU_NSYMPERM_TOTAL]          = 0;
    xstore[BASICLU_TIME_FACTORIZE_TOTAL]    = 0;
    xstore[BASICLU_TIME_SOLVE_TOTAL]        = 0;
    xstore[BASICLU_TIME_UPDATE_TOTAL]       = 0;

    /* lu_reset() and lu_save() initializes the remaining slots */
    lu_load(&this, istore, xstore, NULL, NULL, NULL, NULL, NULL, NULL);
    lu_reset(&this);
    lu_save(&this, istore, xstore, BASICLU_OK);
}
