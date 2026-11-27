#ifndef _LU_DEF_H
#define _LU_DEF_H

/* -------------------------------------------------------------------------- */
/* ANSI standard include files */
/* -------------------------------------------------------------------------- */

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>

#include "ipm/basiclu/basiclu.h"

#define BASICLU_HASH 7743090    /* hash in istore[0], xstore[0] */

enum { NO_TASK, SINGLETONS, SETUP_BUMP, FACTORIZE_BUMP, BUILD_FACTORS };

/* -------------------------------------------------------------------------- */
/* standard macros and inlines */
/* -------------------------------------------------------------------------- */

#define MAX(a,b) ((a)>=(b) ? (a):(b))
#define MIN(a,b) ((a)<=(b) ? (a):(b))

static inline void lu_iswap(lu_int *x, lu_int i, lu_int j)
{
    lu_int t = x[i]; x[i] = x[j]; x[j] = t;
}

static inline void lu_fswap(double *x, lu_int i, lu_int j)
{
    double t = x[i]; x[i] = x[j]; x[j] = t;
}

#endif
