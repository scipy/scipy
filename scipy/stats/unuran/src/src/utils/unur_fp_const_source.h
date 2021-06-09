/** File automatically created by scripts/compute_machine_constants         **/

/*****************************************************************************
 *                                                                           *
 *          unuran -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifndef UNUR_FP_CONST_SOURCE_H_SEEN
#define UNUR_FP_CONST_SOURCE_H_SEEN
/*---------------------------------------------------------------------------*/

/* maximal relative error when testing equality of two doubles */
#define UNUR_EPSILON  2.22044604925031308084726333618e-14

/* square root of DBL_EPSILON */
#define UNUR_SQRT_DBL_EPSILON   1.490116119384765625e-08

/* log of DBL_EPSILON */
#define UNUR_LOG_DBL_EPSILON   -36.0436533891171535515240975656

/* the machine roundoff error */
#define MACHEP  1.11022302462515654042363166809e-16

/* largest argument for exp() */
#define MAXLOG  709.782712893383973096206318587

/* smallest argument for exp() without underflow */
#define MINLOG  -708.396418532264078748994506896

/* the maximal number that pow(x,x-0.5) has no overflow */
/* we use a (very) conservative portable bound          */
#define MAXSTIR  108.116855767857671821730036754

/*---------------------------------------------------------------------------*/
#endif  /* UNUR_FP_CONST_SOURCE_H_SEEN */
/*---------------------------------------------------------------------------*/
