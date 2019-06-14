#ifndef _ELLINT_TEST_COMMON_H_INCLUDED
#define _ELLINT_TEST_COMMON_H_INCLUDED
#include <stddef.h>
#include <stdlib.h>
#include <tgmath.h>
#include <stdbool.h>
#include <stdio.h>
#include "ellint_carlson.h"
#include "ellint_common.h"
#include "ellint_poly.h"


#ifndef ELLINT_TEST_RERR
/* The tests are based on the 14-digit tabulated values. Comparison beyond the
 * last significant digit therefore cannot be guaranteed to be meaningful.
 */
#define ELLINT_TEST_RERR	(5e-14)
#endif

typedef int (*ellint_testfcn_t)(const EllInt_Num_t [], double,
                                EllInt_Num_t * restrict);

typedef struct testarray
{
    size_t ntests;
    size_t nargs;
    double rerr;
    const EllInt_Num_t * res_array;
    const EllInt_Num_t * * args_array;
    ellint_testfcn_t wrapper;
} testarray_t;


extern void testarray_init(testarray_t * array, size_t ntests, size_t nargs,
                           const EllInt_Num_t * restrict results,
	  	           const EllInt_Num_t * * restrict args,
			   ellint_testfcn_t f);

extern int testarray_exec(const testarray_t * array);

#define TESTARRAY_MAIN(VOID)	\
do {				\
    size_t ntests;		\
    int status;			\
    testarray_t thistest;	\
    ntests = (sizeof eres) / (sizeof (EllInt_Num_t));	\
    testarray_init(&thistest, ntests, NARGS, eres,	\
                   (const EllInt_Num_t * *)args, testwrapper);	\
    status = testarray_exec(&thistest);			\
    return status;		\
} while ( 0 )
#endif  /* _ELLINT_TEST_COMMON_H_INCLUDED */
