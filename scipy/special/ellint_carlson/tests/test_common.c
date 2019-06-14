#include "test_common.h"


void testarray_init(testarray_t * array, size_t ntests, size_t nargs,
                    const EllInt_Num_t * restrict results,
		    const EllInt_Num_t * * restrict args, ellint_testfcn_t f)
{
    array->rerr = ELLINT_TEST_RERR;
    array->ntests = ntests;
    array->nargs = nargs;
    array->res_array = results;
    array->args_array = args;
    array->wrapper = f;
}

int testarray_exec(const testarray_t * array)
{
    int npass, nfail;
    size_t i, j;

    npass = 0;
    nfail = 0;
    for ( i = 0; i < array->ntests; i++ )
    {
	int status;
	EllInt_Num_t out;
	EllInt_Num_t expected;
	EllInt_Num_t * argrow;

	argrow = (EllInt_Num_t *)(array->args_array) + i * (array->nargs);
	expected = (array->res_array)[i];
	for ( j = 0; j < array->nargs; j++ )
	{
	    printf("argument %lu: % .4f%+.4fI\n",
		   j, creal(argrow[j]), cimag(argrow[j]));
	}
	printf("expected  : % #0.13a, % #0.13a\n",
               creal(expected), cimag(expected));

	status = (array->wrapper)(argrow, array->rerr, &out);

	if ( status == ELLINT_STATUS_SUCCESS )
	{
	    double diff;
	    bool passed;

	    printf("actual    : % #0.13a, % #0.13a\n", creal(out), cimag(out));
	    if ( fpclassify(fabs(expected)) == FP_ZERO )
	    {
		diff = fabs(out - expected);
	    } else {
		diff = fabs(out - expected) / fabs(expected);
	    }
	    passed = ( diff <= array->rerr );
	    printf("diff      : % .3e %s %.3e (%s)\n",
		    diff, passed ? "<=" : ">", array->rerr,
		    passed ? "PASS" : "FAIL");
	    if ( passed )
	    {
		npass++;
	    } else {
		nfail++;
	    }
	} else {
	    nfail++;
	    printf("actual    : error status %d\n", status);
	}
	printf("\n");
    }
    printf("Tested:% 3d; Pass:% 3d; Fail:% 3d\n",
           (int)(array->ntests), npass, nfail);
    printf("------------------------------\n\n");
    return nfail;
}
