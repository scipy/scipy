#include "scomplex.h"

double scsum1_(int *n, complex *cx, int *incx)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    SCSUM1 takes the sum of the absolute values of a complex   
    vector and returns a single precision result.   

    Based on SCASUM from the Level 1 BLAS.   
    The change is to use the 'genuine' absolute value.   

    Contributed by Nick Higham for use with CLACON.   

    Arguments   
    =========   

    N       (input) INT
            The number of elements in the vector CX.   

    CX      (input) COMPLEX array, dimension (N)   
            The vector whose elements will be summed.   

    INCX    (input) INT
            The spacing between successive values of CX.  INCX > 0.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int i__1, i__2;
    float ret_val;
    /* Builtin functions */
    double c_abs(complex *);
    /* Local variables */
    static int i, nincx;
    static float stemp;


#define CX(I) cx[(I)-1]


    ret_val = 0.f;
    stemp = 0.f;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*     CODE FOR INCREMENT NOT EQUAL TO 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i = 1; *incx < 0 ? i >= nincx : i <= nincx; i += *incx) {

/*        NEXT LINE MODIFIED. */

	stemp += c_abs(&CX(i));
/* L10: */
    }
    ret_val = stemp;
    return ret_val;

/*     CODE FOR INCREMENT EQUAL TO 1 */

L20:
    i__2 = *n;
    for (i = 1; i <= *n; ++i) {

/*        NEXT LINE MODIFIED. */

	stemp += c_abs(&CX(i));
/* L30: */
    }
    ret_val = stemp;
    return ret_val;

/*     End of SCSUM1 */

} /* scsum1_ */

