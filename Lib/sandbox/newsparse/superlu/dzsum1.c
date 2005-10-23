#include "dcomplex.h"

double dzsum1_(int *n, doublecomplex *cx, int *incx)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   

    Purpose   
    =======   

    DZSUM1 takes the sum of the absolute values of a complex   
    vector and returns a double precision result.   

    Based on DZASUM from the Level 1 BLAS.   
    The change is to use the 'genuine' absolute value.   

    Contributed by Nick Higham for use with ZLACON.   

    Arguments   
    =========   

    N       (input) INT   
            The number of elements in the vector CX.   

    CX      (input) COMPLEX*16 array, dimension (N)   
            The vector whose elements will be summed.   

    INCX    (input) INT   
            The spacing between successive values of CX.  INCX > 0.   

    ===================================================================== 
*/  

    /* Builtin functions */
    double z_abs(doublecomplex *);
    
    /* Local variables */
    int i, nincx;
    double stemp;


#define CX(I) cx[(I)-1]

    stemp = 0.;
    if (*n <= 0) {
	return stemp;
    }
    if (*incx == 1) {
	goto L20;
    }

    /*     CODE FOR INCREMENT NOT EQUAL TO 1 */

    nincx = *n * *incx;
    for (i = 1; *incx < 0 ? i >= nincx : i <= nincx; i += *incx) {

	/*        NEXT LINE MODIFIED. */

	stemp += z_abs(&CX(i));
/* L10: */
    }
    
    return stemp;

    /*     CODE FOR INCREMENT EQUAL TO 1 */

L20:
    for (i = 1; i <= *n; ++i) {

	/*        NEXT LINE MODIFIED. */

	stemp += z_abs(&CX(i));
/* L30: */
    }
    
    return stemp;

    /*     End of DZSUM1 */

} /* dzsum1_ */

