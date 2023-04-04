/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/
/*! @file icmax1.c
 * \brief Finds the index of the element whose real part has maximum absolute value
 *
 * <pre>
 *     -- LAPACK auxiliary routine (version 2.0) --   
 *     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
 *     Courant Institute, Argonne National Lab, and Rice University   
 *     October 31, 1992   
 * </pre>
 */
#include <math.h>
#include "slu_scomplex.h"
#include "slu_Cnames.h"

/*! \brief

 <pre>
    Purpose   
    =======   

    ICMAX1 finds the index of the element whose real part has maximum   
    absolute value.   

    Based on ICAMAX from Level 1 BLAS.   
    The change is to use the 'genuine' absolute value.   

    Contributed by Nick Higham for use with CLACON.   

    Arguments   
    =========   

    N       (input) INT   
            The number of elements in the vector CX.   

    CX      (input) COMPLEX array, dimension (N)   
            The vector whose elements will be summed.   

    INCX    (input) INT   
            The spacing between successive values of CX.  INCX >= 1.   

   ===================================================================== 
  </pre>
*/
int icmax1_slu(int *n, complex *cx, int *incx)
{
/*
       NEXT LINE IS THE ONLY MODIFICATION.   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int ret_val, i__1, i__2;
    float r__1;
    /* Local variables */
    float smax;
    int i, ix;


#define CX(I) cx[(I)-1]


    ret_val = 0;
    if (*n < 1) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L30;
    }

/*     CODE FOR INCREMENT NOT EQUAL TO 1 */

    ix = 1;
    smax = (r__1 = CX(1).r, fabs(r__1));
    ix += *incx;
    i__1 = *n;
    for (i = 2; i <= *n; ++i) {
	i__2 = ix;
	if ((r__1 = CX(ix).r, fabs(r__1)) <= smax) {
	    goto L10;
	}
	ret_val = i;
	i__2 = ix;
	smax = (r__1 = CX(ix).r, fabs(r__1));
L10:
	ix += *incx;
/* L20: */
    }
    return ret_val;

/*     CODE FOR INCREMENT EQUAL TO 1 */

L30:
    smax = (r__1 = CX(1).r, fabs(r__1));
    i__1 = *n;
    for (i = 2; i <= *n; ++i) {
	i__2 = i;
	if ((r__1 = CX(i).r, fabs(r__1)) <= smax) {
	    goto L40;
	}
	ret_val = i;
	i__2 = i;
	smax = (r__1 = CX(i).r, fabs(r__1));
L40:
	;
    }
    return ret_val;

/*     End of ICMAX1 */

} /* icmax1_slu */

