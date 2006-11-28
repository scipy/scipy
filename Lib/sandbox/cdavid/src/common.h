/*
 * Last Change: Tue Nov 28 04:00 PM 2006 J
 *
 * Implements FP macros missing in C89
 */
#ifndef _GABSIG_C_COMMON_H
    #define _GABSIG_C_COMMON_H

#include <math.h>
#if defined(fpclassify)

    #if !defined(isnan)
        #define isnan(x) (fpclassify((x)) == FP_NAN)
    #endif
    #if !defined(isinf)
        #define isinf(x) (fpclassify((x)) == FP_INFINITE)
    #endif

#else  /* check to see if already have a function like this */
    #if !defined(HAVE_ISNAN)
        #if !defined(isnan)
            #define isnan(x) ((x) == (x))
        #endif
    #endif /* HAVE_ISNAN */
        
    #if !defined(HAVE_ISINF)
        #if !defined(isinf)
           #define isinf(x) (!isnan((x)) && isnan((x)-(x)))
        #endif
    #endif /* HAVE_ISINF */

#endif /* defined(fpclassify) */

#if !defined(isfinite)
    #define isfinite(x) (!isnan((x)) && !isinf((x)))
#endif

#endif /* nef of recursive header inclusion protection */
