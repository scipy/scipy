#ifndef _ELLINT_POLY_H_INCLUDED
#define _ELLINT_POLY_H_INCLUDED

#if ( (defined ELLINT_POLY_REAL) && (defined ELLINT_POLY_COMPLEX) )
    #error "either ELLINT_POLY_REAL or ELLINT_POLY_COMPLEX must be defined"
#elif ( (!(defined ELLINT_POLY_REAL)) && (!(defined ELLINT_POLY_COMPLEX)) )
    #error "either ELLINT_POLY_REAL or ELLINT_POLY_COMPLEX must be defined"
#endif

#ifdef ELLINT_POLY_REAL
    #define ELLINT_POLY_FCN(fcn)	(f##fcn)
    #define ELLINT_RNEG(x)		( (x) < 0.0 )
#else
    #define ELLINT_POLY_FCN(fcn)	(c##fcn)
    #define ELLINT_RNEG(x)		( creal(x) < 0.0 )
#endif

#endif  /* _ELLINT_POLY_H_INCLUDED */
