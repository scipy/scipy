#ifndef _ELLINT_COMPAT_H_INCLUDED
#define _ELLINT_COMPAT_H_INCLUDED


#include <math.h>
#include <complex.h>


#if ( __STDC_VERSION__ >= 199901L )
    #define ELLINT_HAS_C99	(1)
#else
    #define ELLINT_HAS_C99	(0)
#endif


#if ( ELLINT_HAS_C99 )
    #include <stdbool.h>
#else
    #define bool		unsigned char
    #define true		( (bool)1 )
    #define false		( (bool)0 )
#endif


#ifndef _MSC_VER	/* Not MSVC, fairly straightforward C99/C11 */


typedef double complex	double_complex;


#if ( !(ELLINT_HAS_C99) )
    #ifdef __GNUC__
	#define restrict	__restrict__
	#define inline		__inline__
    #else
	#define restrict	/* nothing */
	#define inline		/* nothing */
    #endif
    #define complex		_Complex
#endif


#ifdef ELLINT_POLY_REAL
    #define EllInt_Num_t	double
    #define ELLINT_NAN		( nan("") )
    #define MKCMPLX(r, IVOID)	( (r) )
    #define CREAL(z)		( (double)(z) )
    #define CIMAG(z)		( (double)(0.0) )
    #define CONJ(z)		( (double)(z) )
    
    #define FABS(z)		( fabs( (z) ) )
    #define SQRT(z)		( sqrt( (z) ) )

    #define CARG(z)		( ( (z) < 0.0 ) ? M_PI : 0.0 )
#else				/* complex */
    #define EllInt_Num_t	double_complex
    /* CMPLX macro is a C99 feature, circumvent it to conform to C89 if that is
     * the C standard being used */
    #ifndef CMPLX
	#if ( defined(__GNUC__) &&	\
	      ( __GNUC__ * 10000 + __GNUC_MINOR__ * 100 >= 40700 ) )
	    #define MKCMPLX(r, i)	__builtin_complex( (r), (i) )
	#else	/* fallback, will choke on non-regular floats such as nan */
	    #define MKCMPLX(r, i)	\
		( (double_complex)((double)(r) + (double)(i) * I) )
	#endif
    #else
	#define MKCMPLX(r, i)		( CMPLX( (r), (i) ) )
    #endif
    #define ELLINT_NAN			( MKCMPLX(nan(""), nan("")) )

    #define CREAL(z)		( creal(z) )
    #define CIMAG(z)		( cimag(z) )
    #define CONJ(z)		( conj(z) )
    
    #define FABS(z)		( cabs( (z) ) )
    #define SQRT(z)		( csqrt( (z) ) )
    #define CARG(z)		( carg( (z) ) )
#endif

#define NEG(x)	( -(x) )

#define ADD(x, y)	( (x) + (y) )
#define ADDcr(x, y)	( (x) + (y) )

#define SUB(x, y)	( (x) - (y) )

#define MULcr(x, y)	( (x) * (y) )
#define MULcc(x, y)	( (x) * (y) )

#define RECI(z)		( 1.0 / (double)(z) )

#define DIVcr(x, y)	( (x) / (y) )
#define DIVcc(x, y)	( (x) / (y) )
#define DIVrc(x, y)	( (x) / (y) )


#else			/* MSVC, type incompatibility, should do special stuff
			   for both real and complex */


#if ( _MSC_VER < 1900 )
#warning("MSVC version may not support required language features.")
#endif


#define	restrict	__restrict
#define inline		__inline


typedef _Dcomplex double_complex;


#ifdef ELLINT_POLY_REAL		/* real */
    #define EllInt_Num_t	double
    #define ELLINT_NAN		( nan("") )
    #define MKCMPLX(r, IVOID)	( (r) )

    #define NEG(x)	( -(x) )

    #define ADD(x, y)	( (x) + (y) )
    #define ADDcr(x, y)	( (x) + (y) )

    #define SUB(x, y)	( (x) - (y) )

    #define MULcc(x, y)	( (x) * (y) )
    #define MULcr(x, y)	( (x) * (y) )

    #define DIVcc(x, y)	( (x) / (y) )
    #define DIVcr(x, y)	( (x) / (double)(y) )
    #define DIVrc(x, y)	( (double)(x) / (y) )

    #define CREAL(z)	( (double)(z) )
    #define CIMAG(z)	( (double)(0.0) )
    #define CONJ(z)	( (double)(z) )

    #define FABS(z)	( fabs( (z) ) )
    #define SQRT(z)	( sqrt( (z) ) )

    #define CARG(z)	( ( (z) < 0.0 ) ? M_PI : 0.0 )
#else				/* complex */
    #define EllInt_Num_t	double_complex
    #define MKCMPLX(r, i)	( _Cbuild( (double)(r) , (double)(i) ) )
    #define ELLINT_NAN		( MKCMPLX(nan(""), nan("")) )

    #define CREAL(z)	( creal( (z) ) )
    #define CIMAG(z)	( cimag( (z) ) )
    #define CONJ(z)	( conj( (z) ) )

    #define NEG(x)	( MKCMPLX( -CREAL( (x) ), -CIMAG( (x) ) ) )

    #define ADD(x, y)	\
	( MKCMPLX( creal( (x) ) + creal( (y) ), cimag( (x) ) + cimag( (y) ) ) )
    #define ADDcr(x, y)	( ADD( (x), MKCMPLX( (y), 0.0) ) )

    #define SUB(x, y)	\
	( MKCMPLX( creal( (x) ) - creal( (y) ), cimag( (x) ) - cimag( (y) ) ) )

    #define MULcc(x, y)	( _Cmulcc( (x), (y) ) )
    #define MULcr(x, y)	( _Cmulcr( (x), (double)(y) ) )

    #define RECI(z)	( conj( (z) ) / norm( (z) ) )
    #define DIVcc(x, y)	( _Cmulcc( (x), RECI( (y) ) ) )
    #define DIVcr(x, y)	( _Cmulcr( (x), (1.0 / (double)(y)) ) )
    #define DIVrc(x, y)	( MULcr(RECI( (y) ), x) )

    #define FABS(z)	( cabs( (z) ) )
    #define SQRT(z)	( csqrt( (z) ) )

    #define CARG(z)	( carg( (z) ) )
#endif


#endif			/* Selective compilation, MSVC vs. C99/C11 compiler */


#define MKCR(r)		( MKCMPLX( (r), 0.0 ) )


#endif /* _ELLINT_COMPAT_H_INCLUDED */
