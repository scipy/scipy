/* ========================================================================== */
/* === umf_version.h ======================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.1 (Apr. 30, 2003), Copyright (c) 2003 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
   Define routine names, depending on version being compiled.

   DINT:	double precision, int's as integers
   DLONG:	double precision, long's as integers
   ZLONG:	complex double precision, long's as integers
   ZINT:	complex double precision, int's as integers
*/

/* Set DINT as the default, if nothing is defined */
#if !defined (DLONG) && !defined (DINT) && !defined (ZLONG) && !defined (ZINT)
#define DINT
#endif

/* Determine if this is a real or complex version */
#if defined (ZLONG) || defined (ZINT)
#define COMPLEX
#endif

/* -------------------------------------------------------------------------- */
/* integer type (Int is int or long) now defined in amd_internal.h */
/* -------------------------------------------------------------------------- */

#if defined (DLONG) || defined (ZLONG)
#define LONG_INTEGER
#endif

/* -------------------------------------------------------------------------- */
/* Numerical relop macros for correctly handling the NaN case */
/* -------------------------------------------------------------------------- */

/*
SCALAR_IS_NAN(x):
    True if x is NaN.  False otherwise.  The commonly-existing isnan(x)
    function could be used, but it's not in Kernighan & Ritchie 2nd edition
    (ANSI C).  It may appear in <math.h>, but I'm not certain about
    portability.  The expression x != x is true if and only if x is NaN,
    according to the IEEE 754 floating-point standard.

SCALAR_IS_ZERO(x):
    True if x is zero.  False if x is nonzero, NaN, or +/- Inf.
    This is (x == 0) if the compiler is IEEE 754 compliant.

SCALAR_IS_NONZERO(x):
    True if x is nonzero, NaN, or +/- Inf.  False if x zero.
    This is (x != 0) if the compiler is IEEE 754 compliant.

SCALAR_IS_LTZERO(x):
    True if x is < zero or -Inf.  False if x is >= 0, NaN, or +Inf.
    This is (x < 0) if the compiler is IEEE 754 compliant.
*/

#if defined (MATHWORKS)

/* The MathWorks has their own macros in util.h that handle NaN's properly. */
#define SCALAR_IS_NAN(x)	(utIsNaN (x))
#define SCALAR_IS_ZERO(x)	(utEQZero (x))
#define SCALAR_IS_NONZERO(x)	(utNEZero (x))
#define SCALAR_IS_LTZERO(x)	(utLTZero (x))

#elif defined (UMF_WINDOWS)

/* Yes, this is exceedingly ugly.  Blame Microsoft, which hopelessly */
/* violates the IEEE 754 floating-point standard in a bizarre way. */
/* If you're using an IEEE 754-compliant compiler, then x != x is true */
/* iff x is NaN.  For Microsoft, (x < x) is true iff x is NaN. */
/* So either way, this macro safely detects a NaN. */
#define SCALAR_IS_NAN(x)	(((x) != (x)) || (((x) < (x))))
#define SCALAR_IS_ZERO(x)	(((x) == 0.) && !SCALAR_IS_NAN(x))
#define SCALAR_IS_NONZERO(x)	(((x) != 0.) || SCALAR_IS_NAN(x))
#define SCALAR_IS_LTZERO(x)	(((x) < 0.) && !SCALAR_IS_NAN(x))

#else

/* These all work properly, according to the IEEE 754 standard ... except on */
/* a PC with windows.  Works fine in Linux on the same PC... */
#define SCALAR_IS_NAN(x)	((x) != (x))
#define SCALAR_IS_ZERO(x)	((x) == 0.)
#define SCALAR_IS_NONZERO(x)	((x) != 0.)
#define SCALAR_IS_LTZERO(x)	((x) < 0.)

#endif

/* scalar absolute value macro. If x is NaN, the result is NaN: */
#define SCALAR_ABS(x) ((SCALAR_IS_LTZERO (x)) ? -(x) : (x))

/* true if an integer (stored in double x) would overflow (or if x is NaN) */
#define INT_OVERFLOW(x) ((!((x) * (1.0+1e-8) <= (double) Int_MAX)) \
			|| SCALAR_IS_NAN (x))

/* print a scalar (avoid printing "-0" for negative zero).  */
#define PRINT_SCALAR(a) \
{ \
    if (SCALAR_IS_NONZERO (a)) \
    { \
	PRINTF ((" (%g)", (a))) ; \
    } \
    else \
    { \
	PRINTF ((" (0)")) ; \
    } \
}

/* -------------------------------------------------------------------------- */
/* Real floating-point arithmetic */
/* -------------------------------------------------------------------------- */

#ifndef COMPLEX

#define Entry double

#define REAL_COMPONENT(c)		(c)
#define IMAG_COMPONENT(c)		(0.)
#define ASSIGN(c,s1,s2)		    { (c) = (s1) ; }
#define CLEAR(c)		    { (c) = 0. ; }
#define CLEAR_AND_INCREMENT(p)	    { *p++ = 0. ; }
#define IS_NAN(a)		    SCALAR_IS_NAN (a)
#define IS_ZERO(a)		    SCALAR_IS_ZERO (a)
#define IS_NONZERO(a)		    SCALAR_IS_NONZERO (a)
#define SCALE_DIV(c,s)		    { (c) /= (s) ; }
#ifndef NRECIPROCAL
#define SCALE_RECIP(c,s)	    { (c) *= (s) ; }
#endif
#define ASSEMBLE(c,a)		    { (c) += (a) ; }
#define ASSEMBLE_AND_INCREMENT(c,p) { (c) += *p++ ; }
#define DECREMENT(c,a)		    { (c) -= (a) ; }
#define MULT(c,a,b)		    { (c) = (a) * (b) ; }
#define MULT_CONJ(c,a,b)	    { (c) = (a) * (b) ; }
#define MULT_SUB(c,a,b)		    { (c) -= (a) * (b) ; }
#define MULT_SUB_CONJ(c,a,b)	    { (c) -= (a) * (b) ; }
#define DIV(c,a,b)		    { (c) = (a) / (b) ; }
#define RECIPROCAL(c)		    { (c) = 1.0 / (c) ; }
#define DIV_CONJ(c,a,b)		    { (c) = (a) / (b) ; }
#define APPROX_ABS(s,a)		    { (s) = SCALAR_ABS (a) ; }
#define ABS(s,a)		    { (s) = SCALAR_ABS (a) ; }
#define PRINT_ENTRY(a)		    PRINT_SCALAR (a)

/* for flop counts */
#define MULTSUB_FLOPS	2.	/* c -= a*b */
#define DIV_FLOPS	1.	/* c = a/b */
#define ABS_FLOPS	0.	/* c = abs (a) */
#define ASSEMBLE_FLOPS	1.	/* c += a */
#define DECREMENT_FLOPS	1.	/* c -= a */
#define MULT_FLOPS	1.	/* c = a*b */
#define SCALE_FLOPS	1.	/* c = a/s */

#else

/* -------------------------------------------------------------------------- */
/* Complex floating-point arithmetic */
/* -------------------------------------------------------------------------- */

/*
    Note:  An alternative to this DoubleComplex type would be to use a
    struct { double r ; double i ; }.  The problem with that method
    (used by the Sun Performance Library, for example) is that ANSI C provides
    no guarantee about the layout of a struct.  It is possible that the sizeof
    the struct above would be greater than 2 * sizeof (double).  This would
    mean that the complex BLAS could not be used.  The method used here avoids
    that possibility.  ANSI C *does* guarantee that an array of structs has
    the same size as n times the size of one struct.

    The ANSI C99 version of the C language includes a "double _Complex" type.
    It should be possible in that case to do the following:

    #define Entry double _Complex

    and remove the DoubleComplex struct.  The macros, below, could then be
    replaced with instrinsic operators.  Note that the #define Real and
    #define Imag should also be removed (they only appear in this file).

    For the MULT, MULT_SUB, MULT_SUB_CONJ, and MULT_CONJ macros,
    the output argument c cannot be the same as any input argument.

*/

typedef struct
{
    double component [2] ;	/* real and imaginary parts */

} DoubleComplex ;

#define Entry DoubleComplex
#define Real component [0]
#define Imag component [1]

/* for flop counts */
#define MULTSUB_FLOPS	8.	/* c -= a*b */
#define DIV_FLOPS	9.	/* c = a/b */
#define ABS_FLOPS	6.	/* c = abs (a), count sqrt as one flop */
#define ASSEMBLE_FLOPS	2.	/* c += a */
#define DECREMENT_FLOPS	2.	/* c -= a */
#define MULT_FLOPS	6.	/* c = a*b */
#define SCALE_FLOPS	2.	/* c = a/s or c = a*s */

/* -------------------------------------------------------------------------- */

/* real part of c */
#define REAL_COMPONENT(c) ((c).Real)

/* -------------------------------------------------------------------------- */

/* imag part of c */
#define IMAG_COMPONENT(c) ((c).Imag)

/* -------------------------------------------------------------------------- */

/* c = (s1) + (s2)i */
#define ASSIGN(c,s1,s2) \
{ \
    (c).Real = (s1) ; \
    (c).Imag = (s2) ; \
}

/* -------------------------------------------------------------------------- */

/* c = 0 */
#define CLEAR(c) \
{ \
    (c).Real = 0. ; \
    (c).Imag = 0. ; \
}

/* -------------------------------------------------------------------------- */

/* *p++ = 0 */
#define CLEAR_AND_INCREMENT(p) \
{ \
    p->Real = 0. ; \
    p->Imag = 0. ; \
    p++ ; \
}

/* -------------------------------------------------------------------------- */

/* True if a == 0 */
#define IS_ZERO(a) \
    (SCALAR_IS_ZERO ((a).Real) && SCALAR_IS_ZERO ((a).Imag))

/* -------------------------------------------------------------------------- */

/* True if a is NaN */
#define IS_NAN(a) \
    (SCALAR_IS_NAN ((a).Real) || SCALAR_IS_NAN ((a).Imag))

/* -------------------------------------------------------------------------- */

/* True if a != 0 */
#define IS_NONZERO(a) \
    (SCALAR_IS_NONZERO ((a).Real) || SCALAR_IS_NONZERO ((a).Imag))

/* -------------------------------------------------------------------------- */

/* c /= s */
#define SCALE_DIV(c,s) \
{ \
    (c).Real /= (s) ; \
    (c).Imag /= (s) ; \
}

/* -------------------------------------------------------------------------- */

/* c *= s, where s is the reciprocal scale factor.  Not used if
 * NRECIPROCAL is defined at compile time. */
#ifndef NRECIPROCAL
#define SCALE_RECIP(c,s) \
{ \
    (c).Real *= (s) ; \
    (c).Imag *= (s) ; \
}
#endif

/* -------------------------------------------------------------------------- */

/* c += a */
#define ASSEMBLE(c,a) \
{ \
    (c).Real += (a).Real ; \
    (c).Imag += (a).Imag ; \
}

/* -------------------------------------------------------------------------- */

/* c += *p++ */
#define ASSEMBLE_AND_INCREMENT(c,p) \
{ \
    (c).Real += p->Real ; \
    (c).Imag += p->Imag ; \
    p++ ; \
}

/* -------------------------------------------------------------------------- */

/* c -= a */
#define DECREMENT(c,a) \
{ \
    (c).Real -= (a).Real ; \
    (c).Imag -= (a).Imag ; \
}

/* -------------------------------------------------------------------------- */

/* c = a*b, assert because c cannot be the same as a or b */
#define MULT(c,a,b) \
{ \
    ASSERT (&(c) != &(a) && &(c) != &(b)) ; \
    (c).Real = (a).Real * (b).Real - (a).Imag * (b).Imag ; \
    (c).Imag = (a).Imag * (b).Real + (a).Real * (b).Imag ; \
}

/* -------------------------------------------------------------------------- */

/* c = a*conjugate(b), assert because c cannot be the same as a or b */
#define MULT_CONJ(c,a,b) \
{ \
    ASSERT (&(c) != &(a) && &(c) != &(b)) ; \
    (c).Real = (a).Real * (b).Real + (a).Imag * (b).Imag ; \
    (c).Imag = (a).Imag * (b).Real - (a).Real * (b).Imag ; \
}

/* -------------------------------------------------------------------------- */

/* c -= a*b, assert because c cannot be the same as a or b */
#define MULT_SUB(c,a,b) \
{ \
    ASSERT (&(c) != &(a) && &(c) != &(b)) ; \
    (c).Real -= (a).Real * (b).Real - (a).Imag * (b).Imag ; \
    (c).Imag -= (a).Imag * (b).Real + (a).Real * (b).Imag ; \
}

/* -------------------------------------------------------------------------- */

/* c -= a*conjugate(b), assert because c cannot be the same as a or b */
#define MULT_SUB_CONJ(c,a,b) \
{ \
    ASSERT (&(c) != &(a) && &(c) != &(b)) ; \
    (c).Real -= (a).Real * (b).Real + (a).Imag * (b).Imag ; \
    (c).Imag -= (a).Imag * (b).Real - (a).Real * (b).Imag ; \
}

/* -------------------------------------------------------------------------- */

/* c = a/b, be careful to avoid underflow and overflow */
#ifdef MATHWORKS
#define DIV(c,a,b) \
{ \
    (void) utDivideComplex ((a).Real, (a).Imag, (b).Real, (b).Imag, \
	&((c).Real), &((c).Imag)) ; \
}
#else
/* This uses ACM Algo 116, by R. L. Smith, 1962. */
/* c can be the same variable as a or b. */
/* Ignore NaN case for double relop br>=bi. */
#define DIV(c,a,b) \
{ \
    double r, den, ar, ai, br, bi ; \
    br = (b).Real ; \
    bi = (b).Imag ; \
    ar = (a).Real ; \
    ai = (a).Imag ; \
    if (SCALAR_ABS (br) >= SCALAR_ABS (bi)) \
    { \
	r = bi / br ; \
	den = br + r * bi ; \
	(c).Real = (ar + ai * r) / den ; \
	(c).Imag = (ai - ar * r) / den ; \
    } \
    else \
    { \
	r = br / bi ; \
	den = r * br + bi ; \
	(c).Real = (ar * r + ai) / den ; \
	(c).Imag = (ai * r - ar) / den ; \
    } \
}
#endif

/* -------------------------------------------------------------------------- */

/* c = 1/c, be careful to avoid underflow and overflow */
/* Not used if MATHWORKS is defined. */
/* This uses ACM Algo 116, by R. L. Smith, 1962. */
/* Ignore NaN case for double relop cr>=ci. */
#define RECIPROCAL(c) \
{ \
    double r, den, cr, ci ; \
    cr = (c).Real ; \
    ci = (c).Imag ; \
    if (SCALAR_ABS (cr) >= SCALAR_ABS (ci)) \
    { \
	r = ci / cr ; \
	den = cr + r * ci ; \
	(c).Real = 1.0 / den ; \
	(c).Imag = - r / den ; \
    } \
    else \
    { \
	r = cr / ci ; \
	den = r * cr + ci ; \
	(c).Real = r / den ; \
	(c).Imag = - 1.0 / den ; \
    } \
}


/* -------------------------------------------------------------------------- */

/* c = a/conjugate(b), be careful to avoid underflow and overflow */
#ifdef MATHWORKS
#define DIV_CONJ(c,a,b) \
{ \
    (void) utDivideComplex ((a).Real, (a).Imag, (b).Real, (-(b).Imag), \
	&((c).Real), &((c).Imag)) ; \
}
#else
/* This uses ACM Algo 116, by R. L. Smith, 1962. */
/* c can be the same variable as a or b. */
/* Ignore NaN case for double relop br>=bi. */
#define DIV_CONJ(c,a,b) \
{ \
    double r, den, ar, ai, br, bi ; \
    br = (b).Real ; \
    bi = (b).Imag ; \
    ar = (a).Real ; \
    ai = (a).Imag ; \
    if (SCALAR_ABS (br) >= SCALAR_ABS (bi)) \
    { \
	r = (-bi) / br ; \
	den = br - r * bi ; \
	(c).Real = (ar + ai * r) / den ; \
	(c).Imag = (ai - ar * r) / den ; \
    } \
    else \
    { \
	r = br / (-bi) ; \
	den =  r * br - bi; \
	(c).Real = (ar * r + ai) / den ; \
	(c).Imag = (ai * r - ar) / den ; \
    } \
}
#endif

/* -------------------------------------------------------------------------- */

/* approximate absolute value, s = |r|+|i| */
#define APPROX_ABS(s,a) \
{ \
    (s) = SCALAR_ABS ((a).Real) + SCALAR_ABS ((a).Imag) ; \
}

/* -------------------------------------------------------------------------- */

/* exact absolute value, s = sqrt (a.real^2 + amag^2) */
#ifdef MATHWORKS
#define ABS(s,a) \
{ \
    (s) = utFdlibm_hypot ((a).Real, (a).Imag) ; \
}
#else
/* Ignore NaN case for the double relops ar>=ai and ar+ai==ar. */
#define ABS(s,a) \
{ \
    double r, ar, ai ; \
    ar = SCALAR_ABS ((a).Real) ; \
    ai = SCALAR_ABS ((a).Imag) ; \
    if (ar >= ai) \
    { \
	if (ar + ai == ar) \
	{ \
	    (s) = ar ; \
	} \
	else \
	{ \
	    r = ai / ar ; \
	    (s) = ar * sqrt (1.0 + r*r) ; \
	} \
    } \
    else \
    { \
	if (ai + ar == ai) \
	{ \
	    (s) = ai ; \
	} \
	else \
	{ \
	    r = ar / ai ; \
	    (s) = ai * sqrt (1.0 + r*r) ; \
	} \
    } \
}
#endif

/* -------------------------------------------------------------------------- */

/* print an entry (avoid printing "-0" for negative zero).  */
#define PRINT_ENTRY(a) \
{ \
    if (SCALAR_IS_NONZERO ((a).Real)) \
    { \
	PRINTF ((" (%g", (a).Real)) ; \
    } \
    else \
    { \
	PRINTF ((" (0")) ; \
    } \
    if (SCALAR_IS_LTZERO ((a).Imag)) \
    { \
	PRINTF ((" - %gi)", -(a).Imag)) ; \
    } \
    else if (SCALAR_IS_ZERO ((a).Imag)) \
    { \
	PRINTF ((" + 0i)")) ; \
    } \
    else \
    { \
	PRINTF ((" + %gi)", (a).Imag)) ; \
    } \
}

/* -------------------------------------------------------------------------- */

#endif	/* #ifndef COMPLEX */

/* -------------------------------------------------------------------------- */
/* Double precision, with int's as integers */
/* -------------------------------------------------------------------------- */

#ifdef DINT

#define UMF_analyze		 umf_i_analyze
#define UMF_apply_order		 umf_i_apply_order
#define UMF_assemble		 umfdi_assemble
#define UMF_assemble_fixq	 umfdi_assemble_fixq
#define UMF_blas3_update	 umfdi_blas3_update
#define UMF_build_tuples	 umfdi_build_tuples
#define UMF_build_tuples_usage	 umfdi_build_tuples_usage
#define UMF_colamd		 umf_i_colamd
#define UMF_colamd_set_defaults	 umf_i_colamd_set_defaults
#define UMF_create_element	 umfdi_create_element
#define UMF_extend_front	 umfdi_extend_front
#define UMF_free		 umf_i_free
#define UMF_fsize		 umf_i_fsize
#define UMF_garbage_collection	 umfdi_garbage_collection
#define UMF_get_memory		 umfdi_get_memory
#define UMF_grow_front		 umfdi_grow_front
#define UMF_init_front		 umfdi_init_front
#define UMF_is_permutation	 umf_i_is_permutation
#define UMF_kernel		 umfdi_kernel
#define UMF_kernel_init		 umfdi_kernel_init
#define UMF_kernel_init_usage	 umfdi_kernel_init_usage
#define UMF_kernel_wrapup	 umfdi_kernel_wrapup
#define UMF_local_search	 umfdi_local_search
#define UMF_lsolve		 umfdi_lsolve
#define UMF_ltsolve		 umfdi_ltsolve
#define UMF_lhsolve		 umfdi_lhsolve
#define UMF_malloc		 umf_i_malloc
#define UMF_mem_alloc_element	 umfdi_mem_alloc_element
#define UMF_mem_alloc_head_block umfdi_mem_alloc_head_block
#define UMF_mem_alloc_tail_block umfdi_mem_alloc_tail_block
#define UMF_mem_free_tail_block	 umfdi_mem_free_tail_block
#define UMF_mem_init_memoryspace umfdi_mem_init_memoryspace
#define UMF_realloc		 umf_i_realloc
#define UMF_report_perm		 umf_i_report_perm
#define UMF_report_vector	 umfdi_report_vector
#define UMF_row_search		 umfdi_row_search
#define UMF_scale		 umfdi_scale
#define UMF_scale_column	 umfdi_scale_column
#define UMF_set_stats		 umf_i_set_stats
#define UMF_singletons		 umf_i_singletons
#define UMF_solve		 umfdi_solve
#define UMF_start_front		 umfdi_start_front
#define UMF_store_lu		 umfdi_store_lu
#define UMF_symbolic_usage	 umfdi_symbolic_usage
#define UMF_transpose		 umfdi_transpose
#define UMF_tuple_lengths	 umfdi_tuple_lengths
#define UMF_usolve		 umfdi_usolve
#define UMF_utsolve		 umfdi_utsolve
#define UMF_uhsolve		 umfdi_uhsolve
#define UMF_valid_numeric	 umfdi_valid_numeric
#define UMF_valid_symbolic	 umfdi_valid_symbolic
#define UMF_triplet_map_x	 umfdi_triplet_map_x
#define UMF_triplet_map_nox	 umfdi_triplet_map_nox
#define UMF_triplet_nomap_x	 umfdi_triplet_nomap_x
#define UMF_triplet_nomap_nox	 umfdi_triplet_nomap_nox
#define UMF_2by2		 umfdi_2by2

#define UMFPACK_col_to_triplet	 umfpack_di_col_to_triplet
#define UMFPACK_defaults	 umfpack_di_defaults
#define UMFPACK_free_numeric	 umfpack_di_free_numeric
#define UMFPACK_free_symbolic	 umfpack_di_free_symbolic
#define UMFPACK_get_lunz	 umfpack_di_get_lunz
#define UMFPACK_get_numeric	 umfpack_di_get_numeric
#define UMFPACK_get_symbolic	 umfpack_di_get_symbolic
#define UMFPACK_numeric		 umfpack_di_numeric
#define UMFPACK_qsymbolic	 umfpack_di_qsymbolic
#define UMFPACK_report_control	 umfpack_di_report_control
#define UMFPACK_report_info	 umfpack_di_report_info
#define UMFPACK_report_matrix	 umfpack_di_report_matrix
#define UMFPACK_report_numeric	 umfpack_di_report_numeric
#define UMFPACK_report_perm	 umfpack_di_report_perm
#define UMFPACK_report_status	 umfpack_di_report_status
#define UMFPACK_report_symbolic	 umfpack_di_report_symbolic
#define UMFPACK_report_triplet	 umfpack_di_report_triplet
#define UMFPACK_report_vector	 umfpack_di_report_vector
#define UMFPACK_save_numeric	 umfpack_di_save_numeric
#define UMFPACK_save_symbolic	 umfpack_di_save_symbolic
#define UMFPACK_load_numeric	 umfpack_di_load_numeric
#define UMFPACK_load_symbolic	 umfpack_di_load_symbolic
#define UMFPACK_scale		 umfpack_di_scale
#define UMFPACK_solve		 umfpack_di_solve
#define UMFPACK_symbolic	 umfpack_di_symbolic
#define UMFPACK_transpose	 umfpack_di_transpose
#define UMFPACK_triplet_to_col	 umfpack_di_triplet_to_col
#define UMFPACK_wsolve		 umfpack_di_wsolve

/* for debugging only: */
#define UMF_malloc_count	 umf_i_malloc_count
#define UMF_debug		 umfdi_debug
#define UMF_allocfail		 umfdi_allocfail
#define UMF_gprob		 umfdi_gprob
#define UMF_dump_dense		 umfdi_dump_dense
#define UMF_dump_element	 umfdi_dump_element
#define UMF_dump_rowcol		 umfdi_dump_rowcol
#define UMF_dump_matrix		 umfdi_dump_matrix
#define UMF_dump_current_front	 umfdi_dump_current_front
#define UMF_dump_lu		 umfdi_dump_lu
#define UMF_dump_memory		 umfdi_dump_memory
#define UMF_dump_packed_memory	 umfdi_dump_packed_memory
#define UMF_dump_col_matrix	 umfdi_dump_col_matrix
#define UMF_dump_chain		 umfdi_dump_chain
#define UMF_dump_start		 umfdi_dump_start
#define UMF_dump_rowmerge	 umfdi_dump_rowmerge
#define UMF_dump_diagonal_map	 umfdi_dump_diagonal_map

#endif

/* -------------------------------------------------------------------------- */
/* Double precision, with long's as integers */
/* -------------------------------------------------------------------------- */

#ifdef DLONG

#define UMF_analyze		 umf_l_analyze
#define UMF_apply_order		 umf_l_apply_order
#define UMF_assemble		 umfdl_assemble
#define UMF_assemble_fixq	 umfdl_assemble_fixq
#define UMF_blas3_update	 umfdl_blas3_update
#define UMF_build_tuples	 umfdl_build_tuples
#define UMF_build_tuples_usage	 umfdl_build_tuples_usage
#define UMF_colamd		 umf_l_colamd
#define UMF_colamd_set_defaults	 umf_l_colamd_set_defaults
#define UMF_create_element	 umfdl_create_element
#define UMF_extend_front	 umfdl_extend_front
#define UMF_free		 umf_l_free
#define UMF_fsize		 umf_l_fsize
#define UMF_garbage_collection	 umfdl_garbage_collection
#define UMF_get_memory		 umfdl_get_memory
#define UMF_grow_front		 umfdl_grow_front
#define UMF_init_front		 umfdl_init_front
#define UMF_is_permutation	 umf_l_is_permutation
#define UMF_kernel		 umfdl_kernel
#define UMF_kernel_init		 umfdl_kernel_init
#define UMF_kernel_init_usage	 umfdl_kernel_init_usage
#define UMF_kernel_wrapup	 umfdl_kernel_wrapup
#define UMF_local_search	 umfdl_local_search
#define UMF_lsolve		 umfdl_lsolve
#define UMF_ltsolve		 umfdl_ltsolve
#define UMF_lhsolve		 umfdl_lhsolve
#define UMF_malloc		 umf_l_malloc
#define UMF_mem_alloc_element	 umfdl_mem_alloc_element
#define UMF_mem_alloc_head_block umfdl_mem_alloc_head_block
#define UMF_mem_alloc_tail_block umfdl_mem_alloc_tail_block
#define UMF_mem_free_tail_block	 umfdl_mem_free_tail_block
#define UMF_mem_init_memoryspace umfdl_mem_init_memoryspace
#define UMF_realloc		 umf_l_realloc
#define UMF_report_perm		 umf_l_report_perm
#define UMF_report_vector	 umfdl_report_vector
#define UMF_row_search		 umfdl_row_search
#define UMF_scale		 umfdl_scale
#define UMF_scale_column	 umfdl_scale_column
#define UMF_set_stats		 umf_l_set_stats
#define UMF_singletons		 umf_l_singletons
#define UMF_solve		 umfdl_solve
#define UMF_start_front		 umfdl_start_front
#define UMF_store_lu		 umfdl_store_lu
#define UMF_symbolic_usage	 umfdl_symbolic_usage
#define UMF_transpose		 umfdl_transpose
#define UMF_tuple_lengths	 umfdl_tuple_lengths
#define UMF_usolve		 umfdl_usolve
#define UMF_utsolve		 umfdl_utsolve
#define UMF_uhsolve		 umfdl_uhsolve
#define UMF_valid_numeric	 umfdl_valid_numeric
#define UMF_valid_symbolic	 umfdl_valid_symbolic
#define UMF_triplet_map_x	 umfdl_triplet_map_x
#define UMF_triplet_map_nox	 umfdl_triplet_map_nox
#define UMF_triplet_nomap_x	 umfdl_triplet_nomap_x
#define UMF_triplet_nomap_nox	 umfdl_triplet_nomap_nox
#define UMF_2by2		 umfdl_2by2

#define UMFPACK_col_to_triplet	 umfpack_dl_col_to_triplet
#define UMFPACK_defaults	 umfpack_dl_defaults
#define UMFPACK_free_numeric	 umfpack_dl_free_numeric
#define UMFPACK_free_symbolic	 umfpack_dl_free_symbolic
#define UMFPACK_get_lunz	 umfpack_dl_get_lunz
#define UMFPACK_get_numeric	 umfpack_dl_get_numeric
#define UMFPACK_get_symbolic	 umfpack_dl_get_symbolic
#define UMFPACK_numeric		 umfpack_dl_numeric
#define UMFPACK_qsymbolic	 umfpack_dl_qsymbolic
#define UMFPACK_report_control	 umfpack_dl_report_control
#define UMFPACK_report_info	 umfpack_dl_report_info
#define UMFPACK_report_matrix	 umfpack_dl_report_matrix
#define UMFPACK_report_numeric	 umfpack_dl_report_numeric
#define UMFPACK_report_perm	 umfpack_dl_report_perm
#define UMFPACK_report_status	 umfpack_dl_report_status
#define UMFPACK_report_symbolic	 umfpack_dl_report_symbolic
#define UMFPACK_report_triplet	 umfpack_dl_report_triplet
#define UMFPACK_report_vector	 umfpack_dl_report_vector
#define UMFPACK_save_numeric	 umfpack_dl_save_numeric
#define UMFPACK_save_symbolic	 umfpack_dl_save_symbolic
#define UMFPACK_load_numeric	 umfpack_dl_load_numeric
#define UMFPACK_load_symbolic	 umfpack_dl_load_symbolic
#define UMFPACK_scale		 umfpack_dl_scale
#define UMFPACK_solve		 umfpack_dl_solve
#define UMFPACK_symbolic	 umfpack_dl_symbolic
#define UMFPACK_transpose	 umfpack_dl_transpose
#define UMFPACK_triplet_to_col	 umfpack_dl_triplet_to_col
#define UMFPACK_wsolve		 umfpack_dl_wsolve

/* for debugging only: */
#define UMF_malloc_count	 umf_l_malloc_count
#define UMF_debug		 umfdl_debug
#define UMF_allocfail		 umfdl_allocfail
#define UMF_gprob		 umfdl_gprob
#define UMF_dump_dense		 umfdl_dump_dense
#define UMF_dump_element	 umfdl_dump_element
#define UMF_dump_rowcol		 umfdl_dump_rowcol
#define UMF_dump_matrix		 umfdl_dump_matrix
#define UMF_dump_current_front	 umfdl_dump_current_front
#define UMF_dump_lu		 umfdl_dump_lu
#define UMF_dump_memory		 umfdl_dump_memory
#define UMF_dump_packed_memory	 umfdl_dump_packed_memory
#define UMF_dump_col_matrix	 umfdl_dump_col_matrix
#define UMF_dump_chain		 umfdl_dump_chain
#define UMF_dump_start		 umfdl_dump_start
#define UMF_dump_rowmerge	 umfdl_dump_rowmerge
#define UMF_dump_diagonal_map	 umfdl_dump_diagonal_map

#endif

/* -------------------------------------------------------------------------- */
/* Complex double precision, with int's as integers */
/* -------------------------------------------------------------------------- */

#ifdef ZINT

#define UMF_analyze		 umf_i_analyze
#define UMF_apply_order		 umf_i_apply_order
#define UMF_assemble		 umfzi_assemble
#define UMF_assemble_fixq	 umfzi_assemble_fixq
#define UMF_blas3_update	 umfzi_blas3_update
#define UMF_build_tuples	 umfzi_build_tuples
#define UMF_build_tuples_usage	 umfzi_build_tuples_usage
#define UMF_colamd		 umf_i_colamd
#define UMF_colamd_set_defaults	 umf_i_colamd_set_defaults
#define UMF_create_element	 umfzi_create_element
#define UMF_extend_front	 umfzi_extend_front
#define UMF_free		 umf_i_free
#define UMF_fsize		 umf_i_fsize
#define UMF_garbage_collection	 umfzi_garbage_collection
#define UMF_get_memory		 umfzi_get_memory
#define UMF_grow_front		 umfzi_grow_front
#define UMF_init_front		 umfzi_init_front
#define UMF_is_permutation	 umf_i_is_permutation
#define UMF_kernel		 umfzi_kernel
#define UMF_kernel_init		 umfzi_kernel_init
#define UMF_kernel_init_usage	 umfzi_kernel_init_usage
#define UMF_kernel_wrapup	 umfzi_kernel_wrapup
#define UMF_local_search	 umfzi_local_search
#define UMF_lsolve		 umfzi_lsolve
#define UMF_ltsolve		 umfzi_ltsolve
#define UMF_lhsolve		 umfzi_lhsolve
#define UMF_malloc		 umf_i_malloc
#define UMF_mem_alloc_element	 umfzi_mem_alloc_element
#define UMF_mem_alloc_head_block umfzi_mem_alloc_head_block
#define UMF_mem_alloc_tail_block umfzi_mem_alloc_tail_block
#define UMF_mem_free_tail_block	 umfzi_mem_free_tail_block
#define UMF_mem_init_memoryspace umfzi_mem_init_memoryspace
#define UMF_realloc		 umf_i_realloc
#define UMF_report_perm		 umf_i_report_perm
#define UMF_report_vector	 umfzi_report_vector
#define UMF_row_search		 umfzi_row_search
#define UMF_scale		 umfzi_scale
#define UMF_scale_column	 umfzi_scale_column
#define UMF_set_stats		 umfzi_set_stats
#define UMF_singletons		 umf_i_singletons
#define UMF_solve		 umfzi_solve
#define UMF_start_front		 umfzi_start_front
#define UMF_store_lu		 umfzi_store_lu
#define UMF_symbolic_usage	 umfzi_symbolic_usage
#define UMF_transpose		 umfzi_transpose
#define UMF_tuple_lengths	 umfzi_tuple_lengths
#define UMF_usolve		 umfzi_usolve
#define UMF_utsolve		 umfzi_utsolve
#define UMF_uhsolve		 umfzi_uhsolve
#define UMF_valid_numeric	 umfzi_valid_numeric
#define UMF_valid_symbolic	 umfzi_valid_symbolic
#define UMF_triplet_map_x	 umfzi_triplet_map_x
#define UMF_triplet_map_nox	 umfzi_triplet_map_nox
#define UMF_triplet_nomap_x	 umfzi_triplet_nomap_x
#define UMF_triplet_nomap_nox	 umfzi_triplet_nomap_nox
#define UMF_2by2		 umfzi_2by2

#define UMFPACK_col_to_triplet	 umfpack_zi_col_to_triplet
#define UMFPACK_defaults	 umfpack_zi_defaults
#define UMFPACK_free_numeric	 umfpack_zi_free_numeric
#define UMFPACK_free_symbolic	 umfpack_zi_free_symbolic
#define UMFPACK_get_lunz	 umfpack_zi_get_lunz
#define UMFPACK_get_numeric	 umfpack_zi_get_numeric
#define UMFPACK_get_symbolic	 umfpack_zi_get_symbolic
#define UMFPACK_numeric		 umfpack_zi_numeric
#define UMFPACK_qsymbolic	 umfpack_zi_qsymbolic
#define UMFPACK_report_control	 umfpack_zi_report_control
#define UMFPACK_report_info	 umfpack_zi_report_info
#define UMFPACK_report_matrix	 umfpack_zi_report_matrix
#define UMFPACK_report_numeric	 umfpack_zi_report_numeric
#define UMFPACK_report_perm	 umfpack_zi_report_perm
#define UMFPACK_report_status	 umfpack_zi_report_status
#define UMFPACK_report_symbolic	 umfpack_zi_report_symbolic
#define UMFPACK_report_triplet	 umfpack_zi_report_triplet
#define UMFPACK_report_vector	 umfpack_zi_report_vector
#define UMFPACK_save_numeric	 umfpack_zi_save_numeric
#define UMFPACK_save_symbolic	 umfpack_zi_save_symbolic
#define UMFPACK_load_numeric	 umfpack_zi_load_numeric
#define UMFPACK_load_symbolic	 umfpack_zi_load_symbolic
#define UMFPACK_scale		 umfpack_zi_scale
#define UMFPACK_solve		 umfpack_zi_solve
#define UMFPACK_symbolic	 umfpack_zi_symbolic
#define UMFPACK_transpose	 umfpack_zi_transpose
#define UMFPACK_triplet_to_col	 umfpack_zi_triplet_to_col
#define UMFPACK_wsolve		 umfpack_zi_wsolve

/* for debugging only: */
#define UMF_malloc_count	 umf_i_malloc_count
#define UMF_debug		 umfzi_debug
#define UMF_allocfail		 umfzi_allocfail
#define UMF_gprob		 umfzi_gprob
#define UMF_dump_dense		 umfzi_dump_dense
#define UMF_dump_element	 umfzi_dump_element
#define UMF_dump_rowcol		 umfzi_dump_rowcol
#define UMF_dump_matrix		 umfzi_dump_matrix
#define UMF_dump_current_front	 umfzi_dump_current_front
#define UMF_dump_lu		 umfzi_dump_lu
#define UMF_dump_memory		 umfzi_dump_memory
#define UMF_dump_packed_memory	 umfzi_dump_packed_memory
#define UMF_dump_col_matrix	 umfzi_dump_col_matrix
#define UMF_dump_chain		 umfzi_dump_chain
#define UMF_dump_start		 umfzi_dump_start
#define UMF_dump_rowmerge	 umfzi_dump_rowmerge
#define UMF_dump_diagonal_map	 umfzi_dump_diagonal_map

#endif

/* -------------------------------------------------------------------------- */
/* Complex double precision, with long's as integers */
/* -------------------------------------------------------------------------- */

#ifdef ZLONG

#define UMF_analyze		 umf_l_analyze
#define UMF_apply_order		 umf_l_apply_order
#define UMF_assemble		 umfzl_assemble
#define UMF_assemble_fixq	 umfzl_assemble_fixq
#define UMF_blas3_update	 umfzl_blas3_update
#define UMF_build_tuples	 umfzl_build_tuples
#define UMF_build_tuples_usage	 umfzl_build_tuples_usage
#define UMF_colamd		 umf_l_colamd
#define UMF_colamd_set_defaults	 umf_l_colamd_set_defaults
#define UMF_create_element	 umfzl_create_element
#define UMF_extend_front	 umfzl_extend_front
#define UMF_free		 umf_l_free
#define UMF_fsize		 umf_l_fsize
#define UMF_garbage_collection	 umfzl_garbage_collection
#define UMF_get_memory		 umfzl_get_memory
#define UMF_grow_front		 umfzl_grow_front
#define UMF_init_front		 umfzl_init_front
#define UMF_is_permutation	 umf_l_is_permutation
#define UMF_kernel		 umfzl_kernel
#define UMF_kernel_init		 umfzl_kernel_init
#define UMF_kernel_init_usage	 umfzl_kernel_init_usage
#define UMF_kernel_wrapup	 umfzl_kernel_wrapup
#define UMF_local_search	 umfzl_local_search
#define UMF_lsolve		 umfzl_lsolve
#define UMF_ltsolve		 umfzl_ltsolve
#define UMF_lhsolve		 umfzl_lhsolve
#define UMF_malloc		 umf_l_malloc
#define UMF_mem_alloc_element	 umfzl_mem_alloc_element
#define UMF_mem_alloc_head_block umfzl_mem_alloc_head_block
#define UMF_mem_alloc_tail_block umfzl_mem_alloc_tail_block
#define UMF_mem_free_tail_block	 umfzl_mem_free_tail_block
#define UMF_mem_init_memoryspace umfzl_mem_init_memoryspace
#define UMF_realloc		 umf_l_realloc
#define UMF_report_perm		 umf_l_report_perm
#define UMF_report_vector	 umfzl_report_vector
#define UMF_row_search		 umfzl_row_search
#define UMF_scale		 umfzl_scale
#define UMF_scale_column	 umfzl_scale_column
#define UMF_set_stats		 umfzl_set_stats
#define UMF_singletons		 umf_l_singletons
#define UMF_solve		 umfzl_solve
#define UMF_start_front		 umfzl_start_front
#define UMF_store_lu		 umfzl_store_lu
#define UMF_symbolic_usage	 umfzl_symbolic_usage
#define UMF_transpose		 umfzl_transpose
#define UMF_tuple_lengths	 umfzl_tuple_lengths
#define UMF_usolve		 umfzl_usolve
#define UMF_utsolve		 umfzl_utsolve
#define UMF_uhsolve		 umfzl_uhsolve
#define UMF_valid_numeric	 umfzl_valid_numeric
#define UMF_valid_symbolic	 umfzl_valid_symbolic
#define UMF_triplet_map_x	 umfzl_triplet_map_x
#define UMF_triplet_map_nox	 umfzl_triplet_map_nox
#define UMF_triplet_nomap_x	 umfzl_triplet_nomap_x
#define UMF_triplet_nomap_nox	 umfzl_triplet_nomap_nox
#define UMF_2by2		 umfzl_2by2

#define UMFPACK_col_to_triplet	 umfpack_zl_col_to_triplet
#define UMFPACK_defaults	 umfpack_zl_defaults
#define UMFPACK_free_numeric	 umfpack_zl_free_numeric
#define UMFPACK_free_symbolic	 umfpack_zl_free_symbolic
#define UMFPACK_get_lunz	 umfpack_zl_get_lunz
#define UMFPACK_get_numeric	 umfpack_zl_get_numeric
#define UMFPACK_get_symbolic	 umfpack_zl_get_symbolic
#define UMFPACK_numeric		 umfpack_zl_numeric
#define UMFPACK_qsymbolic	 umfpack_zl_qsymbolic
#define UMFPACK_report_control	 umfpack_zl_report_control
#define UMFPACK_report_info	 umfpack_zl_report_info
#define UMFPACK_report_matrix	 umfpack_zl_report_matrix
#define UMFPACK_report_numeric	 umfpack_zl_report_numeric
#define UMFPACK_report_perm	 umfpack_zl_report_perm
#define UMFPACK_report_status	 umfpack_zl_report_status
#define UMFPACK_report_symbolic	 umfpack_zl_report_symbolic
#define UMFPACK_report_triplet	 umfpack_zl_report_triplet
#define UMFPACK_report_vector	 umfpack_zl_report_vector
#define UMFPACK_save_numeric	 umfpack_zl_save_numeric
#define UMFPACK_save_symbolic	 umfpack_zl_save_symbolic
#define UMFPACK_load_numeric	 umfpack_zl_load_numeric
#define UMFPACK_load_symbolic	 umfpack_zl_load_symbolic
#define UMFPACK_scale		 umfpack_zl_scale
#define UMFPACK_solve		 umfpack_zl_solve
#define UMFPACK_symbolic	 umfpack_zl_symbolic
#define UMFPACK_transpose	 umfpack_zl_transpose
#define UMFPACK_triplet_to_col	 umfpack_zl_triplet_to_col
#define UMFPACK_wsolve		 umfpack_zl_wsolve

/* for debugging only: */
#define UMF_malloc_count	 umf_l_malloc_count
#define UMF_debug		 umfzl_debug
#define UMF_allocfail		 umfzl_allocfail
#define UMF_gprob		 umfzl_gprob
#define UMF_dump_dense		 umfzl_dump_dense
#define UMF_dump_element	 umfzl_dump_element
#define UMF_dump_rowcol		 umfzl_dump_rowcol
#define UMF_dump_matrix		 umfzl_dump_matrix
#define UMF_dump_current_front	 umfzl_dump_current_front
#define UMF_dump_lu		 umfzl_dump_lu
#define UMF_dump_memory		 umfzl_dump_memory
#define UMF_dump_packed_memory	 umfzl_dump_packed_memory
#define UMF_dump_col_matrix	 umfzl_dump_col_matrix
#define UMF_dump_chain		 umfzl_dump_chain
#define UMF_dump_start		 umfzl_dump_start
#define UMF_dump_rowmerge	 umfzl_dump_rowmerge
#define UMF_dump_diagonal_map	 umfzl_dump_diagonal_map

#endif
