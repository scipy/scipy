
/*! @file slu_dcomplex.h
 * \brief Header file for complex operations
 * <pre> 
 *  -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 * Contains definitions for various complex operations.
 * This header file is to be included in source files z*.c
 * </pre>
 */
#ifndef __SUPERLU_DCOMPLEX /* allow multiple inclusions */
#define __SUPERLU_DCOMPLEX


#ifndef DCOMPLEX_INCLUDE
#define DCOMPLEX_INCLUDE

typedef struct { double r, i; } doublecomplex;


/* Macro definitions */

/*! \brief Complex Addition c = a + b */
#define z_add(c, a, b) { (c)->r = (a)->r + (b)->r; \
			 (c)->i = (a)->i + (b)->i; }

/*! \brief Complex Subtraction c = a - b */
#define z_sub(c, a, b) { (c)->r = (a)->r - (b)->r; \
			 (c)->i = (a)->i - (b)->i; }

/*! \brief Complex-Double Multiplication */
#define zd_mult(c, a, b) { (c)->r = (a)->r * (b); \
                           (c)->i = (a)->i * (b); }

/*! \brief Complex-Complex Multiplication */
#define zz_mult(c, a, b) { \
	double cr, ci; \
    	cr = (a)->r * (b)->r - (a)->i * (b)->i; \
    	ci = (a)->i * (b)->r + (a)->r * (b)->i; \
    	(c)->r = cr; \
    	(c)->i = ci; \
    }

#define zz_conj(a, b) { \
        (a)->r = (b)->r; \
        (a)->i = -((b)->i); \
    }

/*! \brief Complex equality testing */
#define z_eq(a, b)  ( (a)->r == (b)->r && (a)->i == (b)->i )


#ifdef __cplusplus
extern "C" {
#endif

/* Prototypes for functions in dcomplex.c */
void z_div(doublecomplex *, doublecomplex *, doublecomplex *);
double z_abs(doublecomplex *);     /* exact */
double z_abs1(doublecomplex *);    /* approximate */
void z_exp(doublecomplex *, doublecomplex *);
void d_cnjg(doublecomplex *r, doublecomplex *z);
double d_imag(doublecomplex *);
doublecomplex z_sgn(doublecomplex *);
doublecomplex z_sqrt(doublecomplex *);



#ifdef __cplusplus
  }
#endif

#endif

#endif  /* __SUPERLU_DCOMPLEX */
