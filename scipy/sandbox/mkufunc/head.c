
#include <math.h>

/* ================================================== g_prerequisite.h === */

typedef unsigned char bool_t;

/* ================================================== exception.h ======== */

#define RPY_DEBUG_RETURN()        /* nothing */


/* ================================================== int.h ============== */

/*** unary operations ***/

#define OP_INT_IS_TRUE(x,r)   OP_INT_NE(x,0,r)

#define OP_INT_INVERT(x,r)    r = ~((x))

#define OP_INT_NEG(x,r)    r = -(x)

#define OP_INT_NEG_OVF(x,r) \
    if ((x) == LONG_MIN) FAIL_OVF("integer negate"); \
	OP_INT_NEG(x,r)
#define OP_LLONG_NEG_OVF(x,r) \
    if ((x) == LLONG_MIN) FAIL_OVF("integer negate"); \
	OP_LLONG_NEG(x,r)

#define OP_INT_ABS(x,r)    r = (x) >= 0 ? x : -(x)

#define OP_INT_ABS_OVF(x,r) \
    if ((x) == LONG_MIN) FAIL_OVF("integer absolute"); \
	OP_INT_ABS(x,r)
#define OP_LLONG_ABS_OVF(x,r) \
    if ((x) == LLONG_MIN) FAIL_OVF("integer absolute"); \
	OP_LLONG_ABS(x,r)

/***  binary operations ***/

#define OP_INT_EQ(x,y,r)	  r = ((x) == (y))
#define OP_INT_NE(x,y,r)	  r = ((x) != (y))
#define OP_INT_LE(x,y,r)	  r = ((x) <= (y))
#define OP_INT_GT(x,y,r)	  r = ((x) >  (y))
#define OP_INT_LT(x,y,r)	  r = ((x) <  (y))
#define OP_INT_GE(x,y,r)	  r = ((x) >= (y))

/* addition, subtraction */

#define OP_INT_ADD(x,y,r)     r = (x) + (y)

#define OP_INT_ADD_OVF(x,y,r) \
	OP_INT_ADD(x,y,r); \
	if ((r^(x)) >= 0 || (r^(y)) >= 0); \
	else FAIL_OVF("integer addition")

#define OP_INT_ADD_NONNEG_OVF(x,y,r)  /* y can be assumed >= 0 */ \
    OP_INT_ADD(x,y,r); \
    if (r >= (x)); \
    else FAIL_OVF("integer addition")
/* XXX can a C compiler be too clever and think it can "prove" that
 * r >= x always hold above? */

#define OP_INT_SUB(x,y,r)     r = (x) - (y)

#define OP_INT_SUB_OVF(x,y,r) \
	OP_INT_SUB(x,y,r); \
	if ((r^(x)) >= 0 || (r^~(y)) >= 0); \
	else FAIL_OVF("integer subtraction")

#define OP_INT_MUL(x,y,r)     r = (x) * (y)

#if defined(HAVE_LONG_LONG) && SIZE_OF_LONG_LONG < SIZE_OF_LONG
#  define OP_INT_MUL_OVF_LL      1
#lse
#  define OP_INT_MUL_OVF_LL      0
#endif

#if !OP_INT_MUL_OVF_LL

#define OP_INT_MUL_OVF(x,y,r) \
	if (op_int_mul_ovf(x,y,&r)); \
	else FAIL_OVF("integer multiplication")

#else

#define OP_INT_MUL_OVF(x,y,r) \
	{ \
		PY_LONG_LONG lr = (PY_LONG_LONG)(x) * (PY_LONG_LONG)(y); \
		r = (long)lr; \
		if ((PY_LONG_LONG)r == lr); \
		else FAIL_OVF("integer multiplication"); \
	}
#endif

/* shifting */

/* NB. shifting has same limitations as C: the shift count must be
       >= 0 and < LONG_BITS. */
#define OP_INT_RSHIFT(x,y,r)    r = Py_ARITHMETIC_RIGHT_SHIFT(long, x, y)
#define OP_UINT_RSHIFT(x,y,r)   r = (x) >> (y)
#define OP_LLONG_RSHIFT(x,y,r)  r = Py_ARITHMETIC_RIGHT_SHIFT(PY_LONG_LONG,x,y)
#define OP_ULLONG_RSHIFT(x,y,r) r = (x) >> (y)

#define OP_INT_LSHIFT(x,y,r)    r = (x) << (y)
#define OP_UINT_LSHIFT(x,y,r)   r = (x) << (y)
#define OP_LLONG_LSHIFT(x,y,r)  r = (x) << (y)
#define OP_ULLONG_LSHIFT(x,y,r) r = (x) << (y)

#define OP_INT_LSHIFT_OVF(x,y,r) \
	OP_INT_LSHIFT(x,y,r); \
	if ((x) != Py_ARITHMETIC_RIGHT_SHIFT(long, r, (y))) \
		FAIL_OVF("x<<y losing bits or changing sign")

/* the safe value-checking version of the above macros */

#define OP_INT_RSHIFT_VAL(x,y,r) \
	if ((y) >= 0) { OP_INT_RSHIFT(x,y,r); } \
	else FAIL_VAL("negative shift count")
#define OP_LLONG_RSHIFT_VAL(x,y,r) \
	if ((y) >= 0) { OP_LLONG_RSHIFT(x,y,r); } \
	else FAIL_VAL("negative shift count")

#define OP_INT_LSHIFT_VAL(x,y,r) \
	if ((y) >= 0) { OP_INT_LSHIFT(x,y,r); } \
	else FAIL_VAL("negative shift count")
#define OP_LLONG_LSHIFT_VAL(x,y,r) \
	if ((y) >= 0) { OP_LLONG_LSHIFT(x,y,r); } \
	else FAIL_VAL("negative shift count")

#define OP_INT_LSHIFT_OVF_VAL(x,y,r) \
	if ((y) >= 0) { OP_INT_LSHIFT_OVF(x,y,r); } \
	else FAIL_VAL("negative shift count")

/* pff */
#define OP_UINT_LSHIFT_VAL(x,y,r) \
	if ((y) >= 0) { OP_UINT_LSHIFT(x,y,r); } \
	else FAIL_VAL("negative shift count")
#define OP_ULLONG_LSHIFT_VAL(x,y,r) \
	if ((y) >= 0) { OP_ULLONG_LSHIFT(x,y,r); } \
	else FAIL_VAL("negative shift count")

#define OP_UINT_RSHIFT_VAL(x,y,r) \
	if ((y) >= 0) { OP_UINT_RSHIFT(x,y,r); } \
	else FAIL_VAL("negative shift count")
#define OP_ULLONG_RSHIFT_VAL(x,y,r) \
	if ((y) >= 0) { OP_ULLONG_RSHIFT(x,y,r); } \
	else FAIL_VAL("negative shift count")


/* floor division */

#define OP_INT_FLOORDIV(x,y,r)    r = (x) / (y)
#define OP_UINT_FLOORDIV(x,y,r)   r = (x) / (y)
#define OP_LLONG_FLOORDIV(x,y,r)  r = (x) / (y)
#define OP_ULLONG_FLOORDIV(x,y,r) r = (x) / (y)

#define OP_INT_FLOORDIV_OVF(x,y,r) \
	if ((y) == -1 && (x) == LONG_MIN) \
            { FAIL_OVF("integer division"); } \
        else OP_INT_FLOORDIV(x,y,r)

#define OP_INT_FLOORDIV_ZER(x,y,r) \
	if ((y)) { OP_INT_FLOORDIV(x,y,r); } \
	else FAIL_ZER("integer division")
#define OP_UINT_FLOORDIV_ZER(x,y,r) \
	if ((y)) { OP_UINT_FLOORDIV(x,y,r); } \
	else FAIL_ZER("unsigned integer division")
#define OP_LLONG_FLOORDIV_ZER(x,y,r) \
	if ((y)) { OP_LLONG_FLOORDIV(x,y,r); } \
	else FAIL_ZER("integer division")
#define OP_ULLONG_FLOORDIV_ZER(x,y,r) \
	if ((y)) { OP_ULLONG_FLOORDIV(x,y,r); } \
	else FAIL_ZER("unsigned integer division")

#define OP_INT_FLOORDIV_OVF_ZER(x,y,r) \
	if ((y)) { OP_INT_FLOORDIV_OVF(x,y,r); } \
	else FAIL_ZER("integer division")

/* modulus */

#define OP_INT_MOD(x,y,r)     r = (x) % (y)
#define OP_UINT_MOD(x,y,r)    r = (x) % (y)
#define OP_LLONG_MOD(x,y,r)   r = (x) % (y)
#define OP_ULLONG_MOD(x,y,r)  r = (x) % (y)

#define OP_INT_MOD_OVF(x,y,r) \
	if ((y) == -1 && (x) == LONG_MIN) \
            { FAIL_OVF("integer modulo"); }\
        else OP_INT_MOD(x,y,r)

#define OP_INT_MOD_ZER(x,y,r) \
	if ((y)) { OP_INT_MOD(x,y,r); } \
	else FAIL_ZER("integer modulo")
#define OP_UINT_MOD_ZER(x,y,r) \
	if ((y)) { OP_UINT_MOD(x,y,r); } \
	else FAIL_ZER("unsigned integer modulo")
#define OP_LLONG_MOD_ZER(x,y,r) \
	if ((y)) { OP_LLONG_MOD(x,y,r); } \
	else FAIL_ZER("integer modulo")
#define OP_ULLONG_MOD_ZER(x,y,r) \
	if ((y)) { OP_ULLONG_MOD(x,y,r); } \
	else FAIL_ZER("integer modulo")

#define OP_INT_MOD_OVF_ZER(x,y,r) \
	if ((y)) { OP_INT_MOD_OVF(x,y,r); } \
	else FAIL_ZER("integer modulo")

/* bit operations */

#define OP_INT_AND(x,y,r)     r = (x) & (y)
#define OP_INT_OR( x,y,r)     r = (x) | (y)
#define OP_INT_XOR(x,y,r)     r = (x) ^ (y)

/*** conversions ***/

#define OP_CAST_BOOL_TO_INT(x,r)    r = (long)(x)
#define OP_CAST_BOOL_TO_UINT(x,r)   r = (unsigned long)(x)
#define OP_CAST_UINT_TO_INT(x,r)    r = (long)(x)
#define OP_CAST_INT_TO_UINT(x,r)    r = (unsigned long)(x)
#define OP_CAST_INT_TO_LONGLONG(x,r) r = (long long)(x)
#define OP_CAST_CHAR_TO_INT(x,r)    r = (long)((unsigned char)(x))
#define OP_CAST_INT_TO_CHAR(x,r)    r = (char)(x)
#define OP_CAST_PTR_TO_INT(x,r)     r = (long)(x)    /* XXX */

#define OP_TRUNCATE_LONGLONG_TO_INT(x,r) r = (long)(x)

#define OP_CAST_UNICHAR_TO_INT(x,r)    r = (long)((unsigned long)(x)) /*?*/
#define OP_CAST_INT_TO_UNICHAR(x,r)    r = (unsigned int)(x)

/* bool operations */

#define OP_BOOL_NOT(x, r) r = !(x)

/* _________________ certain implementations __________________ */

#if !OP_INT_MUL_OVF_LL
/* adjusted from intobject.c, Python 2.3.3 */

/* prototypes */

int op_int_mul_ovf(long a, long b, long *longprod);

/* implementations */

#ifndef PYPY_NOT_MAIN_FILE

int
op_int_mul_ovf(long a, long b, long *longprod)
{
	double doubled_longprod;	/* (double)longprod */
	double doubleprod;		/* (double)a * (double)b */

	*longprod = a * b;
	doubleprod = (double)a * (double)b;
	doubled_longprod = (double)*longprod;

	/* Fast path for normal case:  small multiplicands, and no info
	   is lost in either method. */
	if (doubled_longprod == doubleprod)
		return 1;

	/* Somebody somewhere lost info.  Close enough, or way off?  Note
	   that a != 0 and b != 0 (else doubled_longprod == doubleprod == 0).
	   The difference either is or isn't significant compared to the
	   true value (of which doubleprod is a good approximation).
	*/
	{
		const double diff = doubled_longprod - doubleprod;
		const double absdiff = diff >= 0.0 ? diff : -diff;
		const double absprod = doubleprod >= 0.0 ? doubleprod :
							  -doubleprod;
		/* absdiff/absprod <= 1/32 iff
		   32 * absdiff <= absprod -- 5 good bits is "close enough" */
		if (32.0 * absdiff <= absprod)
			return 1;
		return 0;
	}
}

#endif /* PYPY_NOT_MAIN_FILE */

#endif /* !OP_INT_MUL_OVF_LL */

/* implementations */

#define OP_UINT_IS_TRUE OP_INT_IS_TRUE
#define OP_UINT_INVERT OP_INT_INVERT
#define OP_UINT_ADD OP_INT_ADD
#define OP_UINT_SUB OP_INT_SUB
#define OP_UINT_MUL OP_INT_MUL
#define OP_UINT_LT OP_INT_LT
#define OP_UINT_LE OP_INT_LE
#define OP_UINT_EQ OP_INT_EQ
#define OP_UINT_NE OP_INT_NE
#define OP_UINT_GT OP_INT_GT
#define OP_UINT_GE OP_INT_GE
#define OP_UINT_AND OP_INT_AND
#define OP_UINT_OR OP_INT_OR
#define OP_UINT_XOR OP_INT_XOR

#define OP_LLONG_IS_TRUE OP_INT_IS_TRUE
#define OP_LLONG_NEG     OP_INT_NEG
#define OP_LLONG_ABS     OP_INT_ABS
#define OP_LLONG_INVERT  OP_INT_INVERT

#define OP_LLONG_ADD OP_INT_ADD
#define OP_LLONG_SUB OP_INT_SUB
#define OP_LLONG_MUL OP_INT_MUL
#define OP_LLONG_LT  OP_INT_LT
#define OP_LLONG_LE  OP_INT_LE
#define OP_LLONG_EQ  OP_INT_EQ
#define OP_LLONG_NE  OP_INT_NE
#define OP_LLONG_GT  OP_INT_GT
#define OP_LLONG_GE  OP_INT_GE
#define OP_LLONG_AND    OP_INT_AND
#define OP_LLONG_OR     OP_INT_OR
#define OP_LLONG_XOR    OP_INT_XOR

#define OP_ULLONG_IS_TRUE OP_LLONG_IS_TRUE
#define OP_ULLONG_INVERT  OP_LLONG_INVERT
#define OP_ULLONG_ADD OP_LLONG_ADD
#define OP_ULLONG_SUB OP_LLONG_SUB
#define OP_ULLONG_MUL OP_LLONG_MUL
#define OP_ULLONG_LT OP_LLONG_LT
#define OP_ULLONG_LE OP_LLONG_LE
#define OP_ULLONG_EQ OP_LLONG_EQ
#define OP_ULLONG_NE OP_LLONG_NE
#define OP_ULLONG_GT OP_LLONG_GT
#define OP_ULLONG_GE OP_LLONG_GE
#define OP_ULLONG_AND OP_LLONG_AND
#define OP_ULLONG_OR OP_LLONG_OR
#define OP_ULLONG_XOR OP_LLONG_XOR

/* ================================================== float.h ============ */

/*** unary operations ***/

#define OP_FLOAT_IS_TRUE(x,r)   OP_FLOAT_NE(x,0.0,r)
#define OP_FLOAT_NEG(x,r)       r = -x
#define OP_FLOAT_ABS(x,r)       r = fabs(x)

/***  binary operations ***/

#define OP_FLOAT_EQ(x,y,r)	  r = (x == y)
#define OP_FLOAT_NE(x,y,r)	  r = (x != y)
#define OP_FLOAT_LE(x,y,r)	  r = (x <= y)
#define OP_FLOAT_GT(x,y,r)	  r = (x >  y)
#define OP_FLOAT_LT(x,y,r)	  r = (x <  y)
#define OP_FLOAT_GE(x,y,r)	  r = (x >= y)

#define OP_FLOAT_CMP(x,y,r) \
	r = ((x > y) - (x < y))

/* addition, subtraction */

#define OP_FLOAT_ADD(x,y,r)     r = x + y
#define OP_FLOAT_SUB(x,y,r)     r = x - y
#define OP_FLOAT_MUL(x,y,r)     r = x * y
#define OP_FLOAT_TRUEDIV(x,y,r) r = x / y
#define OP_FLOAT_POW(x,y,r)     r = pow(x, y) 

/*** conversions ***/

#define OP_CAST_FLOAT_TO_INT(x,r)    r = (long)(x)
#define OP_CAST_FLOAT_TO_UINT(x,r)   r = (unsigned long)(x)
#define OP_CAST_INT_TO_FLOAT(x,r)    r = (double)(x)
#define OP_CAST_UINT_TO_FLOAT(x,r)   r = (double)(x)
#define OP_CAST_LONGLONG_TO_FLOAT(x,r) r = (double)(x)
#define OP_CAST_BOOL_TO_FLOAT(x,r)   r = (double)(x)

#ifdef HAVE_LONG_LONG
#define OP_CAST_FLOAT_TO_LONGLONG(x,r) r = (long long)(x)
#endif




/* ================================================== EOF ================ */
/* ================================================== EOF ================ */
