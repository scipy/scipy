/* #define SIZEOF_SHORT 2 */
/* #define SIZEOF_INT 4 */
/* #define SIZEOF_LONG 4 */
/* #define SIZEOF_FLOAT 4 */
/* #define SIZEOF_DOUBLE 8 */
#define SIZEOF_LONG_DOUBLE 12
#define SIZEOF_PY_INTPTR_T 4
/* #define SIZEOF_LONG_LONG 8 */
#define SIZEOF_PY_LONG_LONG 8
/* #define CHAR_BIT 8 */
#define MATHLIB 
#define HAVE_FLOAT_FUNCS
#define HAVE_LOG1P
#define HAVE_EXPM1
#define HAVE_INVERSE_HYPERBOLIC
#define HAVE_INVERSE_HYPERBOLIC_FLOAT
#define HAVE_RINT
#ifdef WITH_THREAD
#define NPY_ALLOW_THREADS 1
#else
#define NPY_ALLOW_THREADS 0
#endif
