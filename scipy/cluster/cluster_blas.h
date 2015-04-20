/*
 * Handle different Fortran conventions.
 */

#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F
#else
#define F_FUNC(f,F) f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F##_
#else
#define F_FUNC(f,F) f##_
#endif
#endif

#define f_sgemm F_FUNC(sgemm,SGEMM)
#define f_dgemm F_FUNC(dgemm,DGEMM)
