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

extern void F_FUNC(dstevr,DSTEVR)(char *jobz, char *range, int *n, double *d, double *e,
				  double *vl, double *vu, int *il, int *iu, double *abstol,
				  int *m, double *w, double *z, int *ldz, int *isuppz,
				  double *work, int *lwork, int *iwork, int *liwork, int *info);

static void c_dstevr(char *jobz, char *range, int *n, double *d, double *e,
                     double *vl, double *vu, int *il, int *iu, double *abstol,
                     int *m, double *w, double *z, int *ldz, int *isuppz,
                     double *work, int *lwork, int *iwork, int *liwork, int *info) {
    F_FUNC(dstevr,DSTEVR)(jobz, range, n, d, e, vl, vu, il, iu, abstol, m,
			  w, z, ldz, isuppz, work, lwork, iwork, liwork, info);
}
