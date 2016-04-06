/*
 * Handle different Fortran conventions.
 */

#include "qhull_misc_config.h"

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

#define qh_dgetrf F_FUNC(dgetrf,DGETRF)
#define qh_dgecon F_FUNC(dgecon,DGECON)
#define qh_dgetrs F_FUNC(dgetrs,DGETRS)

#define qhull_misc_lib_check() QHULL_LIB_CHECK

#if HAVE_OPEN_MEMSTREAM
FILE *qhull_open_memstream(char **ptr, size_t *sizeloc)
{
    return open_memstream(ptr, sizeloc);
}
#else
FILE *qhull_open_memstream(char **ptr, size_t *sizeloc)
{
    return NULL;
}
#endif
