/*
 * Handle qh_new_qhull_scipy entry point.
 */
#ifndef QHULL_MISC_H_
#define QHULL_MISC_H_

/* for CBLAS_INT only*/
#include "npy_cblas.h"

#define qhull_misc_lib_check() QHULL_LIB_CHECK

#include <libqhull_r/libqhull_r.h>

int qh_new_qhull_scipy(qhT *qh, int dim, int numpoints, coordT *points, boolT ismalloc,
                       char *qhull_cmd, FILE *outfile, FILE *errfile, coordT* feaspoint);

#endif /* QHULL_MISC_H_ */
