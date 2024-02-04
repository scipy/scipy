
#ifndef SUPERLU_CONFIG_H
#define SUPERLU_CONFIG_H

/* Enable metis */
/* #undef HAVE_METIS */

/* Enable colamd */
/* #undef HAVE_COLAMD */

/* enable 64bit index mode */
/* #undef XSDK_INDEX_SIZE */

/* Integer type for indexing sparse matrix meta structure */
#if (XSDK_INDEX_SIZE == 64)
#include <stdint.h>
#define _LONGINT 1
typedef int64_t int_t;
#else
typedef int int_t; /* default */
#endif

#endif /* SUPERLU_CONFIG_H */

