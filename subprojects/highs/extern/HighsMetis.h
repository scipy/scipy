/*
 * Interface to the version of METIS bundled with HiGHS, which is compiled
 * only when HIPO is enabled.
 *
 * This header lives outside extern/metis so that it is available even when
 * extern/metis is absent (e.g. SciPy trims it from its sdist), which keeps
 * HighsExtrasApi layout-compatible with an extras library built with HIPO.
 * It is also included by extern/metis/metis.h, so is the only definition of
 * the METIS interface used by HiGHS.
 */

#ifndef HIGHS_METIS_H
#define HIGHS_METIS_H

#include <stdint.h>

#include "HConfig.h"

/* Width of the integer type holding vertices and adjacency lists */
#ifdef HIGHSINT64
#define IDXTYPEWIDTH 64
typedef int64_t idx_t;
#else
#define IDXTYPEWIDTH 32
typedef int32_t idx_t;
#endif

/* The maximum length of the options[] array */
#define METIS_NOPTIONS 13

#ifdef _WINDLL
#define METIS_API(type) __declspec(dllexport) type __cdecl
#elif defined(__cdecl)
#define METIS_API(type) type __cdecl
#else
#define METIS_API(type) type
#endif

#ifdef __cplusplus
extern "C" {
#endif

METIS_API(int)
Highs_METIS_NodeND(idx_t *nvtxs, const idx_t *xadj, const idx_t *adjncy,
                   idx_t *vwgt, idx_t *options, idx_t *perm, idx_t *iperm);

METIS_API(int) Highs_METIS_SetDefaultOptions(idx_t *options);

#ifdef __cplusplus
}
#endif

/*! Return codes */
typedef enum {
  METIS_OK = 1,           /*!< Returned normally */
  METIS_ERROR_INPUT = -2, /*!< Returned due to erroneous inputs and/or options */
  METIS_ERROR_MEMORY = -3, /*!< Returned due to insufficient memory */
  METIS_ERROR = -4         /*!< Some other errors */
} rstatus_et;

/*! Options codes (i.e., options[]) */
typedef enum {
  METIS_OPTION_CTYPE,
  METIS_OPTION_IPTYPE,
  METIS_OPTION_RTYPE,
  METIS_OPTION_DBGLVL,
  METIS_OPTION_NITER,
  METIS_OPTION_SEED,
  METIS_OPTION_COMPRESS,
  METIS_OPTION_CCORDER,
  METIS_OPTION_PFACTOR,
  METIS_OPTION_NSEPS,
  METIS_OPTION_UFACTOR,
  METIS_OPTION_DROPEDGES,
  METIS_OPTION_NO2HOP,
} moptions_et;

#endif /* HIGHS_METIS_H */
