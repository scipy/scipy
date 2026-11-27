#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL _scipy_sparse_sparsetools_ARRAY_API

#include "sparsetools.h"
#include "csc.h"

#define SPTOOLS_CSC_ROUTINES(X_0, X_I, X_IT) \
  X_IT(csc_matvec)                           \
  X_IT(csc_matvecs)

SPTOOLS_DEFINE_UNIT(csc, SPTOOLS_CSC_ROUTINES)
