#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL _scipy_sparse_sparsetools_ARRAY_API

#include "sparsetools.h"
#include "dia.h"
#include "csgraph.h"
#include "coo.h"

#define SPTOOLS_OTHER_ROUTINES(X_0, X_I, X_IT) \
  X_IT(coo_tocsr)                              \
  X_IT(coo_todense)                            \
  X_IT(coo_todense_nd)                         \
  X_IT(coo_matvec)                             \
  X_IT(coo_matvec_nd)                          \
  X_IT(coo_matmat_dense)                       \
  X_IT(coo_matmat_dense_nd)                    \
  X_IT(dia_tocsr)                              \
  X_IT(dia_matmat)                             \
  X_IT(dia_matvec)                             \
  X_IT(dia_matvecs)                            \
  X_I(cs_graph_components)

SPTOOLS_DEFINE_UNIT(other, SPTOOLS_OTHER_ROUTINES)

