#ifndef __SUPERLU_OBJECT /* allow multiple inclusions */
#define __SUPERLU_OBJECT

#include "Python.h"
#define PY_ARRAY_UNIQUE_SYMBOL scipy_superlu
#include "Numeric/arrayobject.h"
#include "SuperLU/SRC/util.h"
#include "SuperLU/SRC/scomplex.h"
#include "SuperLU/SRC/dcomplex.h"

/*********************************************************************** 
 * SuperLUObject definition
 */

typedef struct SciPyLUObject {
  PyObject_VAR_HEAD
  int m,n;
  SuperMatrix L;
  SuperMatrix U;
  int *perm_r;
  int *perm_c;
  int type;
} SciPyLUObject;

extern PyTypeObject SciPySuperLUType;

int DenseSuper_from_Numeric(SuperMatrix *, PyObject *);
int NRFormat_from_spMatrix(SuperMatrix *, int, int, int, PyArrayObject *, PyArrayObject *, PyArrayObject *, int);
int NCFormat_from_spMatrix(SuperMatrix *, int, int, int, PyArrayObject *, PyArrayObject *, PyArrayObject *, int);
colperm_t superlu_module_getpermc(int);
PyObject *newSciPyLUObject(SuperMatrix *, double, double, int, int, int, int);

void
dgstrf (superlu_options_t *, SuperMatrix *, double,
        int, int, int *, void *, int,
        int *, int *, SuperMatrix *, SuperMatrix *,
        SuperLUStat_t *, int *);

void
sgstrf (superlu_options_t *, SuperMatrix *, float,
        int, int, int *, void *, int,
        int *, int *, SuperMatrix *, SuperMatrix *,
        SuperLUStat_t *, int *);

void
cgstrf (superlu_options_t *, SuperMatrix *, float,
        int, int, int *, void *, int,
        int *, int *, SuperMatrix *, SuperMatrix *,
        SuperLUStat_t *, int *);

void
dgstrs (trans_t, SuperMatrix *, SuperMatrix *,
        int *, int *, SuperMatrix *,
        SuperLUStat_t *, int *);

void
sgstrs (trans_t, SuperMatrix *, SuperMatrix *,
        int *, int *, SuperMatrix *,
        SuperLUStat_t *, int *);

void
cgstrs (trans_t, SuperMatrix *, SuperMatrix *,
        int *, int *, SuperMatrix *,
        SuperLUStat_t *, int *);

void
sCreate_Dense_Matrix(SuperMatrix *, int, int, float *, int, Stype_t, Dtype_t, Mtype_t);
void
dCreate_Dense_Matrix(SuperMatrix *, int, int, double *, int, Stype_t, Dtype_t, Mtype_t);
void
cCreate_Dense_Matrix(SuperMatrix *, int, int, complex *, int, Stype_t, Dtype_t, Mtype_t);

void
sCreate_CompRow_Matrix(SuperMatrix *, int, int, int,
                       float *, int *, int *,
                       Stype_t, Dtype_t, Mtype_t);

void
dCreate_CompRow_Matrix(SuperMatrix *, int, int, int,
                       double *, int *, int *,
                       Stype_t, Dtype_t, Mtype_t);

void
cCreate_CompRow_Matrix(SuperMatrix *, int, int, int,
                       complex *, int *, int *,
                       Stype_t, Dtype_t, Mtype_t);

void
sCreate_CompCol_Matrix(SuperMatrix *, int, int, int,
                       float *, int *, int *,
                       Stype_t, Dtype_t, Mtype_t);
void
dCreate_CompCol_Matrix(SuperMatrix *, int, int, int,
                       double *, int *, int *,
                       Stype_t, Dtype_t, Mtype_t);
void
cCreate_CompCol_Matrix(SuperMatrix *, int, int, int,
                       complex *, int *, int *,
                       Stype_t, Dtype_t, Mtype_t);


#endif  /* __SUPERLU_OBJECT */
