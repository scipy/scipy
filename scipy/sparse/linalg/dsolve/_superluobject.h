/* -*-c-*-  */
/*
 * _superlu object
 *
 * Python object representing SuperLU factorization + some utility functions.
 */

#ifndef __SUPERLU_OBJECT
#define __SUPERLU_OBJECT

#include <setjmp.h>
#include <Python.h>

/* Undef a macro from Python which conflicts with superlu */
#ifdef c_abs
#undef c_abs
#endif

#include "SuperLU/SRC/slu_zdefs.h"
#include "numpy/arrayobject.h"
#include "SuperLU/SRC/slu_util.h"
#include "SuperLU/SRC/slu_dcomplex.h"
#include "SuperLU/SRC/slu_scomplex.h"


#define _CHECK_INTEGER(x) (PyArray_ISINTEGER(x) && (x)->descr->elsize == sizeof(int))

/*
 * SuperLUObject definition
 */
typedef struct {
    PyObject_HEAD
    npy_intp m, n;
    SuperMatrix L;
    SuperMatrix U;
    int *perm_r;
    int *perm_c;
    PyObject *cached_U;
    PyObject *cached_L;
    int type;
} SuperLUObject;

typedef struct {
    PyObject_HEAD
    jmp_buf jmpbuf;
    PyObject *memory_dict;
} SuperLUGlobalObject;

extern PyTypeObject SuperLUType;
extern PyTypeObject SuperLUGlobalType;

int DenseSuper_from_Numeric(SuperMatrix *, PyObject *);
int NRFormat_from_spMatrix(SuperMatrix *, int, int, int, PyArrayObject *,
			   PyArrayObject *, PyArrayObject *, int);
int NCFormat_from_spMatrix(SuperMatrix *, int, int, int, PyArrayObject *,
			   PyArrayObject *, PyArrayObject *, int);
int LU_to_csc_matrix(SuperMatrix *L, SuperMatrix *U,
                     PyObject **L_csc, PyObject **U_csc);
colperm_t superlu_module_getpermc(int);
PyObject *newSuperLUObject(SuperMatrix *, PyObject *, int, int);
int set_superlu_options_from_dict(superlu_options_t * options,
				  int ilu, PyObject * option_dict,
				  int *panel_size, int *relax);

void XDestroy_SuperMatrix_Store(SuperMatrix *);
void XDestroy_SuperNode_Matrix(SuperMatrix *);
void XDestroy_CompCol_Matrix(SuperMatrix *);
void XDestroy_CompCol_Permuted(SuperMatrix *);
void XStatFree(SuperLUStat_t *);

jmp_buf *superlu_python_jmpbuf();

/*
 * Definitions for other SuperLU data types than Z,
 * and type-generic definitions.
 */

#define CHECK_SLU_TYPE(type) \
    (type == NPY_FLOAT || type == NPY_DOUBLE || type == NPY_CFLOAT || type == NPY_CDOUBLE)

#define TYPE_GENERIC_FUNC(name, returntype)                \
    returntype s##name(name##_ARGS);                       \
    returntype d##name(name##_ARGS);                       \
    returntype c##name(name##_ARGS);                       \
    static returntype name(int type, name##_ARGS)          \
    {                                                      \
        switch(type) {                                     \
        case NPY_FLOAT:   s##name(name##_ARGS_REF); break; \
        case NPY_DOUBLE:  d##name(name##_ARGS_REF); break; \
        case NPY_CFLOAT:  c##name(name##_ARGS_REF); break; \
        case NPY_CDOUBLE: z##name(name##_ARGS_REF); break; \
        default: return;                                   \
        }                                                  \
    }

#define SLU_TYPECODE_TO_NPY(s)                    \
    ( ((s) == SLU_S) ? NPY_FLOAT :                \
      ((s) == SLU_D) ? NPY_DOUBLE :               \
      ((s) == SLU_C) ? NPY_CFLOAT :               \
      ((s) == SLU_Z) ? NPY_CDOUBLE : -1)

#define NPY_TYPECODE_TO_SLU(s)                    \
    ( ((s) == NPY_FLOAT) ? SLU_S :                \
      ((s) == NPY_DOUBLE) ? SLU_D :               \
      ((s) == NPY_CFLOAT) ? SLU_C :               \
      ((s) == NPY_CDOUBLE) ? SLU_Z : -1)

#define gstrf_ARGS                                                  \
    superlu_options_t *a, SuperMatrix *b,                           \
    int c, int d, int *e, void *f, int g,                           \
    int *h, int *i, SuperMatrix *j, SuperMatrix *k,                 \
    GlobalLU_t *l, SuperLUStat_t *m, int *n
#define gstrf_ARGS_REF a,b,c,d,e,f,g,h,i,j,k,l,m,n

#define gsitrf_ARGS gstrf_ARGS
#define gsitrf_ARGS_REF gstrf_ARGS_REF

#define gstrs_ARGS                              \
    trans_t a, SuperMatrix *b, SuperMatrix *c,  \
    int *d, int *e, SuperMatrix *f,             \
    SuperLUStat_t *g, int *h
#define gstrs_ARGS_REF a,b,c,d,e,f,g,h

#define gssv_ARGS                                               \
    superlu_options_t *a, SuperMatrix *b, int *c, int *d,       \
    SuperMatrix *e, SuperMatrix *f, SuperMatrix *g,             \
    SuperLUStat_t *h, int *i
#define gssv_ARGS_REF a,b,c,d,e,f,g,h,i

#define Create_Dense_Matrix_ARGS                               \
    SuperMatrix *a, int b, int c, void *d, int e,              \
    Stype_t f, Dtype_t g, Mtype_t h
#define Create_Dense_Matrix_ARGS_REF a,b,c,d,e,f,g,h

#define Create_CompRow_Matrix_ARGS              \
    SuperMatrix *a, int b, int c, int d,        \
    void *e, int *f, int *g,                    \
    Stype_t h, Dtype_t i, Mtype_t j
#define Create_CompRow_Matrix_ARGS_REF a,b,c,d,e,f,g,h,i,j

#define Create_CompCol_Matrix_ARGS Create_CompRow_Matrix_ARGS
#define Create_CompCol_Matrix_ARGS_REF Create_CompRow_Matrix_ARGS_REF

TYPE_GENERIC_FUNC(gstrf, void);
TYPE_GENERIC_FUNC(gsitrf, void);
TYPE_GENERIC_FUNC(gstrs, void);
TYPE_GENERIC_FUNC(gssv, void);
TYPE_GENERIC_FUNC(Create_Dense_Matrix, void);
TYPE_GENERIC_FUNC(Create_CompRow_Matrix, void);
TYPE_GENERIC_FUNC(Create_CompCol_Matrix, void);

#endif				/* __SUPERLU_OBJECT */
