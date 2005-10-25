#ifndef SPMATRIX_API_H
#define SPMATRIX_API_H

/*
 * C API
 */

/* Type definitions */

#define LLMatType_NUM 0
#define CSRMatType_NUM 1
#define SSSMatType_NUM 2

/* Function definitions */

#define SpMatrix_ParseVecOpArgs_NUM 3
#define SpMatrix_ParseVecOpArgs_RET int
#define SpMatrix_ParseVecOpArgs_PROTO (PyObject *args, double **x_data, double **y_data, int n)

#define SpMatrix_GetShape_NUM 4
#define SpMatrix_GetShape_RET int
#define SpMatrix_GetShape_PROTO (PyObject *, int[])

#define SpMatrix_GetOrder_NUM 5
#define SpMatrix_GetOrder_RET int
#define SpMatrix_GetOrder_PROTO (PyObject *, int*)

#define SpMatrix_GetItem_NUM 6
#define SpMatrix_GetItem_RET double
#define SpMatrix_GetItem_PROTO (PyObject *, int, int)

#define SpMatrix_Matvec_NUM 7
#define SpMatrix_Matvec_RET int
#define SpMatrix_Matvec_PROTO (PyObject *matrix, int nx, double *x, int ny, double *y)

#define SpMatrix_Precon_NUM 8
#define SpMatrix_Precon_RET int
#define SpMatrix_Precon_PROTO (PyObject *repc, int n, double *x, double *y)

#define SpMatrix_NewLLMatObject_NUM 9
#define SpMatrix_NewLLMatObject_RET PyObject *
#define SpMatrix_NewLLMatObject_PROTO (int dim[], int sym, int sizeHint)

#define SpMatrix_LLMatGetItem_NUM 10
#define SpMatrix_LLMatGetItem_RET double
#define SpMatrix_LLMatGetItem_PROTO (LLMatObject *a, int i, int j)

#define SpMatrix_LLMatSetItem_NUM 11
#define SpMatrix_LLMatSetItem_RET int
#define SpMatrix_LLMatSetItem_PROTO (LLMatObject *a, int i, int j, double x)

#define SpMatrix_LLMatUpdateItemAdd_NUM 12
#define SpMatrix_LLMatUpdateItemAdd_RET int
#define SpMatrix_LLMatUpdateItemAdd_PROTO (LLMatObject *a, int i, int j, double x)

#define SpMatrix_LLMatBuildColIndex_NUM 13
#define SpMatrix_LLMatBuildColIndex_RET int
#define SpMatrix_LLMatBuildColIndex_PROTO (struct llColIndex **idx, LLMatObject *self, int includeDiagonal)

#define SpMatrix_LLMatDestroyColIndex_NUM 14
#define SpMatrix_LLMatDestroyColIndex_RET void
#define SpMatrix_LLMatDestroyColIndex_PROTO (struct llColIndex **idx)

#define ItSolvers_Solve_NUM 15
#define ItSolvers_Solve_RET int
#define ItSolvers_Solve_PROTO (PyObject *linsolver, PyObject *A, int n, \
                double *b, double *x, double tol, int itmax, PyObject *K, \
                int *info, int *iter, double *relres)


/* Total number of C API pointers */
#define SpMatrix_API_pointers 16 

#ifdef SPMATRIX_MODULE

static SpMatrix_ParseVecOpArgs_RET SpMatrix_ParseVecOpArgs SpMatrix_ParseVecOpArgs_PROTO;
static SpMatrix_GetShape_RET SpMatrix_GetShape SpMatrix_GetShape_PROTO;
static SpMatrix_GetOrder_RET SpMatrix_GetOrder SpMatrix_GetOrder_PROTO;
static SpMatrix_GetItem_RET SpMatrix_GetItem SpMatrix_GetItem_PROTO;
static SpMatrix_Matvec_RET SpMatrix_Matvec SpMatrix_Matvec_PROTO;
static SpMatrix_Precon_RET SpMatrix_Precon SpMatrix_Precon_PROTO;
static SpMatrix_NewLLMatObject_RET SpMatrix_NewLLMatObject SpMatrix_NewLLMatObject_PROTO;
static SpMatrix_LLMatGetItem_RET SpMatrix_LLMatGetItem SpMatrix_LLMatGetItem_PROTO;
static SpMatrix_LLMatSetItem_RET SpMatrix_LLMatSetItem SpMatrix_LLMatSetItem_PROTO;
static SpMatrix_LLMatUpdateItemAdd_RET SpMatrix_LLMatUpdateItemAdd SpMatrix_LLMatUpdateItemAdd_PROTO;
static SpMatrix_LLMatBuildColIndex_RET SpMatrix_LLMatBuildColIndex SpMatrix_LLMatBuildColIndex_PROTO;
static SpMatrix_LLMatDestroyColIndex_RET SpMatrix_LLMatDestroyColIndex SpMatrix_LLMatDestroyColIndex_PROTO;
static ItSolvers_Solve_RET ItSolvers_Solve ItSolvers_Solve_PROTO;

#define init_c_api(MODULE_DICT) \
{ \
  static void *SpMatrix_API[SpMatrix_API_pointers]; \
  PyObject *c_api_object; \
\
  /* initialize C API pointer array */ \
  SpMatrix_API[LLMatType_NUM] = (void *) &LLMatType; \
  SpMatrix_API[CSRMatType_NUM] = (void *) &CSRMatType; \
  SpMatrix_API[SSSMatType_NUM] = (void *) &SSSMatType; \
  SpMatrix_API[SpMatrix_ParseVecOpArgs_NUM] = (void *) SpMatrix_ParseVecOpArgs; \
  SpMatrix_API[SpMatrix_GetShape_NUM] = (void *) SpMatrix_GetShape; \
  SpMatrix_API[SpMatrix_GetOrder_NUM] = (void *) SpMatrix_GetOrder; \
  SpMatrix_API[SpMatrix_GetItem_NUM] = (void *) SpMatrix_GetItem; \
  SpMatrix_API[SpMatrix_Matvec_NUM] = (void *) SpMatrix_Matvec; \
  SpMatrix_API[SpMatrix_Precon_NUM] = (void *) SpMatrix_Precon; \
  SpMatrix_API[SpMatrix_NewLLMatObject_NUM] = (void *) SpMatrix_NewLLMatObject; \
  SpMatrix_API[SpMatrix_LLMatGetItem_NUM] = (void *) SpMatrix_LLMatGetItem; \
  SpMatrix_API[SpMatrix_LLMatSetItem_NUM] = (void *) SpMatrix_LLMatSetItem; \
  SpMatrix_API[SpMatrix_LLMatUpdateItemAdd_NUM] = (void *) SpMatrix_LLMatUpdateItemAdd; \
  SpMatrix_API[SpMatrix_LLMatBuildColIndex_NUM] = (void *) SpMatrix_LLMatBuildColIndex; \
  SpMatrix_API[SpMatrix_LLMatDestroyColIndex_NUM] = (void *) SpMatrix_LLMatDestroyColIndex; \
  SpMatrix_API[ItSolvers_Solve_NUM] = (void *) ItSolvers_Solve; \
\
  /* Create a CObject containing the API pointer array s address */ \
  c_api_object = PyCObject_FromVoidPtr((void *)SpMatrix_API, NULL); \
\
  /* create a name for this object in the module's namespace */ \
  if (c_api_object != NULL) { \
    PyDict_SetItemString((MODULE_DICT), "_C_API", c_api_object); \
    Py_DECREF(c_api_object); \
  } \
}

#else

#ifdef SPMATRIX_UNIQUE_SYMBOL
#define SpMatrix_API SPMATRIX_UNIQUE_SYMBOL
#endif

/* C API address pointer */
#ifdef NO_IMPORT_SPMATRIX
extern void **SpMatrix_API;
#else
#ifdef SPMATRIX_UNIQUE_SYMBOL
void **SpMatrix_API;
#else
static void **SpMatrix_API;
#endif
#endif

#define LLMatType *(PyTypeObject *)SpMatrix_API[LLMatType_NUM]
#define CSRMatType *(PyTypeObject *)SpMatrix_API[CSRMatType_NUM]
#define SSSMatType *(PyTypeObject *)SpMatrix_API[SSSMatType_NUM]

#define SpMatrix_ParseVecOpArgs \
  (*(SpMatrix_ParseVecOpArgs_RET (*)SpMatrix_ParseVecOpArgs_PROTO) \
  SpMatrix_API[SpMatrix_ParseVecOpArgs_NUM])
  
#define SpMatrix_GetShape \
  (*(SpMatrix_GetShape_RET (*)SpMatrix_GetShape_PROTO) \
  SpMatrix_API[SpMatrix_GetShape_NUM])
  
#define SpMatrix_GetOrder \
  (*(SpMatrix_GetOrder_RET (*)SpMatrix_GetOrder_PROTO) \
  SpMatrix_API[SpMatrix_GetOrder_NUM])
  
#define SpMatrix_GetItem \
  (*(SpMatrix_GetItem_RET (*)SpMatrix_GetItem_PROTO) \
  SpMatrix_API[SpMatrix_GetItem_NUM])
  
#define SpMatrix_Matvec \
  (*(SpMatrix_Matvec_RET (*)SpMatrix_Matvec_PROTO) \
  SpMatrix_API[SpMatrix_Matvec_NUM])
  
#define SpMatrix_Precon \
  (*(SpMatrix_Precon_RET (*)SpMatrix_Precon_PROTO) \
  SpMatrix_API[SpMatrix_Precon_NUM])
  
#define SpMatrix_NewLLMatObject \
  (*(SpMatrix_NewLLMatObject_RET (*)SpMatrix_NewLLMatObject_PROTO) \
  SpMatrix_API[SpMatrix_NewLLMatObject_NUM])
  
#define SpMatrix_LLMatGetItem \
  (*(SpMatrix_LLMatGetItem_RET (*)SpMatrix_LLMatGetItem_PROTO) \
  SpMatrix_API[SpMatrix_LLMatGetItem_NUM])
  
#define SpMatrix_LLMatSetItem \
  (*(SpMatrix_LLMatSetItem_RET (*)SpMatrix_LLMatSetItem_PROTO) \
  SpMatrix_API[SpMatrix_LLMatSetItem_NUM])
  
#define SpMatrix_LLMatUpdateItemAdd \
  (*(SpMatrix_LLMatUpdateItemAdd_RET (*)SpMatrix_LLMatUpdateItemAdd_PROTO) \
  SpMatrix_API[SpMatrix_LLMatUpdateItemAdd_NUM])
  
#define SpMatrix_LLMatBuildColIndex \
  (*(SpMatrix_LLMatBuildColIndex_RET (*)SpMatrix_LLMatBuildColIndex_PROTO) \
  SpMatrix_API[SpMatrix_LLMatBuildColIndex_NUM])
  
#define SpMatrix_LLMatDestroyColIndex \
  (*(SpMatrix_LLMatDestroyColIndex_RET (*)SpMatrix_LLMatDestroyColIndex_PROTO) \
  SpMatrix_API[SpMatrix_LLMatDestroyColIndex_NUM])
  
#define ItSolvers_Solve \
  (*(ItSolvers_Solve_RET (*)ItSolvers_Solve_PROTO) \
  SpMatrix_API[ItSolvers_Solve_NUM])
  
#define import_spmatrix() \
{ \
  PyObject *spmatrix = PyImport_ImportModule("spmatrix"); \
  if (spmatrix != NULL) { \
    PyObject *module_dict = PyModule_GetDict(spmatrix); \
    PyObject *c_api_object = PyDict_GetItemString(module_dict, "_C_API"); \
    if (PyCObject_Check(c_api_object)) { \
      SpMatrix_API = (void **)PyCObject_AsVoidPtr(c_api_object); \
    } \
  } \
  assert(SpMatrix_API != NULL); \
}

#endif

#endif

