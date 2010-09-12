/* -*-c-*-  */
/*
 * _superlu object
 *
 * Python object representing SuperLU factorization + some utility functions.
 */

#include <Python.h>

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL _scipy_sparse_superlu_ARRAY_API

#include "_superluobject.h"
#include "numpy/npy_3kcompat.h"
#include <setjmp.h>
#include <ctype.h>

extern jmp_buf _superlu_py_jmpbuf;


/*********************************************************************** 
 * SciPyLUObject methods
 */

static char solve_doc[] = "x = self.solve(b, trans)\n\
\n\
solves linear system of equations with one or sereral right hand sides.\n\
\n\
parameters\n\
----------\n\
\n\
b        array, right hand side(s) of equation\n\
x        array, solution vector(s)\n\
trans    'N': solve A   * x == b\n\
         'T': solve A^T * x == b\n\
         'H': solve A^H * x == b\n\
         (optional, default value 'N')\n\
";

static PyObject *
SciPyLU_solve(SciPyLUObject *self, PyObject *args, PyObject *kwds) {
  PyArrayObject *b, *x=NULL;
  SuperMatrix B;
#ifndef NPY_PY3K
  char itrans = 'N';
#else
  int itrans = 'N';
#endif
  int info;
  trans_t trans;
  SuperLUStat_t stat;

  static char *kwlist[] = {"rhs","trans",NULL};

  if (!CHECK_SLU_TYPE(self->type)) {
      PyErr_SetString(PyExc_ValueError, "unsupported data type");
      return NULL;
  }

#ifndef NPY_PY3K
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|c", kwlist,
#else
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|C", kwlist,
#endif
                                   &PyArray_Type, &b, 
                                   &itrans))
    return NULL;

  /* solve transposed system: matrix was passed row-wise instead of
   * column-wise */
  if (itrans == 'n' || itrans == 'N')
      trans = NOTRANS;
  else if (itrans == 't' || itrans == 'T')
      trans = TRANS;
  else if (itrans == 'h' || itrans == 'H')
      trans = CONJ;
  else {
    PyErr_SetString(PyExc_ValueError, "trans must be N, T, or H");
    return NULL;
  }

  if ((x = (PyArrayObject *) \
       PyArray_CopyFromObject((PyObject *)b,self->type,1,2))==NULL) return NULL;

  if (b->dimensions[0] != self->n) goto fail;


  if (setjmp(_superlu_py_jmpbuf)) goto fail; 

  if (DenseSuper_from_Numeric(&B, (PyObject *)x)) goto fail;

  StatInit(&stat);

  /* Solve the system, overwriting vector x. */
  gstrs(self->type,
        trans, &self->L, &self->U, self->perm_c, self->perm_r, &B,
        &stat, &info);

  if (info) { 
      PyErr_SetString(PyExc_SystemError,
                      "gstrs was called with invalid arguments");
      goto fail;
  }
  
  /* free memory */
  Destroy_SuperMatrix_Store(&B);
  StatFree(&stat);
  return (PyObject *)x;

fail:
  Destroy_SuperMatrix_Store(&B);  
  StatFree(&stat);
  Py_XDECREF(x);
  return NULL;
}

/** table of object methods
 */
PyMethodDef SciPyLU_methods[] = {
  {"solve", (PyCFunction)SciPyLU_solve, METH_VARARGS|METH_KEYWORDS, solve_doc},
  {NULL, NULL}			/* sentinel */
};


/*********************************************************************** 
 * SciPySuperLUType methods
 */

static void
SciPyLU_dealloc(SciPyLUObject *self)
{
  SUPERLU_FREE(self->perm_r);
  SUPERLU_FREE(self->perm_c);
  if (self->L.Store != NULL) {
      Destroy_SuperNode_Matrix(&self->L);
  }
  if (self->U.Store != NULL) {
      Destroy_CompCol_Matrix(&self->U);
  }
  PyObject_Del(self);
}

static PyObject *
SciPyLU_getattr(SciPyLUObject *self, char *name)
{
  if (strcmp(name, "shape") == 0)
    return Py_BuildValue("(i,i)", self->m, self->n);
  if (strcmp(name, "nnz") == 0)
    return Py_BuildValue("i", ((SCformat *)self->L.Store)->nnz + ((SCformat *)self->U.Store)->nnz);
  if (strcmp(name, "perm_r") == 0) {
    PyArrayObject* perm_r = PyArray_SimpleNewFromData(1, (npy_intp*) (&self->n), NPY_INT, (void*)self->perm_r);
    /* For ref counting of the memory */
    PyArray_BASE(perm_r) = self;
    Py_INCREF(self);
    return perm_r ;
  }
  if (strcmp(name, "perm_c") == 0) {
    PyArrayObject* perm_c = PyArray_SimpleNewFromData(1, (npy_intp*) (&self->n), NPY_INT, (void*)self->perm_c);
    /* For ref counting of the memory */
    PyArray_BASE(perm_c) = self;
    Py_INCREF(self);
    return perm_c ;
  }
  if (strcmp(name, "__members__") == 0) {
    char *members[] = {"shape", "nnz", "perm_r", "perm_c"};
    int i;

    PyObject *list = PyList_New(sizeof(members)/sizeof(char *));
    if (list != NULL) {
      for (i = 0; i < sizeof(members)/sizeof(char *); i ++)
	PyList_SetItem(list, i, PyUString_FromString(members[i]));
      if (PyErr_Occurred()) {
	Py_DECREF(list);
	list = NULL;
      }
    }
    return list;
  }
#if PY_VERSION_HEX >= 0x03000000
  if (1) {
      PyObject *str, *ret;
      str = PyUnicode_FromString(name);
      ret = PyObject_GenericGetAttr((PyObject *)self, str);
      Py_DECREF(str);
      return ret;
  }
#else
  return Py_FindMethod(SciPyLU_methods, (PyObject *)self, name);
#endif
}


/***********************************************************************
 * SciPySuperLUType structure
 */
static char factored_lu_doc[] = "\
Object resulting from a factorization of a sparse matrix\n\
\n\
Attributes\n\
-----------\n\
\n\
shape : 2-tuple\n\
    the shape of the orginal matrix factored\n\
nnz : int\n\
    the number of non-zero elements in the matrix\n\
perm_c\n\
    the permutation applied to the colums of the matrix for the LU factorization\n\
perm_r\n\
    the permutation applied to the rows of the matrix for the LU factorization\n\
\n\
Methods\n\
-------\n\
solve\n\
    solves the system for a given right hand side vector\n \
\n\
";

PyTypeObject SciPySuperLUType = {
#if defined(NPY_PY3K)
  PyVarObject_HEAD_INIT(NULL, 0)
#else
  PyObject_HEAD_INIT(NULL)
  0,
#endif
  "factored_lu",
  sizeof(SciPyLUObject),
  0,
  (destructor)SciPyLU_dealloc,   /* tp_dealloc */
  0,				/* tp_print */
  (getattrfunc)SciPyLU_getattr,  /* tp_getattr */
  0,				/* tp_setattr */
  0,				/* tp_compare / tp_reserved */
  0,				/* tp_repr */
  0,				/* tp_as_number*/
  0,				/* tp_as_sequence*/
  0,				/* tp_as_mapping*/
  0,				/* tp_hash */
  0,				/* tp_call */
  0,				/* tp_str */
  0,				/* tp_getattro */
  0,				/* tp_setattro */
  0,				/* tp_as_buffer */
  Py_TPFLAGS_DEFAULT,		/* tp_flags */
  factored_lu_doc,		/* tp_doc */
  0,                            /* tp_traverse */
  0,                            /* tp_clear */
  0,                            /* tp_richcompare */
  0,                            /* tp_weaklistoffset */
  0,                            /* tp_iter */
  0,                            /* tp_iternext */
  SciPyLU_methods,              /* tp_methods */
  0,                            /* tp_members */
  0,                            /* tp_getset */
  0,                            /* tp_base */
  0,                            /* tp_dict */
  0,                            /* tp_descr_get */
  0,                            /* tp_descr_set */
  0,                            /* tp_dictoffset */
  0,                            /* tp_init */
  0,                            /* tp_alloc */
  0,                            /* tp_new */
  0,                            /* tp_free */
  0,                            /* tp_is_gc */
  0,                            /* tp_bases */
  0,                            /* tp_mro */
  0,                            /* tp_cache */
  0,                            /* tp_subclasses */
  0,                            /* tp_weaklist */
  0,                            /* tp_del */
#if PY_VERSION_HEX >= 0x02060000
  0,                            /* tp_version_tag */
#endif
};


int DenseSuper_from_Numeric(SuperMatrix *X, PyObject *PyX)
{
  int m, n, ldx, nd;
  PyArrayObject *aX;
  
  if (!PyArray_Check(PyX)) {
    PyErr_SetString(PyExc_TypeError, "dgssv: Second argument is not an array.");
    return -1;
  }

  aX = (PyArrayObject *)PyX;
  nd = aX->nd;

  if (nd == 1) {
    m = aX->dimensions[0];
    n = 1;
    ldx = m;
  }
  else {  /* nd == 2 */
    m = aX->dimensions[1];
    n = aX->dimensions[0];
    ldx = m;
  }
  
  if (setjmp(_superlu_py_jmpbuf))
      return -1;
  else {
      if (!CHECK_SLU_TYPE(aX->descr->type_num)) {
          PyErr_SetString(PyExc_ValueError, "unsupported data type");
          return -1;
      }
      Create_Dense_Matrix(aX->descr->type_num, X, m, n,
                          aX->data, ldx, SLU_DN,
                          NPY_TYPECODE_TO_SLU(aX->descr->type_num), SLU_GE);
  }
  return 0;
}

/* Natively handles Compressed Sparse Row and CSC */

int NRFormat_from_spMatrix(SuperMatrix *A, int m, int n, int nnz,
                           PyArrayObject *nzvals, PyArrayObject *colind,
                           PyArrayObject *rowptr, int typenum)
{
  int err = 0;
    
  err = (nzvals->descr->type_num != typenum);
  err += (nzvals->nd != 1);
  err += (nnz > nzvals->dimensions[0]);
  if (err) {
    PyErr_SetString(PyExc_TypeError, "Fourth argument must be a 1-D array at least as big as third argument.");
    return -1;
  }

  if (setjmp(_superlu_py_jmpbuf))
      return -1;
  else {
      if (!CHECK_SLU_TYPE(nzvals->descr->type_num)) {
          PyErr_SetString(PyExc_TypeError, "Invalid type for array.");
          return -1;  
      }
      Create_CompRow_Matrix(nzvals->descr->type_num,
                            A, m, n, nnz, nzvals->data, (int *)colind->data,
                            (int *)rowptr->data, SLU_NR,
                            NPY_TYPECODE_TO_SLU(nzvals->descr->type_num),
                            SLU_GE);
  }

  return 0;
}

int NCFormat_from_spMatrix(SuperMatrix *A, int m, int n, int nnz,
                           PyArrayObject *nzvals, PyArrayObject *rowind,
                           PyArrayObject *colptr, int typenum)
{
  int err=0;

  err = (nzvals->descr->type_num != typenum);
  err += (nzvals->nd != 1);
  err += (nnz > nzvals->dimensions[0]);
  if (err) {
    PyErr_SetString(PyExc_TypeError, "Fifth argument must be a 1-D array at least as big as fourth argument.");
    return -1;
  }


  if (setjmp(_superlu_py_jmpbuf))
      return -1;
  else {
      if (!CHECK_SLU_TYPE(nzvals->descr->type_num)) {
          PyErr_SetString(PyExc_TypeError, "Invalid type for array.");
          return -1;  
      }
      Create_CompCol_Matrix(nzvals->descr->type_num,
                            A, m, n, nnz, nzvals->data, (int *)rowind->data,
                            (int *)colptr->data, SLU_NC,
                            NPY_TYPECODE_TO_SLU(nzvals->descr->type_num),
                            SLU_GE);
  }

  return 0;
}

PyObject *
newSciPyLUObject(SuperMatrix *A, PyObject *option_dict, int intype, int ilu)
{

   /* A must be in SLU_NC format used by the factorization routine. */
  SciPyLUObject *self;
  SuperMatrix AC;     /* Matrix postmultiplied by Pc */
  int lwork = 0;
  int *etree=NULL;
  int info;
  int n;
  superlu_options_t options;
  SuperLUStat_t stat;
  int panel_size, relax;
  int trf_finished = 0;

  n = A->ncol;

  if (!set_superlu_options_from_dict(&options, ilu, option_dict,
                                     &panel_size, &relax)) {
      return NULL;
  }

  /* Create SciPyLUObject */
  self = PyObject_New(SciPyLUObject, &SciPySuperLUType);
  if (self == NULL)
    return PyErr_NoMemory();
  self->m = A->nrow;
  self->n = n;
  self->perm_r = NULL;
  self->perm_c = NULL;
  self->type = intype;

  if (setjmp(_superlu_py_jmpbuf)) goto fail;
  
  /* Calculate and apply minimum degree ordering*/
  etree = intMalloc(n);
  self->perm_r = intMalloc(n);
  self->perm_c = intMalloc(n);
  StatInit(&stat);

  get_perm_c(options.ColPerm, A, self->perm_c); /* calc column permutation */
  sp_preorder(&options, A, self->perm_c, etree, &AC); /* apply column
                                                       * permutation */

  /* Perform factorization */
  if (!CHECK_SLU_TYPE(SLU_TYPECODE_TO_NPY(A->Dtype))) {
      PyErr_SetString(PyExc_ValueError, "Invalid type in SuperMatrix.");
      goto fail;
  }
  if (ilu) {
      gsitrf(SLU_TYPECODE_TO_NPY(A->Dtype),
             &options, &AC, relax, panel_size,
             etree, NULL, lwork, self->perm_c, self->perm_r,
             &self->L, &self->U, &stat, &info);
  }
  else {
      gstrf(SLU_TYPECODE_TO_NPY(A->Dtype),
            &options, &AC, relax, panel_size,
            etree, NULL, lwork, self->perm_c, self->perm_r,
            &self->L, &self->U, &stat, &info);
  }
  trf_finished = 1;

  if (info) {
    if (info < 0)
        PyErr_SetString(PyExc_SystemError,
                        "gstrf was called with invalid arguments");
    else {
        if (info <= n) 
            PyErr_SetString(PyExc_RuntimeError, "Factor is exactly singular");
        else
            PyErr_NoMemory();
    }
    goto fail;
  }

  /* free memory */
  SUPERLU_FREE(etree);
  Destroy_CompCol_Permuted(&AC);
  StatFree(&stat);
  
  return (PyObject *)self;

fail:
  if (!trf_finished) {
      /* Avoid trying to free partially initialized matrices;
         might leak some memory, but avoids a crash */
      self->L.Store = NULL;
      self->U.Store = NULL;
  }
  SUPERLU_FREE(etree);
  Destroy_CompCol_Permuted(&AC);
  StatFree(&stat);
  SciPyLU_dealloc(self);
  return NULL;
}


/***********************************************************************
 * Preparing superlu_options_t
 */

#define ENUM_CHECK_INIT                         \
    long i = -1;                                \
    char *s = "";                               \
    PyObject *tmpobj = NULL;                    \
    if (input == Py_None) return 1;             \
    if (PyString_Check(input)) {                \
        s = PyString_AS_STRING(input);          \
    }                                           \
    else if (PyUnicode_Check(input)) {          \
        tmpobj = PyUnicode_AsASCIIString(input);\
        if (tmpobj == NULL) return 0;           \
        s = PyString_AS_STRING(tmpobj);         \
    }                                           \
    else if (PyInt_Check(input)) {              \
        i = PyInt_AsLong(input);                \
    }

#define ENUM_CHECK_FINISH(message)              \
    Py_XDECREF(tmpobj);                         \
    PyErr_SetString(PyExc_ValueError, message); \
    return 0;

#define ENUM_CHECK(name)                                \
    if (my_strxcmp(s, #name) == 0 || i == (long)name) { \
        *value = name;                                  \
        Py_XDECREF(tmpobj);                             \
        return 1;                                       \
    }

/*
 * Compare strings ignoring case, underscores and whitespace
 */
static int my_strxcmp(const char *a, const char *b)
{
    int c;
    while (*a != '\0' && *b != '\0') {
        while (*a == '_' || isspace(*a)) ++a;
        while (*b == '_' || isspace(*b)) ++b;
        c = (int)tolower(*a) - (int)tolower(*b);
        if (c != 0) {
            return c;
        }
        ++a;
        ++b;
    }
    return (int)tolower(*a) - (int)tolower(*b);
}

static int yes_no_cvt(PyObject *input, yes_no_t *value)
{
    if (input == Py_None) {
        return 1;
    }
    else if (input == Py_True) {
        *value = YES;
    } else if (input == Py_False) {
        *value = NO;
    } else {
        PyErr_SetString(PyExc_ValueError, "value not a boolean");
        return 0;
    }
    return 1;
}

static int fact_cvt(PyObject *input, fact_t *value)
{
    ENUM_CHECK_INIT;
    ENUM_CHECK(DOFACT);
    ENUM_CHECK(SamePattern);
    ENUM_CHECK(SamePattern_SameRowPerm);
    ENUM_CHECK(FACTORED);
    ENUM_CHECK_FINISH("invalid value for 'Fact' parameter");
}

static int rowperm_cvt(PyObject *input, rowperm_t *value)
{
    ENUM_CHECK_INIT;
    ENUM_CHECK(NOROWPERM);
    ENUM_CHECK(LargeDiag);
    ENUM_CHECK(MY_PERMR);
    ENUM_CHECK_FINISH("invalid value for 'RowPerm' parameter");
}

static int colperm_cvt(PyObject *input, colperm_t *value)
{
    ENUM_CHECK_INIT;
    ENUM_CHECK(NATURAL);
    ENUM_CHECK(MMD_ATA);
    ENUM_CHECK(MMD_AT_PLUS_A);
    ENUM_CHECK(COLAMD);
    ENUM_CHECK(MY_PERMC);
    ENUM_CHECK_FINISH("invalid value for 'ColPerm' parameter");
}

static int trans_cvt(PyObject *input, trans_t *value)
{
    ENUM_CHECK_INIT;
    ENUM_CHECK(NOTRANS);
    ENUM_CHECK(TRANS);
    ENUM_CHECK(CONJ);
    if (my_strxcmp(s, "N") == 0) { *value = NOTRANS; return 1; }
    if (my_strxcmp(s, "T") == 0) { *value = TRANS; return 1; }
    if (my_strxcmp(s, "H") == 0) { *value = CONJ; return 1; }
    ENUM_CHECK_FINISH("invalid value for 'Trans' parameter");
}

static int iterrefine_cvt(PyObject *input, IterRefine_t *value)
{
    ENUM_CHECK_INIT;
    ENUM_CHECK(NOREFINE);
    ENUM_CHECK(SINGLE);
    ENUM_CHECK(DOUBLE);
    ENUM_CHECK(EXTRA);
    ENUM_CHECK_FINISH("invalid value for 'IterRefine' parameter");
}

static int norm_cvt(PyObject *input, norm_t *value)
{
    ENUM_CHECK_INIT;
    ENUM_CHECK(ONE_NORM);
    ENUM_CHECK(TWO_NORM);
    ENUM_CHECK(INF_NORM);
    ENUM_CHECK_FINISH("invalid value for 'ILU_Norm' parameter");
}

static int milu_cvt(PyObject *input, milu_t *value)
{
    ENUM_CHECK_INIT;
    ENUM_CHECK(SILU);
    ENUM_CHECK(SMILU_1);
    ENUM_CHECK(SMILU_2);
    ENUM_CHECK(SMILU_3);
    ENUM_CHECK_FINISH("invalid value for 'ILU_MILU' parameter");
}

static int droprule_one_cvt(PyObject *input, int *value)
{
    ENUM_CHECK_INIT;
    if (my_strxcmp(s, "BASIC") == 0) { *value = DROP_BASIC; return 1; }
    if (my_strxcmp(s, "PROWS") == 0) { *value = DROP_PROWS; return 1; }
    if (my_strxcmp(s, "COLUMN") == 0) { *value = DROP_COLUMN; return 1; }
    if (my_strxcmp(s, "AREA") == 0) { *value = DROP_AREA; return 1; }
    if (my_strxcmp(s, "SECONDARY") == 0) { *value = DROP_SECONDARY; return 1; }
    if (my_strxcmp(s, "DYNAMIC") == 0) { *value = DROP_DYNAMIC; return 1; }
    if (my_strxcmp(s, "INTERP") == 0) { *value = DROP_INTERP; return 1; }
    ENUM_CHECK_FINISH("invalid value for 'ILU_DropRule' parameter");
}

static int droprule_cvt(PyObject *input, int *value)
{
    PyObject *seq = NULL;
    int i;
    int rule = 0;

    if (input == Py_None) {
        /* Leave as default */
        return 1;
    }
    else if (PyInt_Check(input)) {
        *value = PyInt_AsLong(input);
        return 1;
    }
    else if (PyString_Check(input)) {
        /* Comma-separated string */
        seq = PyObject_CallMethod(input, "split", "s", ",");
        if (seq == NULL || !PySequence_Check(seq))
            goto fail;
    }
    else if (PyUnicode_Check(input)) {
        /* Comma-separated string */
        PyObject *s;
        int ret;
        s = PyUnicode_AsASCIIString(input);
        if (s == NULL) {
            goto fail;
        }
        ret = droprule_cvt(s, value);
        Py_DECREF(s);
        return ret;
    }
    else if (PySequence_Check(input)) {
        /* Sequence of strings or integers */
        seq = input;
        Py_INCREF(seq);
    }
    else {
        PyErr_SetString(PyExc_ValueError, "invalid value for drop rule");
        goto fail;
    }

    /* OR multiple values together */
    for (i = 0; i < PySequence_Size(seq); ++i) {
        PyObject *item;
        int one_value;
        item = PySequence_ITEM(seq, i);
        if (item == NULL) {
            goto fail;
        }
        if (!droprule_one_cvt(item, &one_value)) {
            Py_DECREF(item);
            goto fail;
        }
        Py_DECREF(item);
        rule |= one_value;
    }
    Py_DECREF(seq);

    *value = rule;
    return 1;

fail:
    Py_XDECREF(seq);
    return 0;
}

static int double_cvt(PyObject *input, double *value)
{
    if (input == Py_None) return 1;
    *value = PyFloat_AsDouble(input);
    if (PyErr_Occurred()) return 0;
    return 1;
}

static int int_cvt(PyObject *input, int *value)
{
    if (input == Py_None) return 1;
    *value = PyInt_AsLong(input);
    if (PyErr_Occurred()) return 0;
    return 1;
}

int set_superlu_options_from_dict(superlu_options_t *options,
                                  int ilu, PyObject *option_dict,
                                  int *panel_size, int *relax)
{
    PyObject *args;
    int ret;
    int _relax, _panel_size;

    static char *kwlist[] = {
        "Fact", "Equil", "ColPerm", "Trans", "IterRefine",
        "DiagPivotThresh", "PivotGrowth", "ConditionNumber",
        "RowPerm", "SymmetricMode", "PrintStat", "ReplaceTinyPivot",
        "SolveInitialized", "RefineInitialized", "ILU_Norm",
        "ILU_MILU", "ILU_DropTol", "ILU_FillTol", "ILU_FillFactor",
        "ILU_DropRule", "PanelSize", "Relax", NULL
    };

    if (ilu) {
        ilu_set_default_options(options);
    }
    else {
        set_default_options(options);
    }

    _panel_size = sp_ienv(1);
    _relax = sp_ienv(2);

    if (option_dict == NULL) {
        return 0;
    }
    
    args = PyTuple_New(0);
    ret = PyArg_ParseTupleAndKeywords(
        args, option_dict,
        "|O&O&O&O&O&O&O&O&O&O&O&O&O&O&O&O&O&O&O&O&O&O&", kwlist,
        fact_cvt, &options->Fact,
        yes_no_cvt, &options->Equil,
        colperm_cvt, &options->ColPerm,
        trans_cvt, &options->Trans,
        iterrefine_cvt, &options->IterRefine,
        double_cvt, &options->DiagPivotThresh,
        yes_no_cvt, &options->PivotGrowth,
        yes_no_cvt, &options->ConditionNumber,
        rowperm_cvt, &options->RowPerm,
        yes_no_cvt, &options->SymmetricMode,
        yes_no_cvt, &options->PrintStat,
        yes_no_cvt, &options->ReplaceTinyPivot,
        yes_no_cvt, &options->SolveInitialized,
        yes_no_cvt, &options->RefineInitialized,
        norm_cvt, &options->ILU_Norm,
        milu_cvt, &options->ILU_MILU,
        double_cvt, &options->ILU_DropTol,
        double_cvt, &options->ILU_FillTol,
        double_cvt, &options->ILU_FillFactor,
        droprule_cvt, &options->ILU_DropRule,
        int_cvt, &_panel_size,
        int_cvt, &_relax
        );
    Py_DECREF(args);

    if (panel_size != NULL) {
        *panel_size = _panel_size;
    }
    if (relax != NULL) {
        *relax = _relax;
    }

    return ret;
}
