/* -*-c-*-  */
/*
 * _superlu module
 *
 * Python interface to SuperLU decompositions.
 */

/* Copyright 1999 Travis Oliphant
 *
 * Permision to copy and modified this file is granted under
 * the revised BSD license. No warranty is expressed or IMPLIED
 */

#include <Python.h>
#include <setjmp.h>

#define PY_ARRAY_UNIQUE_SYMBOL _scipy_sparse_superlu_ARRAY_API
#include <numpy/arrayobject.h>

#include "_superluobject.h"
#include "numpy/npy_3kcompat.h"

extern jmp_buf _superlu_py_jmpbuf;

/*
 * NULL-safe deconstruction functions
 */
void XDestroy_SuperMatrix_Store(SuperMatrix * A)
{
    Destroy_SuperMatrix_Store(A);	/* safe as-is */
    A->Store = NULL;
}

void XDestroy_SuperNode_Matrix(SuperMatrix * A)
{
    if (A->Store) {
	Destroy_SuperNode_Matrix(A);
    }
    A->Store = NULL;
}

void XDestroy_CompCol_Matrix(SuperMatrix * A)
{
    if (A->Store) {
	Destroy_CompCol_Matrix(A);
    }
    A->Store = NULL;
}

void XDestroy_CompCol_Permuted(SuperMatrix * A)
{
    if (A->Store) {
	Destroy_CompCol_Permuted(A);
    }
    A->Store = NULL;
}

void XStatFree(SuperLUStat_t * stat)
{
    if (stat->ops) {
	StatFree(stat);
    }
    stat->ops = NULL;
}


/*
 * Data-type dependent implementations for Xgssv and Xgstrf;
 *
 * These have to included from separate files because of SuperLU include
 * structure.
 */

static PyObject *Py_gssv(PyObject * self, PyObject * args,
			 PyObject * kwdict)
{
    PyObject *Py_B = NULL;
    PyArrayObject *Py_X = NULL;
    PyArrayObject *nzvals = NULL;
    PyArrayObject *colind = NULL, *rowptr = NULL;
    int N, nnz;
    int info;
    int csc = 0;
    int *perm_r = NULL, *perm_c = NULL;
    SuperMatrix A = { 0 }, B = {
    0}, L = {
    0}, U = {
    0};
    superlu_options_t options = { 0 };
    SuperLUStat_t stat = { 0 };
    PyObject *option_dict = NULL;
    int type;
    int ssv_finished = 0;

    static char *kwlist[] = {
        "N", "nnz", "nzvals", "colind", "rowptr", "B", "csc",
        "options", NULL
    };

    /* Get input arguments */
    if (!PyArg_ParseTupleAndKeywords(args, kwdict, "iiO!O!O!O|iO", kwlist,
				     &N, &nnz, &PyArray_Type, &nzvals,
				     &PyArray_Type, &colind, &PyArray_Type,
				     &rowptr, &Py_B, &csc, &option_dict)) {
	return NULL;
    }

    if (!_CHECK_INTEGER(colind) || !_CHECK_INTEGER(rowptr)) {
	PyErr_SetString(PyExc_TypeError,
			"colind and rowptr must be of type cint");
	return NULL;
    }

    type = PyArray_TYPE(nzvals);
    if (!CHECK_SLU_TYPE(type)) {
	PyErr_SetString(PyExc_TypeError,
			"nzvals is not of a type supported by SuperLU");
	return NULL;
    }

    if (!set_superlu_options_from_dict
	(&options, 0, option_dict, NULL, NULL)) {
	return NULL;
    }

    /* Create Space for output */
    Py_X = (PyArrayObject*)PyArray_FROMANY(
        Py_B, type, 1, 2,
        NPY_F_CONTIGUOUS | NPY_ENSURECOPY);
    if (Py_X == NULL)
	return NULL;

    if (Py_X->dimensions[0] != N) {
        PyErr_SetString(PyExc_ValueError,
                        "b array has invalid shape");
        Py_DECREF(Py_X);
        return NULL;
    }

    if (csc) {
	if (NCFormat_from_spMatrix(&A, N, N, nnz, nzvals, colind, rowptr,
				   type)) {
	    Py_DECREF(Py_X);
	    return NULL;
	}
    }
    else {
	if (NRFormat_from_spMatrix(&A, N, N, nnz, nzvals, colind, rowptr,
				   type)) {
	    Py_DECREF(Py_X);
	    return NULL;
	}
    }

    if (DenseSuper_from_Numeric(&B, Py_X)) {
	Destroy_SuperMatrix_Store(&A);
	Py_DECREF(Py_X);
	return NULL;
    }

    /* B and Py_X  share same data now but Py_X "owns" it */

    /* Setup options */

    if (setjmp(_superlu_py_jmpbuf)) {
	goto fail;
    }
    else {
	perm_c = intMalloc(N);
	perm_r = intMalloc(N);
	StatInit(&stat);

	/* Compute direct inverse of sparse Matrix */
	gssv(type, &options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
    }
    ssv_finished = 1;

    SUPERLU_FREE(perm_r);
    SUPERLU_FREE(perm_c);
    Destroy_SuperMatrix_Store(&A);	/* holds just a pointer to the data */
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    StatFree(&stat);

    return Py_BuildValue("Ni", Py_X, info);

  fail:
    SUPERLU_FREE(perm_r);
    SUPERLU_FREE(perm_c);
    XDestroy_SuperMatrix_Store(&A);	/* holds just a pointer to the data */
    XDestroy_SuperMatrix_Store(&B);
    XDestroy_SuperNode_Matrix(&L);
    XDestroy_CompCol_Matrix(&U);
    XStatFree(&stat);
    Py_XDECREF(Py_X);
    return NULL;
}

static PyObject *Py_gstrf(PyObject * self, PyObject * args,
			  PyObject * keywds)
{
    /* default value for SuperLU parameters */
    int N, nnz;
    PyArrayObject *rowind, *colptr, *nzvals;
    SuperMatrix A = { 0 };
    PyObject *result;
    PyObject *option_dict = NULL;
    int type;
    int ilu = 0;

    static char *kwlist[] = { "N", "nnz", "nzvals", "colind", "rowptr",
	"options", "ilu",
	NULL
    };

    int res =
	PyArg_ParseTupleAndKeywords(args, keywds, "iiO!O!O!|Oi", kwlist,
				    &N, &nnz,
				    &PyArray_Type, &nzvals,
				    &PyArray_Type, &rowind,
				    &PyArray_Type, &colptr,
				    &option_dict,
				    &ilu);

    if (!res)
	return NULL;

    if (!_CHECK_INTEGER(colptr) || !_CHECK_INTEGER(rowind)) {
	PyErr_SetString(PyExc_TypeError,
			"rowind and colptr must be of type cint");
	return NULL;
    }

    type = PyArray_TYPE(nzvals);
    if (!CHECK_SLU_TYPE(type)) {
	PyErr_SetString(PyExc_TypeError,
			"nzvals is not of a type supported by SuperLU");
	return NULL;
    }

    if (NCFormat_from_spMatrix(&A, N, N, nnz, nzvals, rowind, colptr,
			       type)) {
	goto fail;
    }

    result = newSuperLUObject(&A, option_dict, type, ilu);
    if (result == NULL) {
	goto fail;
    }

    /* arrays of input matrix will not be freed */
    Destroy_SuperMatrix_Store(&A);
    return result;

  fail:
    /* arrays of input matrix will not be freed */
    XDestroy_SuperMatrix_Store(&A);
    return NULL;
}

static char gssv_doc[] =
    "Direct inversion of sparse matrix.\n\nX = gssv(A,B) solves A*X = B for X.";

static char gstrf_doc[] = "gstrf(A, ...)\n\
\n\
performs a factorization of the sparse matrix A=*(N,nnz,nzvals,rowind,colptr) and \n\
returns a factored_lu object.\n\
\n\
arguments\n\
---------\n\
\n\
Matrix to be factorized is represented as N,nnz,nzvals,rowind,colptr\n\
  as separate arguments.  This is compressed sparse column representation.\n\
\n\
N         number of rows and columns \n\
nnz       number of non-zero elements\n\
nzvals    non-zero values \n\
rowind    row-index for this column (same size as nzvals)\n\
colptr    index into rowind for first non-zero value in this column\n\
          size is (N+1).  Last value should be nnz. \n\
\n\
additional keyword arguments:\n\
-----------------------------\n\
options             specifies additional options for SuperLU\n\
                    (same keys and values as in superlu_options_t C structure,\n\
                    and additionally 'Relax' and 'PanelSize')\n\
\n\
ilu                 whether to perform an incomplete LU decomposition\n\
                    (default: false)\n\
";


/*
 * Main SuperLU module
 */

static PyMethodDef SuperLU_Methods[] = {
    {"gssv", (PyCFunction) Py_gssv, METH_VARARGS | METH_KEYWORDS,
     gssv_doc},
    {"gstrf", (PyCFunction) Py_gstrf, METH_VARARGS | METH_KEYWORDS,
     gstrf_doc},
    {NULL, NULL}
};

#if PY_VERSION_HEX >= 0x03000000

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_superlu",
    NULL,
    -1,
    SuperLU_Methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyObject *PyInit__superlu(void)
{
    PyObject *m, *d;

    import_array();

    if (PyType_Ready(&SuperLUType) < 0) {
	return NULL;
    }

    m = PyModule_Create(&moduledef);
    d = PyModule_GetDict(m);

    Py_INCREF(&PyArrayFlags_Type);
    PyDict_SetItemString(d, "SuperLU",
			 (PyObject *) &SuperLUType);

    if (PyErr_Occurred())
	Py_FatalError("can't initialize module _superlu");

    return m;
}

#else

PyMODINIT_FUNC init_superlu(void)
{
    PyObject *m, *d;

    import_array();

    SuperLUType.ob_type = &PyType_Type;
    if (PyType_Ready(&SuperLUType) < 0) {
	return;
    }

    m = Py_InitModule("_superlu", SuperLU_Methods);
    d = PyModule_GetDict(m);

    Py_INCREF(&PyArrayFlags_Type);
    PyDict_SetItemString(d, "SuperLU",
			 (PyObject *) & SuperLUType);
}

#endif
