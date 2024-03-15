/* -*-c-*-  */
/*
 * _superlu module
 *
 * Python interface to SuperLU decompositions.
 */

/* Copyright 1999 Travis Oliphant
 *
 * Permission to copy and modified this file is granted under
 * the revised BSD license. No warranty is expressed or IMPLIED
 */

#include <Python.h>

#define PY_ARRAY_UNIQUE_SYMBOL _scipy_sparse_superlu_ARRAY_API
#include <numpy/ndarrayobject.h>

#include "_superluobject.h"


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
    volatile PyObject *Py_B = NULL;
    volatile PyArrayObject *Py_X = NULL;
    volatile PyArrayObject *nzvals = NULL;
    volatile PyArrayObject *colind = NULL, *rowptr = NULL;
    volatile int N, nnz;
    volatile int info;
    volatile int csc = 0;
    volatile int *perm_r = NULL, *perm_c = NULL;
    volatile SuperMatrix A = { 0 }, B = { 0 }, L = { 0 }, U = { 0 };
    volatile superlu_options_t options = { 0 };
    volatile SuperLUStat_t stat = { 0 };
    volatile PyObject *option_dict = NULL;
    volatile int type;
    volatile jmp_buf *jmpbuf_ptr;
    SLU_BEGIN_THREADS_DEF;

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

    type = PyArray_TYPE((PyArrayObject*)nzvals);
    if (!CHECK_SLU_TYPE(type)) {
	PyErr_SetString(PyExc_TypeError,
			"nzvals is not of a type supported by SuperLU");
	return NULL;
    }

    if (!set_superlu_options_from_dict((superlu_options_t*)&options, 0,
                                       (PyObject*)option_dict, NULL, NULL)) {
	return NULL;
    }

    /* Create Space for output */
    Py_X = (PyArrayObject*)PyArray_FROMANY(
        (PyObject*)Py_B, type, 1, 2,
        NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_ENSURECOPY);
    if (Py_X == NULL)
	return NULL;

    if (PyArray_DIM((PyArrayObject*)Py_X, 0) != N) {
        PyErr_SetString(PyExc_ValueError,
                        "b array has invalid shape");
        Py_DECREF(Py_X);
        return NULL;
    }

    if (csc) {
	if (NCFormat_from_spMatrix((SuperMatrix*)&A, N, N, nnz,
                                   (PyArrayObject *)nzvals, (PyArrayObject *)colind,
                                   (PyArrayObject *)rowptr, type)) {
	    Py_DECREF(Py_X);
	    return NULL;
	}
    }
    else {
	if (NRFormat_from_spMatrix((SuperMatrix*)&A, N, N, nnz, (PyArrayObject *)nzvals,
                                   (PyArrayObject *)colind, (PyArrayObject *)rowptr,
				   type)) {
	    Py_DECREF(Py_X);
	    return NULL;
	}
    }

    if (DenseSuper_from_Numeric((SuperMatrix*)&B, (PyObject*)Py_X)) {
	Destroy_SuperMatrix_Store((SuperMatrix*)&A);
	Py_DECREF(Py_X);
	return NULL;
    }

    /* B and Py_X  share same data now but Py_X "owns" it */

    /* Setup options */

    jmpbuf_ptr = (volatile jmp_buf *)superlu_python_jmpbuf();
    SLU_BEGIN_THREADS;
    if (setjmp(*(jmp_buf*)jmpbuf_ptr)) {
        SLU_END_THREADS;
	goto fail;
    }
    else {
	perm_c = intMalloc(N);
	perm_r = intMalloc(N);
	StatInit((SuperLUStat_t*)&stat);

	/* Compute direct inverse of sparse Matrix */
	gssv(type, (superlu_options_t*)&options, (SuperMatrix*)&A, (int*)perm_c, (int*)perm_r,
             (SuperMatrix*)&L, (SuperMatrix*)&U, (SuperMatrix*)&B, (SuperLUStat_t*)&stat,
             (int*)&info);
        SLU_END_THREADS;
    }

    SUPERLU_FREE((void*)perm_r);
    SUPERLU_FREE((void*)perm_c);
    Destroy_SuperMatrix_Store((SuperMatrix*)&A);	/* holds just a pointer to the data */
    Destroy_SuperMatrix_Store((SuperMatrix*)&B);
    Destroy_SuperNode_Matrix((SuperMatrix*)&L);
    Destroy_CompCol_Matrix((SuperMatrix*)&U);
    StatFree((SuperLUStat_t*)&stat);

    return Py_BuildValue("Ni", Py_X, info);

  fail:
    SUPERLU_FREE((void*)perm_r);
    SUPERLU_FREE((void*)perm_c);
    XDestroy_SuperMatrix_Store((SuperMatrix*)&A);	/* holds just a pointer to the data */
    XDestroy_SuperMatrix_Store((SuperMatrix*)&B);
    XDestroy_SuperNode_Matrix((SuperMatrix*)&L);
    XDestroy_CompCol_Matrix((SuperMatrix*)&U);
    XStatFree((SuperLUStat_t*)&stat);
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
    PyObject *py_csc_construct_func = NULL;
    PyObject *option_dict = NULL;
    int type;
    int ilu = 0;

    static char *kwlist[] = { "N", "nnz", "nzvals", "colind", "rowptr",
        "csc_construct_func", "options", "ilu",
	NULL
    };

    int res =
	PyArg_ParseTupleAndKeywords(args, keywds, "iiO!O!O!O|Oi", kwlist,
				    &N, &nnz,
				    &PyArray_Type, &nzvals,
				    &PyArray_Type, &rowind,
				    &PyArray_Type, &colptr,
                                    &py_csc_construct_func,
				    &option_dict,
				    &ilu);

    if (!res)
	return NULL;

    if (!_CHECK_INTEGER(colptr) || !_CHECK_INTEGER(rowind)) {
	PyErr_SetString(PyExc_TypeError,
			"rowind and colptr must be of type cint");
	return NULL;
    }

    type = PyArray_TYPE((PyArrayObject*)nzvals);
    if (!CHECK_SLU_TYPE(type)) {
	PyErr_SetString(PyExc_TypeError,
			"nzvals is not of a type supported by SuperLU");
	return NULL;
    }

    if (NCFormat_from_spMatrix(&A, N, N, nnz, nzvals, rowind, colptr,
			       type)) {
	goto fail;
    }

    result = newSuperLUObject(&A, option_dict, type, ilu, py_csc_construct_func);
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
Parameters\n\
----------\n\
\n\
Matrix to be factorized is represented as N,nnz,nzvals,rowind,colptr\n\
  as separate arguments.  This is compressed sparse column representation.\n\
\n\
N : int \n\
    number of rows and columns\n\
nnz : int\n\
    number of non-zero elements\n\
nzvals : array\n\
    non-zero values \n\
rowind : array \n\
    row-index for this column (same size as nzvals)\n\
colptr : array \n\
    index into rowind for first non-zero value in this column\n\
    size is (N+1).  Last value should be nnz. \n\
\n\
Other Parameters\n\
----------------\n\
options\n\
    specifies additional options for SuperLU\n\
    (same keys and values as in superlu_options_t C structure,\n\
    and additionally 'Relax' and 'PanelSize')\n\
\n\
ilu : bool\n\
    whether to perform an incomplete LU decomposition\n\
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

PyMODINIT_FUNC
PyInit__superlu(void)
{
    PyObject *module, *mdict;

    import_array();

    if (PyType_Ready(&SuperLUType) < 0) {
        return NULL;
    }
    if (PyType_Ready(&SuperLUGlobalType) < 0) {
    	return NULL;
    }

    module = PyModule_Create(&moduledef);
    if (module == NULL) {
        return NULL;
    }
    mdict = PyModule_GetDict(module);
    if (mdict == NULL) {
        return NULL;
    }

    if (PyDict_SetItemString(mdict, "SuperLU", (PyObject *) &SuperLUType)) {
        return NULL;
    }

    return module;
}
