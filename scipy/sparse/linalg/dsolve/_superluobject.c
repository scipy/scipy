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
#include <ctype.h>


/*********************************************************************** 
 * SuperLUObject methods
 */

static PyObject *SuperLU_solve(SuperLUObject * self, PyObject * args,
			       PyObject * kwds)
{
    volatile PyArrayObject *b, *x = NULL;
    volatile SuperMatrix B = { 0 };
#ifndef NPY_PY3K
    volatile char itrans = 'N';
#else
    volatile int itrans = 'N';
#endif
    volatile int info;
    volatile trans_t trans;
    volatile SuperLUStat_t stat = { 0 };
    static char *kwlist[] = { "rhs", "trans", NULL };
    volatile jmp_buf *jmpbuf_ptr;
    SLU_BEGIN_THREADS_DEF;

    if (!CHECK_SLU_TYPE(self->type)) {
        PyErr_SetString(PyExc_ValueError, "unsupported data type");
        return NULL;
    }

#ifndef NPY_PY3K
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|c", kwlist,
                                     &PyArray_Type, &b, &itrans))
#else
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!|C", kwlist,
                                     &PyArray_Type, &b, &itrans))
#endif
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

    x = (PyArrayObject*)PyArray_FROMANY(
        (PyObject*)b, self->type, 1, 2,
        NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_ENSURECOPY);
    if (x == NULL) {
        goto fail;
    }

    if (PyArray_DIM((PyArrayObject*)x, 0) != self->n) {
        PyErr_SetString(PyExc_ValueError, "b is of incompatible size");
        goto fail;
    }

    if (DenseSuper_from_Numeric((SuperMatrix*)&B, (PyObject *)x))
        goto fail;

    jmpbuf_ptr = (volatile jmp_buf *)superlu_python_jmpbuf();
    if (setjmp(*(jmp_buf*)jmpbuf_ptr)) {
	goto fail;
    }

    StatInit((SuperLUStat_t *)&stat);

    /* Solve the system, overwriting vector x. */
    jmpbuf_ptr = (volatile jmp_buf *)superlu_python_jmpbuf();
    SLU_BEGIN_THREADS;
    if (setjmp(*(jmp_buf*)jmpbuf_ptr)) {
        SLU_END_THREADS;
	goto fail;
    }
    gstrs(self->type,
	  trans, &self->L, &self->U, self->perm_c, self->perm_r,
          (SuperMatrix *)&B, (SuperLUStat_t *)&stat, (int *)&info);
    SLU_END_THREADS;

    if (info) {
	PyErr_SetString(PyExc_SystemError,
			"gstrs was called with invalid arguments");
	goto fail;
    }

    /* free memory */
    Destroy_SuperMatrix_Store((SuperMatrix *)&B);
    StatFree((SuperLUStat_t *)&stat);
    return (PyObject *) x;

  fail:
    XDestroy_SuperMatrix_Store((SuperMatrix *)&B);
    XStatFree((SuperLUStat_t *)&stat);
    Py_XDECREF(x);
    return NULL;
}

/** table of object methods
 */
PyMethodDef SuperLU_methods[] = {
    {"solve", (PyCFunction) SuperLU_solve, METH_VARARGS | METH_KEYWORDS, NULL},
    {NULL, NULL}		/* sentinel */
};


/*********************************************************************** 
 * SuperLUType methods
 */

static void SuperLU_dealloc(SuperLUObject * self)
{
    Py_XDECREF(self->cached_U);
    Py_XDECREF(self->cached_L);
    self->cached_U = NULL;
    self->cached_L = NULL;
    SUPERLU_FREE(self->perm_r);
    SUPERLU_FREE(self->perm_c);
    self->perm_r = NULL;
    self->perm_c = NULL;
    XDestroy_SuperNode_Matrix(&self->L);
    XDestroy_CompCol_Matrix(&self->U);
    PyObject_Del(self);
}

static PyObject *SuperLU_getter(PyObject *selfp, void *data)
{
    SuperLUObject *self = (SuperLUObject *)selfp;
    char *name = (char*)data;

    if (strcmp(name, "shape") == 0) {
	return Py_BuildValue("(i,i)", self->m, self->n);
    }
    else if (strcmp(name, "nnz") == 0)
	return Py_BuildValue("i",
			     ((SCformat *) self->L.Store)->nnz +
			     ((SCformat *) self->U.Store)->nnz);
    else if (strcmp(name, "perm_r") == 0) {
	PyObject *perm_r;
        perm_r = PyArray_SimpleNewFromData(
            1, (npy_intp *) (&self->n), NPY_INT,
            (void *) self->perm_r);
        if (perm_r == NULL) {
            return NULL;
        }

	/* For ref counting of the memory */
	PyArray_SetBaseObject((PyArrayObject*)perm_r, (PyObject*)self);
	Py_INCREF(self);
	return perm_r;
    }
    else if (strcmp(name, "perm_c") == 0) {
	PyObject *perm_c;

        perm_c = PyArray_SimpleNewFromData(
            1, (npy_intp *) (&self->n), NPY_INT,
            (void *) self->perm_c);
        if (perm_c == NULL) {
            return NULL;
        }

	/* For ref counting of the memory */
	PyArray_SetBaseObject((PyArrayObject*)perm_c, (PyObject*)self);
	Py_INCREF(self);
	return perm_c;
    }
    else if (strcmp(name, "U") == 0 || strcmp(name, "L") == 0) {
        int ok;
        if (self->cached_U == NULL) {
            ok = LU_to_csc_matrix(&self->L, &self->U,
                                  &self->cached_L, &self->cached_U);
            if (ok != 0) {
                return NULL;
            }
        }
        if (strcmp(name, "U") == 0) {
            Py_INCREF(self->cached_U);
            return self->cached_U;
        }
        else {
            Py_INCREF(self->cached_L);
            return self->cached_L;
        }
    }
    else {
        PyErr_SetString(PyExc_RuntimeError,
                        "internal error (this is a bug)");
        return NULL;
    }
}


/***********************************************************************
 * SuperLUType structure
 */

PyGetSetDef SuperLU_getset[] = {
    {"shape", SuperLU_getter, (setter)NULL, (char*)NULL, (void*)"shape"},
    {"nnz", SuperLU_getter, (setter)NULL, (char*)NULL, (void*)"nnz"},
    {"perm_r", SuperLU_getter, (setter)NULL, (char*)NULL, (void*)"perm_r"},
    {"perm_c", SuperLU_getter, (setter)NULL, (char*)NULL, (void*)"perm_c"},
    {"U", SuperLU_getter, (setter)NULL, (char*)NULL, (void*)"U"},
    {"L", SuperLU_getter, (setter)NULL, (char*)NULL, (void*)"L"},
    NULL
};


PyTypeObject SuperLUType = {
#if defined(NPY_PY3K)
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(NULL)
    0,
#endif
    "SuperLU",
    sizeof(SuperLUObject),
    0,
    (destructor) SuperLU_dealloc, /* tp_dealloc */
    0,				/* tp_print */
    0,	                        /* tp_getattr */
    0,				/* tp_setattr */
    0,				/* tp_compare / tp_reserved */
    0,				/* tp_repr */
    0,				/* tp_as_number */
    0,				/* tp_as_sequence */
    0,				/* tp_as_mapping */
    0,				/* tp_hash */
    0,				/* tp_call */
    0,				/* tp_str */
    0,				/* tp_getattro */
    0,				/* tp_setattro */
    0,				/* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,		/* tp_flags */
    NULL,                       /* tp_doc */
    0,				/* tp_traverse */
    0,				/* tp_clear */
    0,				/* tp_richcompare */
    0,				/* tp_weaklistoffset */
    0,				/* tp_iter */
    0,				/* tp_iternext */
    SuperLU_methods,		/* tp_methods */
    0,				/* tp_members */
    SuperLU_getset,		/* tp_getset */
    0,				/* tp_base */
    0,				/* tp_dict */
    0,				/* tp_descr_get */
    0,				/* tp_descr_set */
    0,				/* tp_dictoffset */
    0,				/* tp_init */
    0,				/* tp_alloc */
    0,				/* tp_new */
    0,				/* tp_free */
    0,				/* tp_is_gc */
    0,				/* tp_bases */
    0,				/* tp_mro */
    0,				/* tp_cache */
    0,				/* tp_subclasses */
    0,				/* tp_weaklist */
    0,				/* tp_del */
    0,				/* tp_version_tag */
};


int DenseSuper_from_Numeric(SuperMatrix *X, PyObject *PyX)
{
    volatile PyArrayObject *aX;
    volatile int m, n, ldx, nd;
    volatile jmp_buf *jmpbuf_ptr;

    if (!PyArray_Check(PyX)) {
        PyErr_SetString(PyExc_TypeError,
                        "argument is not an array.");
        return -1;
    }

    aX = (PyArrayObject*)PyX;

    if (!CHECK_SLU_TYPE(PyArray_TYPE((PyArrayObject*)aX))) {
        PyErr_SetString(PyExc_ValueError, "unsupported array data type");
        return -1;
    }

    if (!(PyArray_FLAGS((PyArrayObject*)aX) & NPY_ARRAY_F_CONTIGUOUS)) {
        PyErr_SetString(PyExc_ValueError, "array is not fortran contiguous");
        return -1;
    }

    nd = PyArray_NDIM((PyArrayObject*)aX);

    if (nd == 1) {
	m = PyArray_DIM((PyArrayObject*)aX, 0);
	n = 1;
	ldx = m;
    }
    else if (nd == 2) {
	m = PyArray_DIM((PyArrayObject*)aX, 0);
	n = PyArray_DIM((PyArrayObject*)aX, 1);
	ldx = m;
    }
    else {
        PyErr_SetString(PyExc_ValueError, "wrong number of dimensions in array");
        return -1;
    }

    jmpbuf_ptr = (volatile jmp_buf *)superlu_python_jmpbuf();
    if (setjmp(*(jmp_buf*)jmpbuf_ptr)) {
	return -1;
    }
    else {
	Create_Dense_Matrix(PyArray_TYPE((PyArrayObject*)aX), X, m, n,
			    PyArray_DATA((PyArrayObject*)aX), ldx, SLU_DN,
			    NPY_TYPECODE_TO_SLU(PyArray_TYPE((PyArrayObject*)aX)),
			    SLU_GE);
    }
    return 0;
}

/* Natively handles Compressed Sparse Row and CSC */

int NRFormat_from_spMatrix(SuperMatrix * A, int m, int n, int nnz,
			   PyArrayObject * nzvals, PyArrayObject * colind,
			   PyArrayObject * rowptr, int typenum)
{
    volatile int ok = 0;
    volatile jmp_buf *jmpbuf_ptr;

    ok = (PyArray_EquivTypenums(PyArray_TYPE(nzvals), typenum) &&
          PyArray_EquivTypenums(PyArray_TYPE(colind), NPY_INT) &&
          PyArray_EquivTypenums(PyArray_TYPE(rowptr), NPY_INT) &&
          PyArray_NDIM(nzvals) == 1 &&
          PyArray_NDIM(colind) == 1 &&
          PyArray_NDIM(rowptr) == 1 &&
          PyArray_IS_C_CONTIGUOUS(nzvals) &&
          PyArray_IS_C_CONTIGUOUS(colind) &&
          PyArray_IS_C_CONTIGUOUS(rowptr) &&
          nnz <= PyArray_DIM(nzvals, 0) &&
          nnz <= PyArray_DIM(colind, 0) &&
          m+1 <= PyArray_DIM(rowptr, 0));
    if (!ok) {
	PyErr_SetString(PyExc_ValueError,
			"sparse matrix arrays must be 1-D C-contiguous and of proper "
                        "sizes and types");
	return -1;
    }

    jmpbuf_ptr = (volatile jmp_buf *)superlu_python_jmpbuf();
    if (setjmp(*(jmp_buf*)jmpbuf_ptr)) {
	return -1;
    }
    else {
	if (!CHECK_SLU_TYPE(PyArray_TYPE(nzvals))) {
	    PyErr_SetString(PyExc_TypeError, "Invalid type for array.");
	    return -1;
	}
	Create_CompRow_Matrix(PyArray_TYPE(nzvals),
			      A, m, n, nnz, PyArray_DATA((PyArrayObject*)nzvals),
			      (int *) PyArray_DATA((PyArrayObject*)colind),
                              (int *) PyArray_DATA((PyArrayObject*)rowptr),
			      SLU_NR,
			      NPY_TYPECODE_TO_SLU(PyArray_TYPE((PyArrayObject*)nzvals)),
			      SLU_GE);
    }

    return 0;
}

int NCFormat_from_spMatrix(SuperMatrix * A, int m, int n, int nnz,
			   PyArrayObject * nzvals, PyArrayObject * rowind,
			   PyArrayObject * colptr, int typenum)
{
    volatile int ok = 0;
    volatile jmp_buf *jmpbuf_ptr;

    ok = (PyArray_EquivTypenums(PyArray_TYPE(nzvals), typenum) &&
          PyArray_EquivTypenums(PyArray_TYPE(rowind), NPY_INT) &&
          PyArray_EquivTypenums(PyArray_TYPE(colptr), NPY_INT) &&
          PyArray_NDIM(nzvals) == 1 &&
          PyArray_NDIM(rowind) == 1 &&
          PyArray_NDIM(colptr) == 1 &&
          PyArray_IS_C_CONTIGUOUS(nzvals) &&
          PyArray_IS_C_CONTIGUOUS(rowind) &&
          PyArray_IS_C_CONTIGUOUS(colptr) &&
          nnz <= PyArray_DIM(nzvals, 0) &&
          nnz <= PyArray_DIM(rowind, 0) &&
          n+1 <= PyArray_DIM(colptr, 0));
    if (!ok) {
	PyErr_SetString(PyExc_ValueError,
			"sparse matrix arrays must be 1-D C-contiguous and of proper "
                        "sizes and types");
	return -1;
    }

    jmpbuf_ptr = (volatile jmp_buf *)superlu_python_jmpbuf();
    if (setjmp(*(jmp_buf*)jmpbuf_ptr)) {
	return -1;
    }
    else {
	if (!CHECK_SLU_TYPE(PyArray_TYPE(nzvals))) {
	    PyErr_SetString(PyExc_TypeError, "Invalid type for array.");
	    return -1;
	}
	Create_CompCol_Matrix(PyArray_TYPE(nzvals),
			      A, m, n, nnz, PyArray_DATA(nzvals),
			      (int *) PyArray_DATA(rowind), (int *) PyArray_DATA(colptr),
			      SLU_NC,
			      NPY_TYPECODE_TO_SLU(PyArray_TYPE(nzvals)),
			      SLU_GE);
    }

    return 0;
}


/*
 * Create Scipy sparse matrices out from Superlu LU decomposition.
 */

static int LU_to_csc(SuperMatrix *L, SuperMatrix *U,
                     int *U_indices, int *U_indptr, char *U_data,
                     int *L_indices, int *L_indptr, char *L_data,
                     Dtype_t dtype);

int LU_to_csc_matrix(SuperMatrix *L, SuperMatrix *U,
                     PyObject **L_csc, PyObject **U_csc)
{
    SCformat *Lstore;
    NCformat *Ustore;
    PyObject *U_indices = NULL, *U_indptr = NULL, *U_data = NULL;
    PyObject *L_indices = NULL, *L_indptr = NULL, *L_data = NULL;
    PyObject *scipy_sparse = NULL, *datatuple = NULL, *shape = NULL;
    int result = -1, ok;
    int type;
    npy_intp dims[1];

    *L_csc = NULL;
    *U_csc = NULL;

    if (U->Stype != SLU_NC || L->Stype != SLU_SC ||
        U->Mtype != SLU_TRU || L->Mtype != SLU_TRLU ||
        L->nrow != U->nrow || L->ncol != L->nrow ||
        U->ncol != U->nrow || L->Dtype != U->Dtype)
    {
        PyErr_SetString(PyExc_RuntimeError,
                        "internal error: invalid Superlu matrix data");
        return -1;
    }

    Ustore = (NCformat*)U->Store;
    Lstore = (SCformat*)L->Store;

    type = SLU_TYPECODE_TO_NPY(L->Dtype);

    /* Allocate output */
#define CREATE_1D_ARRAY(name, type, size)               \
        do {                                            \
            dims[0] = size;                             \
            name = PyArray_EMPTY(1, dims, type, 0);     \
            if (name == NULL) goto fail;                \
        } while (0)

    CREATE_1D_ARRAY(L_indices, NPY_INT, Lstore->nnz);
    CREATE_1D_ARRAY(L_indptr, NPY_INT, L->ncol + 1);
    CREATE_1D_ARRAY(L_data, type, Lstore->nnz);

    CREATE_1D_ARRAY(U_indices, NPY_INT, Ustore->nnz);
    CREATE_1D_ARRAY(U_indptr, NPY_INT, U->ncol + 1);
    CREATE_1D_ARRAY(U_data, type, Ustore->nnz);

#undef CREATE_1D_ARRAY

    /* Copy data over */
    ok = LU_to_csc(
        L, U,
        (int*)PyArray_DATA((PyArrayObject*)L_indices),
        (int*)PyArray_DATA((PyArrayObject*)L_indptr),
        (void*)PyArray_DATA((PyArrayObject*)L_data),
        (int*)PyArray_DATA((PyArrayObject*)U_indices),
        (int*)PyArray_DATA((PyArrayObject*)U_indptr),
        (void*)PyArray_DATA((PyArrayObject*)U_data),
        L->Dtype
        );

    if (ok != 0) {
        goto fail;
    }

    /* Create sparse matrices */
    scipy_sparse = PyImport_ImportModule("scipy.sparse");
    if (scipy_sparse == NULL) {
        goto fail;
    }

    shape = Py_BuildValue("ii", L->nrow, L->ncol);
    if (shape == NULL) {
        goto fail;
    }

    datatuple = Py_BuildValue("OOO", L_data, L_indices, L_indptr);
    if (datatuple == NULL) {
        goto fail;
    }
    *L_csc = PyObject_CallMethod(scipy_sparse, "csc_matrix",
                                 "OO", datatuple, shape);
    if (*L_csc == NULL) {
        goto fail;
    }

    Py_DECREF(datatuple);
    datatuple = Py_BuildValue("OOO", U_data, U_indices, U_indptr);
    if (datatuple == NULL) {
        Py_DECREF(*L_csc);
        *L_csc = NULL;
        goto fail;
    }
    *U_csc = PyObject_CallMethod(scipy_sparse, "csc_matrix",
                                 "OO", datatuple, shape);
    if (*U_csc == NULL) {
        Py_DECREF(*L_csc);
        *L_csc = NULL;
        goto fail;
    }

    result = 0;
    
fail:
    Py_XDECREF(U_indices);
    Py_XDECREF(U_indptr);
    Py_XDECREF(U_data);
    Py_XDECREF(L_indices);
    Py_XDECREF(L_indptr);
    Py_XDECREF(L_data);
    Py_XDECREF(shape);
    Py_XDECREF(scipy_sparse);
    Py_XDECREF(datatuple);

    return result;
}


/*
 * Convert SuperLU L and U matrices to CSC format.
 *
 * The LU decomposition U factor is partly stored in U and partly in the upper
 * diagonal of L.  The L matrix is stored in column-addressable rectangular
 * superblock format.
 *
 * This routine is partly adapted from SuperLU MATLAB wrappers and the
 * SuperLU Print_SuperNode_Matrix routine.
 */
static int
LU_to_csc(SuperMatrix *L, SuperMatrix *U,
          int *L_rowind, int *L_colptr, char *L_data,
          int *U_rowind, int *U_colptr, char *U_data,
          Dtype_t dtype)
{
    SCformat *Lstore;
    NCformat *Ustore;
    npy_intp elsize;
    int isup, icol, icolstart, icolend, iptr, istart, iend;
    char *src, *dst;
    int U_nnz, L_nnz;

    Ustore = (NCformat*)U->Store;
    Lstore = (SCformat*)L->Store;

    switch (dtype) {
    case SLU_S: elsize = 4; break;
    case SLU_D: elsize = 8; break;
    case SLU_C: elsize = 8; break;
    case SLU_Z: elsize = 16; break;
    default:
        /* shouldn't occur */
        PyErr_SetString(PyExc_ValueError, "unknown dtype");
        return -1;
    }
    
#define IS_ZERO(p)                                                      \
    ((dtype == SLU_S) ? (*(float*)(p) == 0) :                           \
     ((dtype == SLU_D) ? (*(double*)(p) == 0) :                         \
      ((dtype == SLU_C) ? (*(float*)(p) == 0 && *((float*)(p)+1) == 0) : \
       (*(double*)(p) == 0 && *((double*)(p)+1) == 0))))

    U_colptr[0] = 0;
    L_colptr[0] = 0;
    U_nnz = 0;
    L_nnz = 0;

    /* For each supernode */
    for (isup = 0; isup <= Lstore->nsuper; ++isup) {
        icolstart = Lstore->sup_to_col[isup];
        icolend = Lstore->sup_to_col[isup+1];
        istart = Lstore->rowind_colptr[icolstart];
        iend = Lstore->rowind_colptr[icolstart+1];

        /* For each column in supernode */
        for (icol = icolstart; icol < icolend; ++icol) {

            /* Process data in Ustore */
            for (iptr = Ustore->colptr[icol]; iptr < Ustore->colptr[icol+1]; ++iptr) {
                src = (char*)Ustore->nzval + elsize * iptr;
                if (!IS_ZERO(src)) {
                    if (U_nnz >= Ustore->nnz)
                        goto size_error;
                    U_rowind[U_nnz] = Ustore->rowind[iptr];
                    /* "U_data[U_nnz] = Ustore->nzvals[iptr]" */
                    dst = U_data + elsize * U_nnz;
                    memcpy(dst, src, elsize);
                    ++U_nnz;
                }
            }

            /* Process data in Lstore */
            src = (char*)Lstore->nzval + elsize * Lstore->nzval_colptr[icol];
            iptr = istart;

            /* Upper triangle part */
            for (; iptr < iend; ++iptr) {
                if (Lstore->rowind[iptr] > icol) {
                    break;
                }
                if (!IS_ZERO(src)) {
                    if (U_nnz >= Ustore->nnz)
                        goto size_error;
                    U_rowind[U_nnz] = Lstore->rowind[iptr];
                    dst = U_data + elsize * U_nnz;
                    memcpy(dst, src, elsize);
                    ++U_nnz;
                }
                src += elsize;
            }

            /* Add unit diagonal in L */
            if (L_nnz >= Lstore->nnz) return -1;
            dst = L_data + elsize * L_nnz;
            switch (dtype) {
            case SLU_S: *(float*)dst = 1.0; break;
            case SLU_D: *(double*)dst = 1.0; break;
            case SLU_C: *(float*)dst = 1.0; *((float*)dst+1) = 0.0; break;
            case SLU_Z: *(double*)dst = 1.0; *((double*)dst+1) = 0.0; break;
            }
            L_rowind[L_nnz] = icol;
            ++L_nnz;

            /* Lower triangle part */
            for (; iptr < iend; ++iptr) {
                if (!IS_ZERO(src)) {
                    if (L_nnz >= Lstore->nnz)
                         goto size_error;
                    L_rowind[L_nnz] = Lstore->rowind[iptr];
                    dst = L_data + elsize * L_nnz;
                    memcpy(dst, src, elsize);
                    ++L_nnz;
                }
                src += elsize;
            }

            /* Record column pointers */
            U_colptr[icol+1] = U_nnz;
            L_colptr[icol+1] = L_nnz;
        }
    }

    return 0;

size_error:
    PyErr_SetString(PyExc_RuntimeError,
                    "internal error: superlu matrixes have wrong nnz");
    return -1;
}


PyObject *newSuperLUObject(SuperMatrix * A, PyObject * option_dict,
                           int intype, int ilu)
{

    /* A must be in SLU_NC format used by the factorization routine. */
    volatile SuperLUObject *self;
    volatile SuperMatrix AC = { 0 };	/* Matrix postmultiplied by Pc */
    volatile int lwork = 0;
    volatile int *etree = NULL;
    volatile int info;
    volatile int n;
    volatile superlu_options_t options;
    volatile SuperLUStat_t stat = { 0 };
    volatile int panel_size, relax;
    volatile GlobalLU_t Glu;
    static volatile GlobalLU_t static_Glu;
    volatile GlobalLU_t *Glu_ptr;
    volatile jmp_buf *jmpbuf_ptr;
    SLU_BEGIN_THREADS_DEF;

    n = A->ncol;

    if (!set_superlu_options_from_dict((superlu_options_t*)&options, ilu, option_dict,
				       (int*)&panel_size, (int*)&relax)) {
	return NULL;
    }

    /* Create SLUObject */
    self = PyObject_New(SuperLUObject, &SuperLUType);
    if (self == NULL)
	return PyErr_NoMemory();
    self->m = A->nrow;
    self->n = n;
    self->perm_r = NULL;
    self->perm_c = NULL;
    self->L.Store = NULL;
    self->U.Store = NULL;
    self->cached_U = NULL;
    self->cached_L = NULL;
    self->type = intype;

    jmpbuf_ptr = (volatile jmp_buf *)superlu_python_jmpbuf();
    if (setjmp(*(jmp_buf*)jmpbuf_ptr)) {
	goto fail;
    }

    /* Calculate and apply minimum degree ordering */
    etree = intMalloc(n);
    self->perm_r = intMalloc(n);
    self->perm_c = intMalloc(n);
    StatInit((SuperLUStat_t *)&stat);

    /* calc column permutation */
    get_perm_c(options.ColPerm, A, self->perm_c);	

    /* apply column permutation */
    sp_preorder((superlu_options_t*)&options, A, self->perm_c, (int*)etree,
                (SuperMatrix*)&AC);

    /* Perform factorization */
    if (!CHECK_SLU_TYPE(SLU_TYPECODE_TO_NPY(A->Dtype))) {
	PyErr_SetString(PyExc_ValueError, "Invalid type in SuperMatrix.");
	goto fail;
    }

    if (options.Fact == SamePattern || options.Fact == SamePattern_SameRowPerm) {
        /* XXX: not nice, a better new API should be introduced for this */
        Glu_ptr = &static_Glu;
    }
    else {
        Glu_ptr = &Glu;
        jmpbuf_ptr = (volatile jmp_buf *)superlu_python_jmpbuf();
        SLU_BEGIN_THREADS;
        if (setjmp(*(jmp_buf*)jmpbuf_ptr)) {
            SLU_END_THREADS;
            goto fail;
        }
    }

    if (ilu) {
        gsitrf(SLU_TYPECODE_TO_NPY(A->Dtype),
               (superlu_options_t*)&options, (SuperMatrix*)&AC, relax, panel_size,
               (int*)etree, NULL, lwork, self->perm_c, self->perm_r,
               (SuperMatrix*)&self->L, (SuperMatrix*)&self->U, (GlobalLU_t*)Glu_ptr,
               (SuperLUStat_t*)&stat, (int*)&info);
    }
    else {
	gstrf(SLU_TYPECODE_TO_NPY(A->Dtype),
	      (superlu_options_t*)&options, (SuperMatrix*)&AC, relax, panel_size,
	      (int*)etree, NULL, lwork, self->perm_c, self->perm_r,
	      (SuperMatrix*)&self->L, (SuperMatrix*)&self->U, (GlobalLU_t*)Glu_ptr,
              (SuperLUStat_t*)&stat, (int*)&info);
    }

    SLU_END_THREADS;

    if (info) {
	if (info < 0)
	    PyErr_SetString(PyExc_SystemError,
			    "gstrf was called with invalid arguments");
	else {
	    if (info <= n)
		PyErr_SetString(PyExc_RuntimeError,
				"Factor is exactly singular");
	    else
		PyErr_NoMemory();
	}
	goto fail;
    }

    /* free memory */
    SUPERLU_FREE((void*)etree);
    Destroy_CompCol_Permuted((SuperMatrix*)&AC);
    StatFree((SuperLUStat_t*)&stat);

    return (PyObject *) self;

  fail:
    SUPERLU_FREE((void*)etree);
    XDestroy_CompCol_Permuted((SuperMatrix*)&AC);
    XStatFree((SuperLUStat_t*)&stat);
    Py_DECREF(self);
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

#define ENUM_CHECK_NAME(name, sname)                    \
    if (my_strxcmp(s, sname) == 0 || i == (long)name) { \
        *value = name;                                  \
        Py_XDECREF(tmpobj);                             \
        return 1;                                       \
    }

#define ENUM_CHECK(name) ENUM_CHECK_NAME(name, #name)

/*
 * Compare strings ignoring case, underscores and whitespace
 */
static int my_strxcmp(const char *a, const char *b)
{
    int c;

    while (*a != '\0' && *b != '\0') {
	while (*a == '_' || isspace(*a))
	    ++a;
	while (*b == '_' || isspace(*b))
	    ++b;
	c = (int) tolower(*a) - (int) tolower(*b);
	if (c != 0) {
	    return c;
	}
	++a;
	++b;
    }
    return (int) tolower(*a) - (int) tolower(*b);
}

static int yes_no_cvt(PyObject * input, yes_no_t * value)
{
    if (input == Py_None) {
	return 1;
    }
    else if (input == Py_True) {
	*value = YES;
    }
    else if (input == Py_False) {
	*value = NO;
    }
    else {
	PyErr_SetString(PyExc_ValueError, "value not a boolean");
	return 0;
    }
    return 1;
}

static int fact_cvt(PyObject * input, fact_t * value)
{
    ENUM_CHECK_INIT;
    ENUM_CHECK(DOFACT);
    ENUM_CHECK(SamePattern);
    ENUM_CHECK(SamePattern_SameRowPerm);
    ENUM_CHECK(FACTORED);
    ENUM_CHECK_FINISH("invalid value for 'Fact' parameter");
}

static int rowperm_cvt(PyObject * input, rowperm_t * value)
{
    ENUM_CHECK_INIT;
    ENUM_CHECK(NOROWPERM);
    ENUM_CHECK(LargeDiag);
    ENUM_CHECK(MY_PERMR);
    ENUM_CHECK_FINISH("invalid value for 'RowPerm' parameter");
}

static int colperm_cvt(PyObject * input, colperm_t * value)
{
    ENUM_CHECK_INIT;
    ENUM_CHECK(NATURAL);
    ENUM_CHECK(MMD_ATA);
    ENUM_CHECK(MMD_AT_PLUS_A);
    ENUM_CHECK(COLAMD);
    ENUM_CHECK(MY_PERMC);
    ENUM_CHECK_FINISH("invalid value for 'ColPerm' parameter");
}

static int trans_cvt(PyObject * input, trans_t * value)
{
    ENUM_CHECK_INIT;
    ENUM_CHECK(NOTRANS);
    ENUM_CHECK(TRANS);
    ENUM_CHECK(CONJ);
    if (my_strxcmp(s, "N") == 0) {
	*value = NOTRANS;
	return 1;
    }
    if (my_strxcmp(s, "T") == 0) {
	*value = TRANS;
	return 1;
    }
    if (my_strxcmp(s, "H") == 0) {
	*value = CONJ;
	return 1;
    }
    ENUM_CHECK_FINISH("invalid value for 'Trans' parameter");
}

static int iterrefine_cvt(PyObject * input, IterRefine_t * value)
{
    ENUM_CHECK_INIT;
    ENUM_CHECK(NOREFINE);
    ENUM_CHECK(SLU_SINGLE);
    ENUM_CHECK_NAME(SLU_SINGLE, "SINGLE");
    ENUM_CHECK(SLU_DOUBLE);
    ENUM_CHECK_NAME(SLU_DOUBLE, "DOUBLE");
    ENUM_CHECK(SLU_EXTRA);
    ENUM_CHECK_NAME(SLU_EXTRA, "EXTRA");
    ENUM_CHECK_FINISH("invalid value for 'IterRefine' parameter");
}

static int norm_cvt(PyObject * input, norm_t * value)
{
    ENUM_CHECK_INIT;
    ENUM_CHECK(ONE_NORM);
    ENUM_CHECK(TWO_NORM);
    ENUM_CHECK(INF_NORM);
    ENUM_CHECK_FINISH("invalid value for 'ILU_Norm' parameter");
}

static int milu_cvt(PyObject * input, milu_t * value)
{
    ENUM_CHECK_INIT;
    ENUM_CHECK(SILU);
    ENUM_CHECK(SMILU_1);
    ENUM_CHECK(SMILU_2);
    ENUM_CHECK(SMILU_3);
    ENUM_CHECK_FINISH("invalid value for 'ILU_MILU' parameter");
}

static int droprule_one_cvt(PyObject * input, int *value)
{
    ENUM_CHECK_INIT;
    if (my_strxcmp(s, "BASIC") == 0) {
	*value = DROP_BASIC;
	return 1;
    }
    if (my_strxcmp(s, "PROWS") == 0) {
	*value = DROP_PROWS;
	return 1;
    }
    if (my_strxcmp(s, "COLUMN") == 0) {
	*value = DROP_COLUMN;
	return 1;
    }
    if (my_strxcmp(s, "AREA") == 0) {
	*value = DROP_AREA;
	return 1;
    }
    if (my_strxcmp(s, "SECONDARY") == 0) {
	*value = DROP_SECONDARY;
	return 1;
    }
    if (my_strxcmp(s, "DYNAMIC") == 0) {
	*value = DROP_DYNAMIC;
	return 1;
    }
    if (my_strxcmp(s, "INTERP") == 0) {
	*value = DROP_INTERP;
	return 1;
    }
    ENUM_CHECK_FINISH("invalid value for 'ILU_DropRule' parameter");
}

static int droprule_cvt(PyObject * input, int *value)
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
    else if (PyString_Check(input) || PyUnicode_Check(input)) {
        /* Comma-separated string */
        char *fmt = "s";
#if PY_MAJOR_VERSION >= 3
        if (PyBytes_Check(input)) {
            fmt = "y";
        }
#endif
	seq = PyObject_CallMethod(input, "split", fmt, ",");
	if (seq == NULL || !PySequence_Check(seq))
	    goto fail;
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
	int one_value = 0;

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

static int double_cvt(PyObject * input, double *value)
{
    if (input == Py_None)
	return 1;
    *value = PyFloat_AsDouble(input);
    if (PyErr_Occurred())
	return 0;
    return 1;
}

static int int_cvt(PyObject * input, int *value)
{
    if (input == Py_None)
	return 1;
    *value = PyInt_AsLong(input);
    if (PyErr_Occurred())
	return 0;
    return 1;
}

int set_superlu_options_from_dict(superlu_options_t * options,
				  int ilu, PyObject * option_dict,
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
        /* Proceed with default options */
        ret = 1;
    }
    else {
        args = PyTuple_New(0);
        ret = PyArg_ParseTupleAndKeywords(args, option_dict,
                                          "|O&O&O&O&O&O&O&O&O&O&O&O&O&O&O&O&O&O&O&O&O&O&",
                                          kwlist, fact_cvt, &options->Fact,
                                          yes_no_cvt, &options->Equil,
                                          colperm_cvt, &options->ColPerm,
                                          trans_cvt, &options->Trans,
                                          iterrefine_cvt, &options->IterRefine,
                                          double_cvt,
                                          &options->DiagPivotThresh,
                                          yes_no_cvt, &options->PivotGrowth,
                                          yes_no_cvt,
                                          &options->ConditionNumber,
                                          rowperm_cvt, &options->RowPerm,
                                          yes_no_cvt, &options->SymmetricMode,
                                          yes_no_cvt, &options->PrintStat,
                                          yes_no_cvt,
                                          &options->ReplaceTinyPivot,
                                          yes_no_cvt,
                                          &options->SolveInitialized,
                                          yes_no_cvt,
                                          &options->RefineInitialized,
                                          norm_cvt, &options->ILU_Norm,
                                          milu_cvt, &options->ILU_MILU,
                                          double_cvt, &options->ILU_DropTol,
                                          double_cvt, &options->ILU_FillTol,
                                          double_cvt, &options->ILU_FillFactor,
                                          droprule_cvt, &options->ILU_DropRule,
                                          int_cvt, &_panel_size, int_cvt,
                                          &_relax);
        Py_DECREF(args);
    }

    if (panel_size != NULL) {
	*panel_size = _panel_size;
    }

    if (relax != NULL) {
	*relax = _relax;
    }

    return ret;
}
