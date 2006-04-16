#ifndef SPMATRIX_H
#define SPMATRIX_H

#include "ll_mat.h"
#include "csr_mat.h"
#include "sss_mat.h"

#include "spmatrix_api.h"

/*
 * Macro definitions
 */

/** SPMATRIX_PARSE_ARGS_ARR_ARR
 * 
 * Macro for parsing arguments for matvec and precon type operations
 */
#define SPMATRIX_PARSE_ARGS_ARR_ARR(args, arg1, arg2, n1, n2) \
  if (!PyArg_ParseTuple((args), "O!O!", &PyArray_Type, &(arg1), &PyArray_Type, &(arg2))) \
    return NULL; \
  \
  if ((arg1)->nd != 1 || \
      (arg1)->descr->type_num != PyArray_DOUBLE || \
      (arg1)->dimensions[0] != (n1) || \
      !((arg1)->flags & CONTIGUOUS)) { \
    PyErr_SetString(PyExc_ValueError, "arg 1 must be a contiguous 1-dimensional double array of appropriate size."); \
    return NULL; \
  } \
  \
  if ((arg2)->nd != 1 || \
      (arg2)->descr->type_num != PyArray_DOUBLE || \
      (arg2)->dimensions[0] != (n2) || \
      !((arg2)->flags & CONTIGUOUS)) { \
    PyErr_SetString(PyExc_ValueError, "arg 2 must be a contiguous 1-dimensional double array of appropriate size."); \
    return NULL; \
  }

#define SPMATRIX_PARSE_ARGS_ARR_ARR_STRIDE(args, arg1, arg2, n1, n2) \
  if (!PyArg_ParseTuple((args), "O!O!", &PyArray_Type, &(arg1), &PyArray_Type, &(arg2))) \
    return NULL; \
  \
  if ((arg1)->nd != 1 || \
      (arg1)->descr->type_num != PyArray_DOUBLE || \
      (arg1)->dimensions[0] != (n1)) { \
    PyErr_SetString(PyExc_ValueError, "arg 1 must be a 1-dimensional double array of appropriate size."); \
    return NULL; \
  } \
  \
  if ((arg2)->nd != 1 || \
      (arg2)->descr->type_num != PyArray_DOUBLE || \
      (arg2)->dimensions[0] != (n2)) { \
    PyErr_SetString(PyExc_ValueError, "arg 2 must be a 1-dimensional double array of appropriate size."); \
    return NULL; \
  }

/** SPMATRIX_CHECK_ARR_DIM_SIZE
 * 
 * Macro for checking the type and size of vector arguments
 */
#define SPMATRIX_CHECK_ARR_DIM_SIZE(arr, dim, size) \
  if ((arr)->nd != (dim) || \
      (arr)->descr->type_num != PyArray_DOUBLE || \
      (arr)->dimensions[0] != (size) || \
      !((arr)->flags & CONTIGUOUS)) { \
    PyErr_SetString(PyExc_ValueError, "argument must be a contiguous 1-dimensional double array of appropriate size."); \
    return NULL; \
  }

#ifdef SPMATRIX_MODULE
extern PyObject *SpMatrix_ErrorObject;
#endif

#endif
