/**
 * Author: Damian Eads
 * Date:   September 22, 2007 (moved to new file on June 8, 2008)
 * Adapted for incorporation into Scipy, April 9, 2008.
 *
 * Copyright (c) 2007, Damian Eads. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   - Redistributions of source code must retain the above
 *     copyright notice, this list of conditions and the
 *     following disclaimer.
 *   - Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer
 *     in the documentation and/or other materials provided with the
 *     distribution.
 *   - Neither the name of the author nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#if !defined(__clang__) && defined(__GNUC__) && defined(__GNUC_MINOR__)
#if __GNUC__ >= 5 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 4)
/* enable auto-vectorizer */
#pragma GCC optimize("tree-vectorize")
/* float associativity required to vectorize reductions */
#pragma GCC optimize("unsafe-math-optimizations")
/* maybe 5% gain, manual unrolling with more accumulators would be better */
#pragma GCC optimize("unroll-loops")
#endif
#endif
#include <Python.h>
#include <numpy/arrayobject.h>

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "distance_impl.h"

#define DEFINE_WRAP_CDIST(name, type)                                   \
    static PyObject *                                                   \
    cdist_ ## name ## _ ## type ## _wrap(PyObject *self, PyObject *args)\
    {                                                                   \
        PyArrayObject *XA_, *XB_, *dm_;                                 \
        Py_ssize_t mA, mB, n;                                           \
        double *dm;                                                     \
        const type *XA, *XB;                                            \
        if (!PyArg_ParseTuple(args, "O!O!O!",                           \
                              &PyArray_Type, &XA_, &PyArray_Type, &XB_, \
                              &PyArray_Type, &dm_)) {                   \
            return NULL;                                                \
        }                                                               \
        else {                                                          \
            NPY_BEGIN_ALLOW_THREADS;                                    \
            XA = (const type *)PyArray_DATA(XA_);                               \
            XB = (const type *)PyArray_DATA(XB_);                               \
            dm = (double *)PyArray_DATA(dm_);                                   \
            mA = PyArray_DIMS(XA_)[0];                                    \
            mB = PyArray_DIMS(XB_)[0];                                    \
            n = PyArray_DIMS(XA_)[1];                                     \
            cdist_ ## name ## _ ## type(XA, XB, dm, mA, mB, n);         \
            NPY_END_ALLOW_THREADS;                                      \
        }                                                               \
        return Py_BuildValue("d", 0.);                                  \
    }

DEFINE_WRAP_CDIST(bray_curtis, double)
DEFINE_WRAP_CDIST(canberra, double)
DEFINE_WRAP_CDIST(chebyshev, double)
DEFINE_WRAP_CDIST(city_block, double)
DEFINE_WRAP_CDIST(euclidean, double)
DEFINE_WRAP_CDIST(jaccard, double)
DEFINE_WRAP_CDIST(jensenshannon, double)
DEFINE_WRAP_CDIST(sqeuclidean, double)

DEFINE_WRAP_CDIST(dice, char)
DEFINE_WRAP_CDIST(jaccard, char)
DEFINE_WRAP_CDIST(kulczynski1, char)
DEFINE_WRAP_CDIST(rogerstanimoto, char)
DEFINE_WRAP_CDIST(russellrao, char)
DEFINE_WRAP_CDIST(sokalmichener, char)
DEFINE_WRAP_CDIST(sokalsneath, char)
DEFINE_WRAP_CDIST(yule, char)

static PyObject *cdist_hamming_double_wrap(
                            PyObject *self, PyObject *args, PyObject *kwargs) 
{
  PyArrayObject *XA_, *XB_, *dm_, *w_;
  int mA, mB, n;
  double *dm;
  const double *XA, *XB, *w;
  static char *kwlist[] = {"XA", "XB", "dm", "w", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, 
            "O!O!O!O!:cdist_hamming_double_wrap", kwlist,
            &PyArray_Type, &XA_, &PyArray_Type, &XB_, 
            &PyArray_Type, &dm_,
            &PyArray_Type, &w_)) {
    return 0;
  }
  else {
    NPY_BEGIN_ALLOW_THREADS;
    XA = (const double*)PyArray_DATA(XA_);
    XB = (const double*)PyArray_DATA(XB_);
    w = (const double*)PyArray_DATA(w_);
    dm = (double*)PyArray_DATA(dm_);
    mA = PyArray_DIMS(XA_)[0];
    mB = PyArray_DIMS(XB_)[0];
    n = PyArray_DIMS(XA_)[1];
    cdist_hamming_double(XA, XB, dm, mA, mB, n, w);
    NPY_END_ALLOW_THREADS;
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *cdist_hamming_char_wrap(
                            PyObject *self, PyObject *args, PyObject *kwargs) 
{
  PyArrayObject *XA_, *XB_, *dm_, *w_;
  int mA, mB, n;
  double *dm;
  const char *XA, *XB;
  const double *w;
  static char *kwlist[] = {"XA", "XB", "dm", "w", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, 
            "O!O!O!O!:cdist_hamming_char_wrap", kwlist,
            &PyArray_Type, &XA_, &PyArray_Type, &XB_, 
            &PyArray_Type, &dm_,
            &PyArray_Type, &w_)) {
    return 0;
  }
  else {
    NPY_BEGIN_ALLOW_THREADS;
    XA = (const char*)PyArray_DATA(XA_);
    XB = (const char*)PyArray_DATA(XB_);
    w = (const double*)PyArray_DATA(w_);
    dm = (double*)PyArray_DATA(dm_);
    mA = PyArray_DIMS(XA_)[0];
    mB = PyArray_DIMS(XB_)[0];
    n = PyArray_DIMS(XA_)[1];
    cdist_hamming_char(XA, XB, dm, mA, mB, n, w);
    NPY_END_ALLOW_THREADS;
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *cdist_cosine_double_wrap(PyObject *self, PyObject *args, 
                                               PyObject *kwargs) {
  PyArrayObject *XA_, *XB_, *dm_;
  int mA, mB, n, status;
  double *dm;
  const double *XA, *XB;
  static char *kwlist[] = {"XA", "XB", "dm", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, 
            "O!O!O!:cdist_cosine_double_wrap", kwlist,
            &PyArray_Type, &XA_, &PyArray_Type, &XB_, 
            &PyArray_Type, &dm_)) 
  {
    return 0;
  }
  else {
    NPY_BEGIN_THREADS_DEF;
    NPY_BEGIN_THREADS;
    XA = (const double*)PyArray_DATA(XA_);
    XB = (const double*)PyArray_DATA(XB_);
    dm = (double*)PyArray_DATA(dm_);
    mA = PyArray_DIMS(XA_)[0];
    mB = PyArray_DIMS(XB_)[0];
    n = PyArray_DIMS(XA_)[1];

    status = cdist_cosine(XA, XB, dm, mA, mB, n);
    NPY_END_THREADS;
    if(status < 0)
        return PyErr_NoMemory();
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *cdist_mahalanobis_double_wrap(PyObject *self, PyObject *args, 
                                               PyObject *kwargs) {
  PyArrayObject *XA_, *XB_, *covinv_, *dm_;
  int mA, mB, n, status;
  double *dm;
  const double *XA, *XB;
  const double *covinv;
  static char *kwlist[] = {"XA", "XB", "dm", "VI", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, 
            "O!O!O!O!:cdist_mahalanobis_double_wrap", kwlist,
            &PyArray_Type, &XA_, &PyArray_Type, &XB_, 
            &PyArray_Type, &dm_, &PyArray_Type, &covinv_)) 
  {
    return 0;
  }
  else {
    NPY_BEGIN_THREADS_DEF;
    NPY_BEGIN_THREADS;
    XA = (const double*)PyArray_DATA(XA_);
    XB = (const double*)PyArray_DATA(XB_);
    covinv = (const double*)PyArray_DATA(covinv_);
    dm = (double*)PyArray_DATA(dm_);
    mA = PyArray_DIMS(XA_)[0];
    mB = PyArray_DIMS(XB_)[0];
    n = PyArray_DIMS(XA_)[1];

    status = cdist_mahalanobis(XA, XB, dm, mA, mB, n, covinv);
    NPY_END_THREADS;
    if(status < 0)
        return PyErr_NoMemory();
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *cdist_minkowski_double_wrap(PyObject *self, PyObject *args, 
                                             PyObject *kwargs) 
{
  PyArrayObject *XA_, *XB_, *dm_;
  int mA, mB, n;
  double *dm;
  const double *XA, *XB;
  double p;
  static char *kwlist[] = {"XA", "XB", "dm", "p", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, 
            "O!O!O!d:cdist_minkowski_double_wrap", kwlist,
            &PyArray_Type, &XA_, &PyArray_Type, &XB_, 
            &PyArray_Type, &dm_,
            &p)) {
    return 0;
  }
  else {
    NPY_BEGIN_ALLOW_THREADS;
    XA = (const double*)PyArray_DATA(XA_);
    XB = (const double*)PyArray_DATA(XB_);
    dm = (double*)PyArray_DATA(dm_);
    mA = PyArray_DIMS(XA_)[0];
    mB = PyArray_DIMS(XB_)[0];
    n = PyArray_DIMS(XA_)[1];
    cdist_minkowski(XA, XB, dm, mA, mB, n, p);
    NPY_END_ALLOW_THREADS;
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *cdist_seuclidean_double_wrap(PyObject *self, PyObject *args, 
                                              PyObject *kwargs) 
{
  PyArrayObject *XA_, *XB_, *dm_, *var_;
  int mA, mB, n;
  double *dm;
  const double *XA, *XB, *var;
  static char *kwlist[] = {"XA", "XB", "dm", "V", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, 
            "O!O!O!O!:cdist_seuclidean_double_wrap", kwlist,
            &PyArray_Type, &XA_, &PyArray_Type, &XB_, 
            &PyArray_Type, &dm_, &PyArray_Type, &var_)) {
    return 0;
  }
  else {
    NPY_BEGIN_ALLOW_THREADS;
    XA = (const double*)PyArray_DATA(XA_);
    XB = (const double*)PyArray_DATA(XB_);
    dm = (double*)PyArray_DATA(dm_);
    var = (double*)PyArray_DATA(var_);
    mA = PyArray_DIMS(XA_)[0];
    mB = PyArray_DIMS(XB_)[0];
    n = PyArray_DIMS(XA_)[1];

    cdist_seuclidean(XA, XB, var, dm, mA, mB, n);
    NPY_END_ALLOW_THREADS;
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *cdist_weighted_chebyshev_double_wrap(
    PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyArrayObject *XA_, *XB_, *dm_, *w_;
  int mA, mB, n;
  double *dm;
  const double *XA, *XB, *w;
  static char *kwlist[] = {"XA", "XB", "dm", "w", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs,
            "O!O!O!O!:cdist_weighted_chebyshev_double_wrap", kwlist,
            &PyArray_Type, &XA_, &PyArray_Type, &XB_,
            &PyArray_Type, &dm_,
            &PyArray_Type, &w_)) {
    return 0;
  }
  else {
    NPY_BEGIN_ALLOW_THREADS;
    XA = (const double*)PyArray_DATA(XA_);
    XB = (const double*)PyArray_DATA(XB_);
    w = (const double*)PyArray_DATA(w_);
    dm = (double*)PyArray_DATA(dm_);
    mA = PyArray_DIMS(XA_)[0];
    mB = PyArray_DIMS(XB_)[0];
    n = PyArray_DIMS(XA_)[1];
    cdist_weighted_chebyshev(XA, XB, dm, mA, mB, n, w);
    NPY_END_ALLOW_THREADS;
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *cdist_weighted_minkowski_double_wrap(
                            PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyArrayObject *XA_, *XB_, *dm_, *w_;
  int mA, mB, n;
  double *dm;
  const double *XA, *XB, *w;
  double p;
  static char *kwlist[] = {"XA", "XB", "dm", "p", "w", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs,
            "O!O!O!dO!:cdist_weighted_minkowski_double_wrap", kwlist,
            &PyArray_Type, &XA_, &PyArray_Type, &XB_,
            &PyArray_Type, &dm_,
            &p,
            &PyArray_Type, &w_)) {
    return 0;
  }
  else {
    NPY_BEGIN_ALLOW_THREADS;
    XA = (const double*)PyArray_DATA(XA_);
    XB = (const double*)PyArray_DATA(XB_);
    w = (const double*)PyArray_DATA(w_);
    dm = (double*)PyArray_DATA(dm_);
    mA = PyArray_DIMS(XA_)[0];
    mB = PyArray_DIMS(XB_)[0];
    n = PyArray_DIMS(XA_)[1];
    cdist_weighted_minkowski(XA, XB, dm, mA, mB, n, p, w);
    NPY_END_ALLOW_THREADS;
  }
  return Py_BuildValue("d", 0.0);
}

/***************************** pdist ***/

#define DEFINE_WRAP_PDIST(name, type)                                   \
    static PyObject *                                                   \
    pdist_ ## name ## _ ## type ## _wrap(PyObject *self, PyObject *args)\
    {                                                                   \
        PyArrayObject *X_, *dm_;                                        \
        Py_ssize_t m, n;                                                \
        double *dm;                                                     \
        const type *X;                                                  \
        if (!PyArg_ParseTuple(args, "O!O!", &PyArray_Type, &X_,         \
                                            &PyArray_Type, &dm_)) {     \
            return NULL;                                                \
        }                                                               \
        else {                                                          \
            NPY_BEGIN_ALLOW_THREADS;                                    \
            X = (const type *)PyArray_DATA(X_);                                 \
            dm = (double *)PyArray_DATA(dm_);                                   \
            m = PyArray_DIMS(X_)[0];                                      \
            n = PyArray_DIMS(X_)[1];                                      \
            pdist_ ## name ## _ ## type(X, dm, m, n);                   \
            NPY_END_ALLOW_THREADS;                                      \
        }                                                               \
        return Py_BuildValue("d", 0.);                                  \
    }

DEFINE_WRAP_PDIST(bray_curtis, double)
DEFINE_WRAP_PDIST(canberra, double)
DEFINE_WRAP_PDIST(chebyshev, double)
DEFINE_WRAP_PDIST(city_block, double)
DEFINE_WRAP_PDIST(euclidean, double)
DEFINE_WRAP_PDIST(jaccard, double)
DEFINE_WRAP_PDIST(jensenshannon, double)
DEFINE_WRAP_PDIST(sqeuclidean, double)

DEFINE_WRAP_PDIST(dice, char)
DEFINE_WRAP_PDIST(kulczynski1, char)
DEFINE_WRAP_PDIST(jaccard, char)
DEFINE_WRAP_PDIST(rogerstanimoto, char)
DEFINE_WRAP_PDIST(russellrao, char)
DEFINE_WRAP_PDIST(sokalmichener, char)
DEFINE_WRAP_PDIST(sokalsneath, char)
DEFINE_WRAP_PDIST(yule, char)

static PyObject *pdist_hamming_double_wrap(
                            PyObject *self, PyObject *args, PyObject *kwargs) 
{
  PyArrayObject *X_, *dm_, *w_;
  int m, n;
  double *dm;
  const double *X, *w;
  static char *kwlist[] = {"X", "dm", "w", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, 
            "O!O!O!:pdist_hamming_double_wrap", kwlist,
            &PyArray_Type, &X_,
            &PyArray_Type, &dm_,
            &PyArray_Type, &w_)) {
    return 0;
  }
  else {
    NPY_BEGIN_ALLOW_THREADS;
    X = (const double*)PyArray_DATA(X_);
    dm = (double*)PyArray_DATA(dm_);
    w = (const double*)PyArray_DATA(w_);
    m = PyArray_DIMS(X_)[0];
    n = PyArray_DIMS(X_)[1];

    pdist_hamming_double(X, dm, m, n, w);
    NPY_END_ALLOW_THREADS;
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *pdist_hamming_char_wrap(
                            PyObject *self, PyObject *args, PyObject *kwargs) 
{
  PyArrayObject *X_, *dm_, *w_;
  int m, n;
  const char *X;
  const double *w;
  double *dm;
  static char *kwlist[] = {"X", "dm", "w", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, 
            "O!O!O!:pdist_hamming_char_wrap", kwlist,
            &PyArray_Type, &X_,
            &PyArray_Type, &dm_,
            &PyArray_Type, &w_)) {
    return 0;
  }
  else {
    NPY_BEGIN_ALLOW_THREADS;
    X = (const char*)PyArray_DATA(X_);
    dm = (double*)PyArray_DATA(dm_);
    w = (const double*)PyArray_DATA(w_);
    m = PyArray_DIMS(X_)[0];
    n = PyArray_DIMS(X_)[1];

    pdist_hamming_char(X, dm, m, n, w);
    NPY_END_ALLOW_THREADS;
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *pdist_cosine_double_wrap(PyObject *self, PyObject *args, 
                                          PyObject *kwargs) 
{
  PyArrayObject *X_, *dm_;
  int m, n, status;
  double *dm;
  const double *X;
  static char *kwlist[] = {"X", "dm", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, 
            "O!O!:pdist_cosine_double_wrap", kwlist,
            &PyArray_Type, &X_,
            &PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    NPY_BEGIN_THREADS_DEF;
    NPY_BEGIN_THREADS;
    X = (const double*)PyArray_DATA(X_);
    dm = (double*)PyArray_DATA(dm_);
    m = PyArray_DIMS(X_)[0];
    n = PyArray_DIMS(X_)[1];
    
    status = pdist_cosine(X, dm, m, n);
    NPY_END_THREADS;
    if(status < 0)
        return PyErr_NoMemory();
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *pdist_mahalanobis_double_wrap(PyObject *self, PyObject *args, 
                                               PyObject *kwargs) {
  PyArrayObject *X_, *covinv_, *dm_;
  int m, n, status;
  double *dm;
  const double *X;
  const double *covinv;
  static char *kwlist[] = {"X", "dm", "VI", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, 
            "O!O!O!:pdist_mahalanobis_double_wrap", kwlist,
            &PyArray_Type, &X_,
            &PyArray_Type, &dm_, 
            &PyArray_Type, &covinv_)) {
    return 0;
  }
  else {
    NPY_BEGIN_THREADS_DEF;
    NPY_BEGIN_THREADS;
    X = (const double*)PyArray_DATA(X_);
    covinv = (const double*)PyArray_DATA(covinv_);
    dm = (double*)PyArray_DATA(dm_);
    m = PyArray_DIMS(X_)[0];
    n = PyArray_DIMS(X_)[1];

    status = pdist_mahalanobis(X, dm, m, n, covinv);
    NPY_END_THREADS;
    if(status < 0)
        return PyErr_NoMemory();
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *pdist_minkowski_double_wrap(PyObject *self, PyObject *args, 
                                             PyObject *kwargs) 
{
  PyArrayObject *X_, *dm_;
  int m, n;
  double *dm, *X;
  double p;
  static char *kwlist[] = {"X", "dm", "p", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, 
            "O!O!d:pdist_minkowski_double_wrap", kwlist,
            &PyArray_Type, &X_,
            &PyArray_Type, &dm_,
            &p)) {
    return 0;
  }
  else {
    NPY_BEGIN_ALLOW_THREADS;
    X = (double*)PyArray_DATA(X_);
    dm = (double*)PyArray_DATA(dm_);
    m = PyArray_DIMS(X_)[0];
    n = PyArray_DIMS(X_)[1];

    pdist_minkowski(X, dm, m, n, p);
    NPY_END_ALLOW_THREADS;
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *pdist_seuclidean_double_wrap(PyObject *self, PyObject *args, 
                                              PyObject *kwargs) 
{
  PyArrayObject *X_, *dm_, *var_;
  int m, n;
  double *dm;
  const double *X, *var;
  static char *kwlist[] = {"X", "dm", "V", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, 
            "O!O!O!:pdist_seuclidean_double_wrap", kwlist,
            &PyArray_Type, &X_,
            &PyArray_Type, &dm_,
            &PyArray_Type, &var_)) {
    return 0;
  }
  else {
    NPY_BEGIN_ALLOW_THREADS;
    X = (double*)PyArray_DATA(X_);
    dm = (double*)PyArray_DATA(dm_);
    var = (double*)PyArray_DATA(var_);
    m = PyArray_DIMS(X_)[0];
    n = PyArray_DIMS(X_)[1];

    pdist_seuclidean(X, var, dm, m, n);
    NPY_END_ALLOW_THREADS;
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *pdist_weighted_chebyshev_double_wrap(
  PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyArrayObject *X_, *dm_, *w_;
  int m, n;
  double *dm, *X, *w;
  static char *kwlist[] = {"X", "dm", "w", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs,
                                   "O!O!O!:pdist_weighted_minkowski_double_wrap", kwlist,
                                   &PyArray_Type, &X_,
                                   &PyArray_Type, &dm_,
                                   &PyArray_Type, &w_)) {
    return 0;
  }
  else {
    NPY_BEGIN_ALLOW_THREADS;
    X = (double*)PyArray_DATA(X_);
    dm = (double*)PyArray_DATA(dm_);
    w = (double*)PyArray_DATA(w_);
    m = PyArray_DIMS(X_)[0];
    n = PyArray_DIMS(X_)[1];

    pdist_weighted_chebyshev(X, dm, m, n, w);
    NPY_END_ALLOW_THREADS;
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *pdist_weighted_minkowski_double_wrap(
                            PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyArrayObject *X_, *dm_, *w_;
  int m, n;
  double *dm, *X, *w;
  double p;
  static char *kwlist[] = {"X", "dm", "p", "w", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs,
            "O!O!dO!:pdist_weighted_minkowski_double_wrap", kwlist,
            &PyArray_Type, &X_,
            &PyArray_Type, &dm_,
            &p,
            &PyArray_Type, &w_)) {
    return 0;
  }
  else {
    NPY_BEGIN_ALLOW_THREADS;
    X = (double*)PyArray_DATA(X_);
    dm = (double*)PyArray_DATA(dm_);
    w = (double*)PyArray_DATA(w_);
    m = PyArray_DIMS(X_)[0];
    n = PyArray_DIMS(X_)[1];

    pdist_weighted_minkowski(X, dm, m, n, p, w);
    NPY_END_ALLOW_THREADS;
  }
  return Py_BuildValue("d", 0.0);
}


static PyObject *to_squareform_from_vector_wrap(PyObject *self, PyObject *args) 
{
  PyArrayObject *M_, *v_;
  int n, elsize;
  if (!PyArg_ParseTuple(args, "O!O!",
            &PyArray_Type, &M_,
            &PyArray_Type, &v_)) {
    return 0;
  }
  NPY_BEGIN_ALLOW_THREADS;
  n = PyArray_DIMS(M_)[0];
  elsize = PyArray_DESCR(M_)->elsize;
  if (elsize == 8) {
    dist_to_squareform_from_vector_double(
        (double*)PyArray_DATA(M_), (const double*)PyArray_DATA(v_), n);
  } else {
    dist_to_squareform_from_vector_generic(
        (char*)PyArray_DATA(M_), (const char*)PyArray_DATA(v_), n, elsize);
  }
  NPY_END_ALLOW_THREADS;
  return Py_BuildValue("");
}

static PyObject *to_vector_from_squareform_wrap(PyObject *self, PyObject *args) 
{
  PyArrayObject *M_, *v_;
  int n, s;
  char *v;
  const char *M;
  if (!PyArg_ParseTuple(args, "O!O!",
            &PyArray_Type, &M_,
            &PyArray_Type, &v_)) {
    return 0;
  }
  else {
    NPY_BEGIN_ALLOW_THREADS;
    M = (const char*)PyArray_DATA(M_);
    v = (char*)PyArray_DATA(v_);
    n = PyArray_DIMS(M_)[0];
    s = PyArray_DESCR(M_)->elsize;
    dist_to_vector_from_squareform(M, v, n, s);
    NPY_END_ALLOW_THREADS;
  }
  return Py_BuildValue("");
}


static PyMethodDef _distanceWrapMethods[] = {
  {"cdist_braycurtis_double_wrap",
   cdist_bray_curtis_double_wrap,
   METH_VARARGS},
  {"cdist_canberra_double_wrap",
   cdist_canberra_double_wrap,
   METH_VARARGS},
  {"cdist_chebyshev_double_wrap",
   cdist_chebyshev_double_wrap,
   METH_VARARGS},
  {"cdist_cityblock_double_wrap",
   cdist_city_block_double_wrap,
   METH_VARARGS},
  {"cdist_cosine_double_wrap",
   (PyCFunction) cdist_cosine_double_wrap,
   METH_VARARGS | METH_KEYWORDS},
  {"cdist_dice_bool_wrap",
   cdist_dice_char_wrap,
   METH_VARARGS},
  {"cdist_euclidean_double_wrap",
   cdist_euclidean_double_wrap,
   METH_VARARGS},
  {"cdist_sqeuclidean_double_wrap",
   cdist_sqeuclidean_double_wrap,
   METH_VARARGS},
  {"cdist_hamming_double_wrap",
   (PyCFunction) cdist_hamming_double_wrap,
   METH_VARARGS | METH_KEYWORDS},
  {"cdist_hamming_bool_wrap",
   (PyCFunction) cdist_hamming_char_wrap,
   METH_VARARGS | METH_KEYWORDS},
  {"cdist_jaccard_double_wrap",
   cdist_jaccard_double_wrap,
   METH_VARARGS},
  {"cdist_jaccard_bool_wrap",
   cdist_jaccard_char_wrap,
   METH_VARARGS},
  {"cdist_jensenshannon_double_wrap",
   cdist_jensenshannon_double_wrap,
   METH_VARARGS},
   {"cdist_kulczynski1_bool_wrap",
   cdist_kulczynski1_char_wrap,
   METH_VARARGS},
  {"cdist_mahalanobis_double_wrap",
   (PyCFunction) cdist_mahalanobis_double_wrap,
   METH_VARARGS | METH_KEYWORDS},
  {"cdist_minkowski_double_wrap",
   (PyCFunction) cdist_minkowski_double_wrap,
   METH_VARARGS | METH_KEYWORDS},
  {"cdist_weighted_chebyshev_double_wrap",
   (PyCFunction) cdist_weighted_chebyshev_double_wrap,
   METH_VARARGS | METH_KEYWORDS},
  {"cdist_weighted_minkowski_double_wrap",
   (PyCFunction) cdist_weighted_minkowski_double_wrap,
   METH_VARARGS | METH_KEYWORDS},
  {"cdist_rogerstanimoto_bool_wrap",
   cdist_rogerstanimoto_char_wrap,
   METH_VARARGS},
  {"cdist_russellrao_bool_wrap",
   cdist_russellrao_char_wrap,
   METH_VARARGS},
  {"cdist_seuclidean_double_wrap",
   (PyCFunction) cdist_seuclidean_double_wrap,
   METH_VARARGS | METH_KEYWORDS},
  {"cdist_sokalmichener_bool_wrap",
   cdist_sokalmichener_char_wrap,
   METH_VARARGS},
  {"cdist_sokalsneath_bool_wrap",
   cdist_sokalsneath_char_wrap,
   METH_VARARGS},
  {"cdist_yule_bool_wrap",
   cdist_yule_char_wrap,
   METH_VARARGS},
  {"pdist_braycurtis_double_wrap",
   pdist_bray_curtis_double_wrap,
   METH_VARARGS},
  {"pdist_canberra_double_wrap",
   pdist_canberra_double_wrap,
   METH_VARARGS},
  {"pdist_chebyshev_double_wrap",
   pdist_chebyshev_double_wrap,
   METH_VARARGS},
  {"pdist_cityblock_double_wrap",
   pdist_city_block_double_wrap,
   METH_VARARGS},
  {"pdist_cosine_double_wrap",
   (PyCFunction) pdist_cosine_double_wrap,
   METH_VARARGS | METH_KEYWORDS},
  {"pdist_dice_bool_wrap",
   pdist_dice_char_wrap,
   METH_VARARGS},
  {"pdist_euclidean_double_wrap",
   pdist_euclidean_double_wrap,
   METH_VARARGS},
  {"pdist_sqeuclidean_double_wrap",
   pdist_sqeuclidean_double_wrap,
   METH_VARARGS},
  {"pdist_hamming_double_wrap",
   (PyCFunction) pdist_hamming_double_wrap,
   METH_VARARGS | METH_KEYWORDS},
  {"pdist_hamming_bool_wrap",
   (PyCFunction) pdist_hamming_char_wrap,
   METH_VARARGS | METH_KEYWORDS},
  {"pdist_jaccard_double_wrap",
   pdist_jaccard_double_wrap,
   METH_VARARGS},
  {"pdist_jaccard_bool_wrap",
   pdist_jaccard_char_wrap,
   METH_VARARGS},
  {"pdist_jensenshannon_double_wrap",
   pdist_jensenshannon_double_wrap,
   METH_VARARGS},
   {"pdist_kulczynski1_bool_wrap",
   pdist_kulczynski1_char_wrap,
   METH_VARARGS},
  {"pdist_mahalanobis_double_wrap",
   (PyCFunction) pdist_mahalanobis_double_wrap,
   METH_VARARGS | METH_KEYWORDS},
  {"pdist_minkowski_double_wrap",
   (PyCFunction) pdist_minkowski_double_wrap,
   METH_VARARGS | METH_KEYWORDS},
  {"pdist_weighted_chebyshev_double_wrap",
   (PyCFunction) pdist_weighted_chebyshev_double_wrap,
   METH_VARARGS | METH_KEYWORDS},
  {"pdist_weighted_minkowski_double_wrap",
   (PyCFunction) pdist_weighted_minkowski_double_wrap,
   METH_VARARGS | METH_KEYWORDS},
  {"pdist_rogerstanimoto_bool_wrap",
   pdist_rogerstanimoto_char_wrap,
   METH_VARARGS},
  {"pdist_russellrao_bool_wrap",
   pdist_russellrao_char_wrap,
   METH_VARARGS},
  {"pdist_seuclidean_double_wrap",
   (PyCFunction) pdist_seuclidean_double_wrap,
   METH_VARARGS | METH_KEYWORDS},
  {"pdist_sokalmichener_bool_wrap",
   pdist_sokalmichener_char_wrap,
   METH_VARARGS},
  {"pdist_sokalsneath_bool_wrap",
   pdist_sokalsneath_char_wrap,
   METH_VARARGS},
  {"pdist_yule_bool_wrap",
   pdist_yule_char_wrap,
   METH_VARARGS},
  {"to_squareform_from_vector_wrap",
   to_squareform_from_vector_wrap,
   METH_VARARGS},
  {"to_vector_from_squareform_wrap",
   to_vector_from_squareform_wrap,
   METH_VARARGS},
  {NULL, NULL}     /* Sentinel - marks the end of this structure */
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_distance_wrap",
    NULL,
    -1,
    _distanceWrapMethods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit__distance_wrap(void)
{
    import_array();
    return PyModule_Create(&moduledef);
}
