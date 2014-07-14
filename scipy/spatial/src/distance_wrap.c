/**
 * distance_wrap.c
 *
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

#include <Python.h>
#include <numpy/arrayobject.h>

#include "distance.h"

static PyObject *cdist_euclidean_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *XA_, *XB_, *dm_;
  int mA, mB, n;
  double *dm;
  const double *XA, *XB;
  if (!PyArg_ParseTuple(args, "O!O!O!",
			&PyArray_Type, &XA_, &PyArray_Type, &XB_, 
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    XA = (const double*)XA_->data;
    XB = (const double*)XB_->data;
    dm = (double*)dm_->data;
    mA = XA_->dimensions[0];
    mB = XB_->dimensions[0];
    n = XA_->dimensions[1];

    cdist_euclidean(XA, XB, dm, mA, mB, n);
  }
  return Py_BuildValue("d", 0.0);
}


static PyObject *cdist_sqeuclidean_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *XA_, *XB_, *dm_;
  int mA, mB, n;
  double *dm;
  const double *XA, *XB;
  if (!PyArg_ParseTuple(args, "O!O!O!",
			&PyArray_Type, &XA_, &PyArray_Type, &XB_, 
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    XA = (const double*)XA_->data;
    XB = (const double*)XB_->data;
    dm = (double*)dm_->data;
    mA = XA_->dimensions[0];
    mB = XB_->dimensions[0];
    n = XA_->dimensions[1];

    cdist_sqeuclidean(XA, XB, dm, mA, mB, n);
  }
  return Py_BuildValue("d", 0.0);
}


static PyObject *cdist_canberra_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *XA_, *XB_, *dm_;
  int mA, mB, n;
  double *dm;
  const double *XA, *XB;
  if (!PyArg_ParseTuple(args, "O!O!O!",
			&PyArray_Type, &XA_, &PyArray_Type, &XB_, 
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    XA = (const double*)XA_->data;
    XB = (const double*)XB_->data;
    dm = (double*)dm_->data;
    mA = XA_->dimensions[0];
    mB = XB_->dimensions[0];
    n = XA_->dimensions[1];

    cdist_canberra(XA, XB, dm, mA, mB, n);
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *cdist_bray_curtis_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *XA_, *XB_, *dm_;
  int mA, mB, n;
  double *dm;
  const double *XA, *XB;
  if (!PyArg_ParseTuple(args, "O!O!O!",
			&PyArray_Type, &XA_, &PyArray_Type, &XB_, 
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    XA = (const double*)XA_->data;
    XB = (const double*)XB_->data;
    dm = (double*)dm_->data;
    mA = XA_->dimensions[0];
    mB = XB_->dimensions[0];
    n = XA_->dimensions[1];

    cdist_bray_curtis(XA, XB, dm, mA, mB, n);
  }
  return Py_BuildValue("d", 0.0);
}


static PyObject *cdist_mahalanobis_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *XA_, *XB_, *covinv_, *dm_;
  int mA, mB, n;
  double *dm;
  const double *XA, *XB;
  const double *covinv;
  if (!PyArg_ParseTuple(args, "O!O!O!O!",
			&PyArray_Type, &XA_, &PyArray_Type, &XB_, 
			&PyArray_Type, &covinv_,
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    XA = (const double*)XA_->data;
    XB = (const double*)XB_->data;
    covinv = (const double*)covinv_->data;
    dm = (double*)dm_->data;
    mA = XA_->dimensions[0];
    mB = XB_->dimensions[0];
    n = XA_->dimensions[1];

    cdist_mahalanobis(XA, XB, covinv, dm, mA, mB, n);
  }
  return Py_BuildValue("d", 0.0);
}


static PyObject *cdist_chebyshev_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *XA_, *XB_, *dm_;
  int mA, mB, n;
  double *dm;
  const double *XA, *XB;
  if (!PyArg_ParseTuple(args, "O!O!O!",
			&PyArray_Type, &XA_, &PyArray_Type, &XB_, 
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    XA = (const double*)XA_->data;
    XB = (const double*)XB_->data;
    dm = (double*)dm_->data;
    mA = XA_->dimensions[0];
    mB = XB_->dimensions[0];
    n = XA_->dimensions[1];

    cdist_chebyshev(XA, XB, dm, mA, mB, n);
  }
  return Py_BuildValue("d", 0.0);
}


static PyObject *cdist_cosine_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *XA_, *XB_, *dm_, *normsA_, *normsB_;
  int mA, mB, n;
  double *dm;
  const double *XA, *XB, *normsA, *normsB;
  if (!PyArg_ParseTuple(args, "O!O!O!O!O!",
			&PyArray_Type, &XA_, &PyArray_Type, &XB_, 
			&PyArray_Type, &dm_,
			&PyArray_Type, &normsA_,
			&PyArray_Type, &normsB_)) {
    return 0;
  }
  else {
    XA = (const double*)XA_->data;
    XB = (const double*)XB_->data;
    dm = (double*)dm_->data;
    normsA = (const double*)normsA_->data;
    normsB = (const double*)normsB_->data;
    mA = XA_->dimensions[0];
    mB = XB_->dimensions[0];
    n = XA_->dimensions[1];

    cdist_cosine(XA, XB, dm, mA, mB, n, normsA, normsB);
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *cdist_seuclidean_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *XA_, *XB_, *dm_, *var_;
  int mA, mB, n;
  double *dm;
  const double *XA, *XB, *var;
  if (!PyArg_ParseTuple(args, "O!O!O!O!",
			&PyArray_Type, &XA_, &PyArray_Type, &XB_, 
			&PyArray_Type, &var_,
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    XA = (const double*)XA_->data;
    XB = (const double*)XB_->data;
    dm = (double*)dm_->data;
    var = (double*)var_->data;
    mA = XA_->dimensions[0];
    mB = XB_->dimensions[0];
    n = XA_->dimensions[1];

    cdist_seuclidean(XA, XB, var, dm, mA, mB, n);
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *cdist_city_block_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *XA_, *XB_, *dm_;
  int mA, mB, n;
  double *dm;
  const double *XA, *XB;
  if (!PyArg_ParseTuple(args, "O!O!O!",
			&PyArray_Type, &XA_, &PyArray_Type, &XB_, 
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    XA = (const double*)XA_->data;
    XB = (const double*)XB_->data;
    dm = (double*)dm_->data;
    mA = XA_->dimensions[0];
    mB = XB_->dimensions[0];
    n = XA_->dimensions[1];

    cdist_city_block(XA, XB, dm, mA, mB, n);
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *cdist_hamming_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *XA_, *XB_, *dm_;
  int mA, mB, n;
  double *dm;
  const double *XA, *XB;
  if (!PyArg_ParseTuple(args, "O!O!O!",
			&PyArray_Type, &XA_, &PyArray_Type, &XB_, 
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    XA = (const double*)XA_->data;
    XB = (const double*)XB_->data;
    dm = (double*)dm_->data;
    mA = XA_->dimensions[0];
    mB = XB_->dimensions[0];
    n = XA_->dimensions[1];

    cdist_hamming(XA, XB, dm, mA, mB, n);
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *cdist_hamming_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *XA_, *XB_, *dm_;
  int mA, mB, n;
  double *dm;
  const char *XA, *XB;
  if (!PyArg_ParseTuple(args, "O!O!O!",
			&PyArray_Type, &XA_, &PyArray_Type, &XB_, 
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    XA = (const char*)XA_->data;
    XB = (const char*)XB_->data;
    dm = (double*)dm_->data;
    mA = XA_->dimensions[0];
    mB = XB_->dimensions[0];
    n = XA_->dimensions[1];

    cdist_hamming_bool(XA, XB, dm, mA, mB, n);
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *cdist_jaccard_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *XA_, *XB_, *dm_;
  int mA, mB, n;
  double *dm;
  const double *XA, *XB;
  if (!PyArg_ParseTuple(args, "O!O!O!",
			&PyArray_Type, &XA_, &PyArray_Type, &XB_, 
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    XA = (const double*)XA_->data;
    XB = (const double*)XB_->data;
    dm = (double*)dm_->data;
    mA = XA_->dimensions[0];
    mB = XB_->dimensions[0];
    n = XA_->dimensions[1];

    cdist_jaccard(XA, XB, dm, mA, mB, n);
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *cdist_jaccard_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *XA_, *XB_, *dm_;
  int mA, mB, n;
  double *dm;
  const char *XA, *XB;
  if (!PyArg_ParseTuple(args, "O!O!O!",
			&PyArray_Type, &XA_, &PyArray_Type, &XB_, 
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    XA = (const char*)XA_->data;
    XB = (const char*)XB_->data;
    dm = (double*)dm_->data;
    mA = XA_->dimensions[0];
    mB = XB_->dimensions[0];
    n = XA_->dimensions[1];

    cdist_jaccard_bool(XA, XB, dm, mA, mB, n);
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *cdist_minkowski_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *XA_, *XB_, *dm_;
  int mA, mB, n;
  double *dm;
  const double *XA, *XB;
  double p;
  if (!PyArg_ParseTuple(args, "O!O!O!d",
			&PyArray_Type, &XA_, &PyArray_Type, &XB_, 
			&PyArray_Type, &dm_,
			&p)) {
    return 0;
  }
  else {
    XA = (const double*)XA_->data;
    XB = (const double*)XB_->data;
    dm = (double*)dm_->data;
    mA = XA_->dimensions[0];
    mB = XB_->dimensions[0];
    n = XA_->dimensions[1];
    cdist_minkowski(XA, XB, dm, mA, mB, n, p);
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *cdist_weighted_minkowski_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *XA_, *XB_, *dm_, *w_;
  int mA, mB, n;
  double *dm;
  const double *XA, *XB, *w;
  double p;
  if (!PyArg_ParseTuple(args, "O!O!O!dO!",
			&PyArray_Type, &XA_, &PyArray_Type, &XB_, 
			&PyArray_Type, &dm_,
			&p,
			&PyArray_Type, &w_)) {
    return 0;
  }
  else {
    XA = (const double*)XA_->data;
    XB = (const double*)XB_->data;
    w = (const double*)w_->data;
    dm = (double*)dm_->data;
    mA = XA_->dimensions[0];
    mB = XB_->dimensions[0];
    n = XA_->dimensions[1];
    cdist_weighted_minkowski(XA, XB, dm, mA, mB, n, p, w);
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *cdist_yule_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *XA_, *XB_, *dm_;
  int mA, mB, n;
  double *dm;
  const char *XA, *XB;
  if (!PyArg_ParseTuple(args, "O!O!O!",
			&PyArray_Type, &XA_, &PyArray_Type, &XB_, 
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    XA = (const char*)XA_->data;
    XB = (const char*)XB_->data;
    dm = (double*)dm_->data;
    mA = XA_->dimensions[0];
    mB = XB_->dimensions[0];
    n = XA_->dimensions[1];

    cdist_yule_bool(XA, XB, dm, mA, mB, n);
  }
  return Py_BuildValue("");
}

static PyObject *cdist_matching_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *XA_, *XB_, *dm_;
  int mA, mB, n;
  double *dm;
  const char *XA, *XB;
  if (!PyArg_ParseTuple(args, "O!O!O!",
			&PyArray_Type, &XA_, &PyArray_Type, &XB_, 
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    XA = (const char*)XA_->data;
    XB = (const char*)XB_->data;
    dm = (double*)dm_->data;
    mA = XA_->dimensions[0];
    mB = XB_->dimensions[0];
    n = XA_->dimensions[1];

    cdist_matching_bool(XA, XB, dm, mA, mB, n);
  }
  return Py_BuildValue("");
}

static PyObject *cdist_dice_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *XA_, *XB_, *dm_;
  int mA, mB, n;
  double *dm;
  const char *XA, *XB;
  if (!PyArg_ParseTuple(args, "O!O!O!",
			&PyArray_Type, &XA_, &PyArray_Type, &XB_, 
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    XA = (const char*)XA_->data;
    XB = (const char*)XB_->data;
    dm = (double*)dm_->data;
    mA = XA_->dimensions[0];
    mB = XB_->dimensions[0];
    n = XA_->dimensions[1];

    cdist_dice_bool(XA, XB, dm, mA, mB, n);
  }
  return Py_BuildValue("");
}

static PyObject *cdist_rogerstanimoto_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *XA_, *XB_, *dm_;
  int mA, mB, n;
  double *dm;
  const char *XA, *XB;
  if (!PyArg_ParseTuple(args, "O!O!O!",
			&PyArray_Type, &XA_, &PyArray_Type, &XB_, 
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    XA = (const char*)XA_->data;
    XB = (const char*)XB_->data;
    dm = (double*)dm_->data;
    mA = XA_->dimensions[0];
    mB = XB_->dimensions[0];
    n = XA_->dimensions[1];

    cdist_rogerstanimoto_bool(XA, XB, dm, mA, mB, n);
  }
  return Py_BuildValue("");
}

static PyObject *cdist_russellrao_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *XA_, *XB_, *dm_;
  int mA, mB, n;
  double *dm;
  const char *XA, *XB;
  if (!PyArg_ParseTuple(args, "O!O!O!",
			&PyArray_Type, &XA_, &PyArray_Type, &XB_, 
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    XA = (const char*)XA_->data;
    XB = (const char*)XB_->data;
    dm = (double*)dm_->data;
    mA = XA_->dimensions[0];
    mB = XB_->dimensions[0];
    n = XA_->dimensions[1];

    cdist_russellrao_bool(XA, XB, dm, mA, mB, n);
  }
  return Py_BuildValue("");
}

static PyObject *cdist_kulsinski_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *XA_, *XB_, *dm_;
  int mA, mB, n;
  double *dm;
  const char *XA, *XB;
  if (!PyArg_ParseTuple(args, "O!O!O!",
			&PyArray_Type, &XA_, &PyArray_Type, &XB_, 
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    XA = (const char*)XA_->data;
    XB = (const char*)XB_->data;
    dm = (double*)dm_->data;
    mA = XA_->dimensions[0];
    mB = XB_->dimensions[0];
    n = XA_->dimensions[1];

    cdist_kulsinski_bool(XA, XB, dm, mA, mB, n);
  }
  return Py_BuildValue("");
}

static PyObject *cdist_sokalmichener_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *XA_, *XB_, *dm_;
  int mA, mB, n;
  double *dm;
  const char *XA, *XB;
  if (!PyArg_ParseTuple(args, "O!O!O!",
			&PyArray_Type, &XA_, &PyArray_Type, &XB_, 
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    XA = (const char*)XA_->data;
    XB = (const char*)XB_->data;
    dm = (double*)dm_->data;
    mA = XA_->dimensions[0];
    mB = XB_->dimensions[0];
    n = XA_->dimensions[1];

    cdist_sokalmichener_bool(XA, XB, dm, mA, mB, n);
  }
  return Py_BuildValue("");
}

static PyObject *cdist_sokalsneath_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *XA_, *XB_, *dm_;
  int mA, mB, n;
  double *dm;
  const char *XA, *XB;
  if (!PyArg_ParseTuple(args, "O!O!O!",
			&PyArray_Type, &XA_, &PyArray_Type, &XB_, 
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    XA = (const char*)XA_->data;
    XB = (const char*)XB_->data;
    dm = (double*)dm_->data;
    mA = XA_->dimensions[0];
    mB = XB_->dimensions[0];
    n = XA_->dimensions[1];

    cdist_sokalsneath_bool(XA, XB, dm, mA, mB, n);
  }
  return Py_BuildValue("");
}

/***************************** pdist ***/

static PyObject *pdist_euclidean_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *X_, *dm_;
  int m, n;
  double *dm;
  const double *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &X_,
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    X = (const double*)X_->data;
    dm = (double*)dm_->data;
    m = X_->dimensions[0];
    n = X_->dimensions[1];

    pdist_euclidean(X, dm, m, n);
  }
  return Py_BuildValue("d", 0.0);
}


static PyObject *pdist_sqeuclidean_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *X_, *dm_;
  int m, n;
  double *dm;
  const double *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &X_,
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    X = (const double*)X_->data;
    dm = (double*)dm_->data;
    m = X_->dimensions[0];
    n = X_->dimensions[1];

    pdist_sqeuclidean(X, dm, m, n);
  }
  return Py_BuildValue("d", 0.0);
}


static PyObject *pdist_canberra_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *X_, *dm_;
  int m, n;
  double *dm;
  const double *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &X_,
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    X = (const double*)X_->data;
    dm = (double*)dm_->data;
    m = X_->dimensions[0];
    n = X_->dimensions[1];

    pdist_canberra(X, dm, m, n);
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *pdist_bray_curtis_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *X_, *dm_;
  int m, n;
  double *dm;
  const double *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &X_,
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    X = (const double*)X_->data;
    dm = (double*)dm_->data;
    m = X_->dimensions[0];
    n = X_->dimensions[1];

    pdist_bray_curtis(X, dm, m, n);
  }
  return Py_BuildValue("d", 0.0);
}


static PyObject *pdist_mahalanobis_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *X_, *covinv_, *dm_;
  int m, n;
  double *dm;
  const double *X;
  const double *covinv;
  if (!PyArg_ParseTuple(args, "O!O!O!",
			&PyArray_Type, &X_,
			&PyArray_Type, &covinv_,
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    X = (const double*)X_->data;
    covinv = (const double*)covinv_->data;
    dm = (double*)dm_->data;
    m = X_->dimensions[0];
    n = X_->dimensions[1];

    pdist_mahalanobis(X, covinv, dm, m, n);
  }
  return Py_BuildValue("d", 0.0);
}


static PyObject *pdist_chebyshev_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *X_, *dm_;
  int m, n;
  double *dm;
  const double *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &X_,
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    X = (const double*)X_->data;
    dm = (double*)dm_->data;
    m = X_->dimensions[0];
    n = X_->dimensions[1];

    pdist_chebyshev(X, dm, m, n);
  }
  return Py_BuildValue("d", 0.0);
}


static PyObject *pdist_cosine_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *X_, *dm_, *norms_;
  int m, n;
  double *dm;
  const double *X, *norms;
  if (!PyArg_ParseTuple(args, "O!O!O!",
			&PyArray_Type, &X_,
			&PyArray_Type, &dm_,
			&PyArray_Type, &norms_)) {
    return 0;
  }
  else {
    X = (const double*)X_->data;
    dm = (double*)dm_->data;
    norms = (const double*)norms_->data;
    m = X_->dimensions[0];
    n = X_->dimensions[1];

    pdist_cosine(X, dm, m, n, norms);
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *pdist_seuclidean_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *X_, *dm_, *var_;
  int m, n;
  double *dm;
  const double *X, *var;
  if (!PyArg_ParseTuple(args, "O!O!O!",
			&PyArray_Type, &X_,
			&PyArray_Type, &var_,
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    X = (double*)X_->data;
    dm = (double*)dm_->data;
    var = (double*)var_->data;
    m = X_->dimensions[0];
    n = X_->dimensions[1];

    pdist_seuclidean(X, var, dm, m, n);
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *pdist_city_block_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *X_, *dm_;
  int m, n;
  double *dm;
  const double *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &X_,
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    X = (const double*)X_->data;
    dm = (double*)dm_->data;
    m = X_->dimensions[0];
    n = X_->dimensions[1];

    pdist_city_block(X, dm, m, n);
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *pdist_hamming_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *X_, *dm_;
  int m, n;
  double *dm;
  const double *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &X_,
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    X = (const double*)X_->data;
    dm = (double*)dm_->data;
    m = X_->dimensions[0];
    n = X_->dimensions[1];

    pdist_hamming(X, dm, m, n);
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *pdist_hamming_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *X_, *dm_;
  int m, n;
  double *dm;
  const char *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &X_,
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    X = (const char*)X_->data;
    dm = (double*)dm_->data;
    m = X_->dimensions[0];
    n = X_->dimensions[1];

    pdist_hamming_bool(X, dm, m, n);
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *pdist_jaccard_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *X_, *dm_;
  int m, n;
  double *dm;
  const double *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &X_,
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    X = (const double*)X_->data;
    dm = (double*)dm_->data;
    m = X_->dimensions[0];
    n = X_->dimensions[1];

    pdist_jaccard(X, dm, m, n);
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *pdist_jaccard_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *X_, *dm_;
  int m, n;
  double *dm;
  const char *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &X_,
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    X = (const char*)X_->data;
    dm = (double*)dm_->data;
    m = X_->dimensions[0];
    n = X_->dimensions[1];

    pdist_jaccard_bool(X, dm, m, n);
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *pdist_minkowski_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *X_, *dm_;
  int m, n;
  double *dm, *X;
  double p;
  if (!PyArg_ParseTuple(args, "O!O!d",
			&PyArray_Type, &X_,
			&PyArray_Type, &dm_,
			&p)) {
    return 0;
  }
  else {
    X = (double*)X_->data;
    dm = (double*)dm_->data;
    m = X_->dimensions[0];
    n = X_->dimensions[1];

    pdist_minkowski(X, dm, m, n, p);
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *pdist_weighted_minkowski_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *X_, *dm_, *w_;
  int m, n;
  double *dm, *X, *w;
  double p;
  if (!PyArg_ParseTuple(args, "O!O!dO!",
			&PyArray_Type, &X_,
			&PyArray_Type, &dm_,
			&p,
			&PyArray_Type, &w_)) {
    return 0;
  }
  else {
    X = (double*)X_->data;
    dm = (double*)dm_->data;
    w = (const double*)w_->data;
    m = X_->dimensions[0];
    n = X_->dimensions[1];

    pdist_weighted_minkowski(X, dm, m, n, p, w);
  }
  return Py_BuildValue("d", 0.0);
}


static PyObject *pdist_yule_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *X_, *dm_;
  int m, n;
  double *dm;
  const char *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &X_,
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    X = (const char*)X_->data;
    dm = (double*)dm_->data;
    m = X_->dimensions[0];
    n = X_->dimensions[1];

    pdist_yule_bool(X, dm, m, n);
  }
  return Py_BuildValue("");
}

static PyObject *pdist_matching_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *X_, *dm_;
  int m, n;
  double *dm;
  const char *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &X_,
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    X = (const char*)X_->data;
    dm = (double*)dm_->data;
    m = X_->dimensions[0];
    n = X_->dimensions[1];

    pdist_matching_bool(X, dm, m, n);
  }
  return Py_BuildValue("");
}

static PyObject *pdist_dice_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *X_, *dm_;
  int m, n;
  double *dm;
  const char *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &X_,
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    X = (const char*)X_->data;
    dm = (double*)dm_->data;
    m = X_->dimensions[0];
    n = X_->dimensions[1];

    pdist_dice_bool(X, dm, m, n);
  }
  return Py_BuildValue("");
}

static PyObject *pdist_rogerstanimoto_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *X_, *dm_;
  int m, n;
  double *dm;
  const char *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &X_,
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    X = (const char*)X_->data;
    dm = (double*)dm_->data;
    m = X_->dimensions[0];
    n = X_->dimensions[1];

    pdist_rogerstanimoto_bool(X, dm, m, n);
  }
  return Py_BuildValue("");
}

static PyObject *pdist_russellrao_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *X_, *dm_;
  int m, n;
  double *dm;
  const char *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &X_,
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    X = (const char*)X_->data;
    dm = (double*)dm_->data;
    m = X_->dimensions[0];
    n = X_->dimensions[1];

    pdist_russellrao_bool(X, dm, m, n);
  }
  return Py_BuildValue("");
}

static PyObject *pdist_kulsinski_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *X_, *dm_;
  int m, n;
  double *dm;
  const char *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &X_,
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    X = (const char*)X_->data;
    dm = (double*)dm_->data;
    m = X_->dimensions[0];
    n = X_->dimensions[1];

    pdist_kulsinski_bool(X, dm, m, n);
  }
  return Py_BuildValue("");
}

static PyObject *pdist_sokalmichener_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *X_, *dm_;
  int m, n;
  double *dm;
  const char *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &X_,
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    X = (const char*)X_->data;
    dm = (double*)dm_->data;
    m = X_->dimensions[0];
    n = X_->dimensions[1];

    pdist_sokalmichener_bool(X, dm, m, n);
  }
  return Py_BuildValue("");
}

static PyObject *pdist_sokalsneath_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *X_, *dm_;
  int m, n;
  double *dm;
  const char *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &X_,
			&PyArray_Type, &dm_)) {
    return 0;
  }
  else {
    X = (const char*)X_->data;
    dm = (double*)dm_->data;
    m = X_->dimensions[0];
    n = X_->dimensions[1];

    pdist_sokalsneath_bool(X, dm, m, n);
  }
  return Py_BuildValue("");
}

static PyObject *to_squareform_from_vector_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *M_, *v_;
  int n;
  const double *v;
  double *M;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &M_,
			&PyArray_Type, &v_)) {
    return 0;
  }
  else {
    M = (double*)M_->data;
    v = (const double*)v_->data;
    n = M_->dimensions[0];
    dist_to_squareform_from_vector(M, v, n);
  }
  return Py_BuildValue("d", 0.0);
}

static PyObject *to_vector_from_squareform_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *M_, *v_;
  int n;
  double *v;
  const double *M;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &M_,
			&PyArray_Type, &v_)) {
    return 0;
  }
  else {
    M = (const double*)M_->data;
    v = (double*)v_->data;
    n = M_->dimensions[0];
    dist_to_vector_from_squareform(M, v, n);
  }
  return Py_BuildValue("d", 0.0);
}


static PyMethodDef _distanceWrapMethods[] = {
  {"cdist_bray_curtis_wrap", cdist_bray_curtis_wrap, METH_VARARGS},
  {"cdist_canberra_wrap", cdist_canberra_wrap, METH_VARARGS},
  {"cdist_chebyshev_wrap", cdist_chebyshev_wrap, METH_VARARGS},
  {"cdist_city_block_wrap", cdist_city_block_wrap, METH_VARARGS},
  {"cdist_cosine_wrap", cdist_cosine_wrap, METH_VARARGS},
  {"cdist_dice_bool_wrap", cdist_dice_bool_wrap, METH_VARARGS},
  {"cdist_euclidean_wrap", cdist_euclidean_wrap, METH_VARARGS},
  {"cdist_sqeuclidean_wrap", cdist_sqeuclidean_wrap, METH_VARARGS},
  {"cdist_hamming_wrap", cdist_hamming_wrap, METH_VARARGS},
  {"cdist_hamming_bool_wrap", cdist_hamming_bool_wrap, METH_VARARGS},
  {"cdist_jaccard_wrap", cdist_jaccard_wrap, METH_VARARGS},
  {"cdist_jaccard_bool_wrap", cdist_jaccard_bool_wrap, METH_VARARGS},
  {"cdist_kulsinski_bool_wrap", cdist_kulsinski_bool_wrap, METH_VARARGS},
  {"cdist_mahalanobis_wrap", cdist_mahalanobis_wrap, METH_VARARGS},
  {"cdist_matching_bool_wrap", cdist_matching_bool_wrap, METH_VARARGS},
  {"cdist_minkowski_wrap", cdist_minkowski_wrap, METH_VARARGS},
  {"cdist_weighted_minkowski_wrap", cdist_weighted_minkowski_wrap, METH_VARARGS},
  {"cdist_rogerstanimoto_bool_wrap", cdist_rogerstanimoto_bool_wrap, METH_VARARGS},
  {"cdist_russellrao_bool_wrap", cdist_russellrao_bool_wrap, METH_VARARGS},
  {"cdist_seuclidean_wrap", cdist_seuclidean_wrap, METH_VARARGS},
  {"cdist_sokalmichener_bool_wrap", cdist_sokalmichener_bool_wrap, METH_VARARGS},
  {"cdist_sokalsneath_bool_wrap", cdist_sokalsneath_bool_wrap, METH_VARARGS},
  {"cdist_yule_bool_wrap", cdist_yule_bool_wrap, METH_VARARGS},
  {"pdist_bray_curtis_wrap", pdist_bray_curtis_wrap, METH_VARARGS},
  {"pdist_canberra_wrap", pdist_canberra_wrap, METH_VARARGS},
  {"pdist_chebyshev_wrap", pdist_chebyshev_wrap, METH_VARARGS},
  {"pdist_city_block_wrap", pdist_city_block_wrap, METH_VARARGS},
  {"pdist_cosine_wrap", pdist_cosine_wrap, METH_VARARGS},
  {"pdist_dice_bool_wrap", pdist_dice_bool_wrap, METH_VARARGS},
  {"pdist_euclidean_wrap", pdist_euclidean_wrap, METH_VARARGS},
  {"pdist_sqeuclidean_wrap", pdist_sqeuclidean_wrap, METH_VARARGS},
  {"pdist_hamming_wrap", pdist_hamming_wrap, METH_VARARGS},
  {"pdist_hamming_bool_wrap", pdist_hamming_bool_wrap, METH_VARARGS},
  {"pdist_jaccard_wrap", pdist_jaccard_wrap, METH_VARARGS},
  {"pdist_jaccard_bool_wrap", pdist_jaccard_bool_wrap, METH_VARARGS},
  {"pdist_kulsinski_bool_wrap", pdist_kulsinski_bool_wrap, METH_VARARGS},
  {"pdist_mahalanobis_wrap", pdist_mahalanobis_wrap, METH_VARARGS},
  {"pdist_matching_bool_wrap", pdist_matching_bool_wrap, METH_VARARGS},
  {"pdist_minkowski_wrap", pdist_minkowski_wrap, METH_VARARGS},
  {"pdist_weighted_minkowski_wrap", pdist_weighted_minkowski_wrap, METH_VARARGS},
  {"pdist_rogerstanimoto_bool_wrap", pdist_rogerstanimoto_bool_wrap, METH_VARARGS},
  {"pdist_russellrao_bool_wrap", pdist_russellrao_bool_wrap, METH_VARARGS},
  {"pdist_seuclidean_wrap", pdist_seuclidean_wrap, METH_VARARGS},
  {"pdist_sokalmichener_bool_wrap", pdist_sokalmichener_bool_wrap, METH_VARARGS},
  {"pdist_sokalsneath_bool_wrap", pdist_sokalsneath_bool_wrap, METH_VARARGS},
  {"pdist_yule_bool_wrap", pdist_yule_bool_wrap, METH_VARARGS},
  {"to_squareform_from_vector_wrap",
   to_squareform_from_vector_wrap, METH_VARARGS},
  {"to_vector_from_squareform_wrap",
   to_vector_from_squareform_wrap, METH_VARARGS},
  {NULL, NULL}     /* Sentinel - marks the end of this structure */
};

#if PY_VERSION_HEX >= 0x03000000
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

PyObject *PyInit__distance_wrap(void)
{
    PyObject *m;

    m = PyModule_Create(&moduledef);
    import_array();

    return m;
}
#else
PyMODINIT_FUNC init_distance_wrap(void)
{
  (void) Py_InitModule("_distance_wrap", _distanceWrapMethods);
  import_array();  // Must be present for NumPy.  Called first after above line.
}
#endif
