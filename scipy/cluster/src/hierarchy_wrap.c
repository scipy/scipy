/**
 * hierarchy_wrap.c
 *
 * Author: Damian Eads
 * Date:   September 22, 2007
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

#include "hierarchy.h"
#include "Python.h"
#include <numpy/arrayobject.h>
#include <stdio.h>

extern PyObject *linkage_wrap(PyObject *self, PyObject *args) {
  int method, n;
  PyArrayObject *dm, *Z;
  distfunc *df;
  if (!PyArg_ParseTuple(args, "O!O!ii",
			&PyArray_Type, &dm,
			&PyArray_Type, &Z,
			&n,
			&method)) {
    return 0;
  }
  else {
    switch (method) {
    case CPY_LINKAGE_SINGLE:
      df = dist_single;
      break;
    case CPY_LINKAGE_COMPLETE:
      df = dist_complete;
      break;
    case CPY_LINKAGE_AVERAGE:
      df = dist_average;
      break;
    case CPY_LINKAGE_WEIGHTED:
      df = dist_weighted;
      break;
    default:
      /** Report an error. */
      df = 0;
      break;
    }
    linkage((double*)dm->data, (double*)Z->data, 0, 0, n, 0, 0, df, method);
  }
  return Py_BuildValue("d", 0.0);
}

extern PyObject *linkage_euclid_wrap(PyObject *self, PyObject *args) {
  int method, m, n, ml;
  PyArrayObject *dm, *Z, *X;
  distfunc *df;
  if (!PyArg_ParseTuple(args, "O!O!O!iii",
			&PyArray_Type, &dm,
			&PyArray_Type, &Z,
			&PyArray_Type, &X,
			&m,
			&n,
			&method)) {
    return 0;
  }
  else {
    ml = 0;
    /**    fprintf(stderr, "m: %d, n: %d\n", m, n);**/
    switch (method) {
    case CPY_LINKAGE_CENTROID:
      df = dist_centroid;
      break;
    case CPY_LINKAGE_MEDIAN:
      df = dist_centroid;
      break;
    case CPY_LINKAGE_WARD:
      df = dist_ward;
      //      ml = 1;
      break;
    default:
      /** Report an error. */
      df = 0;
      break;
    }
    linkage((double*)dm->data, (double*)Z->data, (double*)X->data,
	    m, n, 1, 1, df, method);
  }
  return Py_BuildValue("d", 0.0);
}

extern PyObject *calculate_cluster_sizes_wrap(PyObject *self, PyObject *args) {
  int n;
  PyArrayObject *Z, *CS;
  if (!PyArg_ParseTuple(args, "O!O!i",
			&PyArray_Type, &Z,
			&PyArray_Type, &CS,
			&n)) {
    return 0;
  }
  calculate_cluster_sizes((const double*)Z->data, (double*)CS->data, n);
  return Py_BuildValue("d", 0.0);
}

extern PyObject *get_max_dist_for_each_cluster_wrap(PyObject *self,
						    PyObject *args) {
  int n;
  PyArrayObject *Z, *md;
  if (!PyArg_ParseTuple(args, "O!O!i",
			&PyArray_Type, &Z,
			&PyArray_Type, &md,
			&n)) {
    return 0;
  }
  get_max_dist_for_each_cluster((const double*)Z->data, (double*)md->data, n);
  return Py_BuildValue("");
}

extern PyObject *get_max_Rfield_for_each_cluster_wrap(PyObject *self,
						      PyObject *args) {
  int n, rf;
  PyArrayObject *Z, *R, *max_rfs;
  if (!PyArg_ParseTuple(args, "O!O!O!ii",
			&PyArray_Type, &Z,
			&PyArray_Type, &R,
			&PyArray_Type, &max_rfs,
			&n, &rf)) {
    return 0;
  }
  get_max_Rfield_for_each_cluster((const double *)Z->data,
				  (const double *)R->data,
				  (double *)max_rfs->data, n, rf);
  return Py_BuildValue("");
}

extern PyObject *prelist_wrap(PyObject *self, PyObject *args) {
  int n;
  PyArrayObject *Z, *ML;
  if (!PyArg_ParseTuple(args, "O!O!i",
			&PyArray_Type, &Z,
			&PyArray_Type, &ML,
			&n)) {
    return 0;
  }
  form_member_list((const double *)Z->data, (int *)ML->data, n);
  return Py_BuildValue("d", 0.0);
}

extern PyObject *cluster_in_wrap(PyObject *self, PyObject *args) {
  int n;
  double cutoff;
  PyArrayObject *Z, *R, *T;
  if (!PyArg_ParseTuple(args, "O!O!O!di",
			&PyArray_Type, &Z,
			&PyArray_Type, &R,
			&PyArray_Type, &T,
			&cutoff,
			&n)) {
    return 0;
  }
  form_flat_clusters_from_in((const double *)Z->data, (const double *)R->data,
			     (int *)T->data, cutoff, n);

  return Py_BuildValue("d", 0.0);
}

extern PyObject *cluster_dist_wrap(PyObject *self, PyObject *args) {
  int n;
  double cutoff;
  PyArrayObject *Z, *T;
  if (!PyArg_ParseTuple(args, "O!O!di",
			&PyArray_Type, &Z,
			&PyArray_Type, &T,
			&cutoff,
			&n)) {
    return 0;
  }
  form_flat_clusters_from_dist((const double *)Z->data,
			       (int *)T->data, cutoff, n);

  return Py_BuildValue("d", 0.0);
}

extern PyObject *cluster_monocrit_wrap(PyObject *self, PyObject *args) {
  int n;
  double cutoff;
  PyArrayObject *Z, *MV, *T;
  if (!PyArg_ParseTuple(args, "O!O!O!di",
			&PyArray_Type, &Z,
			&PyArray_Type, &MV,
			&PyArray_Type, &T,
			&cutoff,
			&n)) {
    return 0;
  }
  form_flat_clusters_from_monotonic_criterion((const double *)Z->data,
					      (const double *)MV->data,
					      (int *)T->data,
					      cutoff,
					      n);

  form_flat_clusters_from_dist((const double *)Z->data,
			       (int *)T->data, cutoff, n);

  return Py_BuildValue("d", 0.0);
}



extern PyObject *cluster_maxclust_dist_wrap(PyObject *self, PyObject *args) {
  int n, mc;
  PyArrayObject *Z, *T;
  if (!PyArg_ParseTuple(args, "O!O!ii",
			&PyArray_Type, &Z,
			&PyArray_Type, &T,
			&n, &mc)) {
    return 0;
  }
  form_flat_clusters_maxclust_dist((const double*)Z->data, (int *)T->data,
				   n, mc);

  return Py_BuildValue("");
}


extern PyObject *cluster_maxclust_monocrit_wrap(PyObject *self, PyObject *args) {
  int n, mc;
  PyArrayObject *Z, *MC, *T;
  if (!PyArg_ParseTuple(args, "O!O!O!ii",
			&PyArray_Type, &Z,
			&PyArray_Type, &MC,
			&PyArray_Type, &T,
			&n, &mc)) {
    return 0;
  }
  form_flat_clusters_maxclust_monocrit((const double *)Z->data,
				       (const double *)MC->data,
				       (int *)T->data, n, mc);

  return Py_BuildValue("");
}


extern PyObject *inconsistent_wrap(PyObject *self, PyObject *args) {
  int n, d;
  PyArrayObject *Z, *R;
  if (!PyArg_ParseTuple(args, "O!O!ii",
			&PyArray_Type, &Z,
			&PyArray_Type, &R,
			&n, &d)) {
    return 0;
  }
  inconsistency_calculation_alt((const double*)Z->data, (double*)R->data, n, d);
  return Py_BuildValue("d", 0.0);
}

extern PyObject *cophenetic_distances_wrap(PyObject *self, PyObject *args) {
  int n;
  PyArrayObject *Z, *d;
  if (!PyArg_ParseTuple(args, "O!O!i",
			&PyArray_Type, &Z,
			&PyArray_Type, &d,
			&n)) {
    return 0;
  }
  cophenetic_distances((const double*)Z->data, (double*)d->data, n);
  return Py_BuildValue("d", 0.0);
}

extern PyObject *chopmin_ns_ij_wrap(PyObject *self, PyObject *args) {
  int mini, minj, n;
  PyArrayObject *row;
  if (!PyArg_ParseTuple(args, "O!iii",
			&PyArray_Type, &row,
			&mini,
			&minj,
			&n)) {
    return 0;
  }
  chopmins_ns_ij((double*)row->data, mini, minj, n);
  return Py_BuildValue("d", 0.0);
}


extern PyObject *chopmin_ns_i_wrap(PyObject *self, PyObject *args) {
  int mini, n;
  PyArrayObject *row;
  if (!PyArg_ParseTuple(args, "O!ii",
			&PyArray_Type, &row,
			&mini,
			&n)) {
    return 0;
  }
  chopmins_ns_i((double*)row->data, mini, n);
  return Py_BuildValue("d", 0.0);
}

extern PyObject *chopmins_wrap(PyObject *self, PyObject *args) {
  int mini, minj, n;
  PyArrayObject *row;
  if (!PyArg_ParseTuple(args, "O!iii",
			&PyArray_Type, &row,
			&mini,
			&minj,
			&n)) {
    return 0;
  }
  chopmins((int*)row->data, mini, minj, n);
  return Py_BuildValue("d", 0.0);
}

extern PyObject *dot_product_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *_d1, *_d2;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &_d1,
			&PyArray_Type, &_d2)) {
    return 0;
  }
  return Py_BuildValue("d", dot_product((const double*)_d1->data,
					(const double*)_d2->data,
					_d1->dimensions[0]));
}

extern PyObject *to_squareform_from_vector_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *_M, *_v;
  int n;
  const double *v;
  double *M;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &_M,
			&PyArray_Type, &_v)) {
    return 0;
  }
  else {
    M = (double*)_M->data;
    v = (const double*)_v->data;
    n = _M->dimensions[0];
    dist_to_squareform_from_vector(M, v, n);
  }
  return Py_BuildValue("d", 0.0);
}

extern PyObject *to_vector_from_squareform_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *_M, *_v;
  int n;
  double *v;
  const double *M;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &_M,
			&PyArray_Type, &_v)) {
    return 0;
  }
  else {
    M = (const double*)_M->data;
    v = (double*)_v->data;
    n = _M->dimensions[0];
    dist_to_vector_from_squareform(M, v, n);
  }
  return Py_BuildValue("d", 0.0);
}

extern PyObject *pdist_euclidean_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *_X, *_dm;
  int m, n;
  double *dm;
  const double *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &_X,
			&PyArray_Type, &_dm)) {
    return 0;
  }
  else {
    X = (const double*)_X->data;
    dm = (double*)_dm->data;
    m = _X->dimensions[0];
    n = _X->dimensions[1];
    
    pdist_euclidean(X, dm, m, n);
  }
  return Py_BuildValue("d", 0.0);
}

extern PyObject *pdist_canberra_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *_X, *_dm;
  int m, n;
  double *dm;
  const double *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &_X,
			&PyArray_Type, &_dm)) {
    return 0;
  }
  else {
    X = (const double*)_X->data;
    dm = (double*)_dm->data;
    m = _X->dimensions[0];
    n = _X->dimensions[1];
    
    pdist_canberra(X, dm, m, n);
  }
  return Py_BuildValue("d", 0.0);
}

extern PyObject *pdist_bray_curtis_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *_X, *_dm;
  int m, n;
  double *dm;
  const double *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &_X,
			&PyArray_Type, &_dm)) {
    return 0;
  }
  else {
    X = (const double*)_X->data;
    dm = (double*)_dm->data;
    m = _X->dimensions[0];
    n = _X->dimensions[1];
    
    pdist_bray_curtis(X, dm, m, n);
  }
  return Py_BuildValue("d", 0.0);
}


extern PyObject *pdist_mahalanobis_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *_X, *_covinv, *_dm;
  int m, n;
  double *dm;
  const double *X;
  const double *covinv;
  if (!PyArg_ParseTuple(args, "O!O!O!",
			&PyArray_Type, &_X,
			&PyArray_Type, &_covinv,
			&PyArray_Type, &_dm)) {
    return 0;
  }
  else {
    X = (const double*)_X->data;
    covinv = (const double*)_covinv->data;
    dm = (double*)_dm->data;
    m = _X->dimensions[0];
    n = _X->dimensions[1];
    
    pdist_mahalanobis(X, covinv, dm, m, n);
  }
  return Py_BuildValue("d", 0.0);
}


extern PyObject *pdist_chebyshev_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *_X, *_dm;
  int m, n;
  double *dm;
  const double *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &_X,
			&PyArray_Type, &_dm)) {
    return 0;
  }
  else {
    X = (const double*)_X->data;
    dm = (double*)_dm->data;
    m = _X->dimensions[0];
    n = _X->dimensions[1];
    
    pdist_chebyshev(X, dm, m, n);
  }
  return Py_BuildValue("d", 0.0);
}


extern PyObject *pdist_cosine_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *_X, *_dm, *_norms;
  int m, n;
  double *dm;
  const double *X, *norms;
  if (!PyArg_ParseTuple(args, "O!O!O!",
			&PyArray_Type, &_X,
			&PyArray_Type, &_dm,
			&PyArray_Type, &_norms)) {
    return 0;
  }
  else {
    X = (const double*)_X->data;
    dm = (double*)_dm->data;
    norms = (const double*)_norms->data;
    m = _X->dimensions[0];
    n = _X->dimensions[1];
    
    pdist_cosine(X, dm, m, n, norms);
  }
  return Py_BuildValue("d", 0.0);
}

extern PyObject *pdist_seuclidean_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *_X, *_dm, *_var;
  int m, n;
  double *dm;
  const double *X, *var;
  if (!PyArg_ParseTuple(args, "O!O!O!",
			&PyArray_Type, &_X,
			&PyArray_Type, &_var,
			&PyArray_Type, &_dm)) {
    return 0;
  }
  else {
    X = (double*)_X->data;
    dm = (double*)_dm->data;
    var = (double*)_var->data;
    m = _X->dimensions[0];
    n = _X->dimensions[1];
    
    pdist_seuclidean(X, var, dm, m, n);
  }
  return Py_BuildValue("d", 0.0);
}

extern PyObject *pdist_city_block_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *_X, *_dm;
  int m, n;
  double *dm;
  const double *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &_X,
			&PyArray_Type, &_dm)) {
    return 0;
  }
  else {
    X = (const double*)_X->data;
    dm = (double*)_dm->data;
    m = _X->dimensions[0];
    n = _X->dimensions[1];
    
    pdist_city_block(X, dm, m, n);
  }
  return Py_BuildValue("d", 0.0);
}

extern PyObject *pdist_hamming_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *_X, *_dm;
  int m, n;
  double *dm;
  const double *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &_X,
			&PyArray_Type, &_dm)) {
    return 0;
  }
  else {
    X = (const double*)_X->data;
    dm = (double*)_dm->data;
    m = _X->dimensions[0];
    n = _X->dimensions[1];
    
    pdist_hamming(X, dm, m, n);
  }
  return Py_BuildValue("d", 0.0);
}

extern PyObject *pdist_hamming_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *_X, *_dm;
  int m, n;
  double *dm;
  const char *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &_X,
			&PyArray_Type, &_dm)) {
    return 0;
  }
  else {
    X = (const char*)_X->data;
    dm = (double*)_dm->data;
    m = _X->dimensions[0];
    n = _X->dimensions[1];
    
    pdist_hamming_bool(X, dm, m, n);
  }
  return Py_BuildValue("d", 0.0);
}

extern PyObject *pdist_jaccard_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *_X, *_dm;
  int m, n;
  double *dm;
  const double *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &_X,
			&PyArray_Type, &_dm)) {
    return 0;
  }
  else {
    X = (const double*)_X->data;
    dm = (double*)_dm->data;
    m = _X->dimensions[0];
    n = _X->dimensions[1];
    
    pdist_jaccard(X, dm, m, n);
  }
  return Py_BuildValue("d", 0.0);
}

extern PyObject *pdist_jaccard_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *_X, *_dm;
  int m, n;
  double *dm;
  const char *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &_X,
			&PyArray_Type, &_dm)) {
    return 0;
  }
  else {
    X = (const char*)_X->data;
    dm = (double*)_dm->data;
    m = _X->dimensions[0];
    n = _X->dimensions[1];
    
    pdist_jaccard_bool(X, dm, m, n);
  }
  return Py_BuildValue("d", 0.0);
}

extern PyObject *pdist_minkowski_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *_X, *_dm;
  int m, n;
  double *dm, *X;
  double p;
  if (!PyArg_ParseTuple(args, "O!O!d",
			&PyArray_Type, &_X,
			&PyArray_Type, &_dm,
			&p)) {
    return 0;
  }
  else {
    X = (double*)_X->data;
    dm = (double*)_dm->data;
    m = _X->dimensions[0];
    n = _X->dimensions[1];
    
    pdist_minkowski(X, dm, m, n, p);
  }
  return Py_BuildValue("d", 0.0);
}


extern PyObject *pdist_yule_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *_X, *_dm;
  int m, n;
  double *dm;
  const char *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &_X,
			&PyArray_Type, &_dm)) {
    return 0;
  }
  else {
    X = (const char*)_X->data;
    dm = (double*)_dm->data;
    m = _X->dimensions[0];
    n = _X->dimensions[1];
    
    pdist_yule_bool(X, dm, m, n);
  }
  return Py_BuildValue("");
}

extern PyObject *pdist_matching_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *_X, *_dm;
  int m, n;
  double *dm;
  const char *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &_X,
			&PyArray_Type, &_dm)) {
    return 0;
  }
  else {
    X = (const char*)_X->data;
    dm = (double*)_dm->data;
    m = _X->dimensions[0];
    n = _X->dimensions[1];
    
    pdist_matching_bool(X, dm, m, n);
  }
  return Py_BuildValue("");
}

extern PyObject *pdist_dice_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *_X, *_dm;
  int m, n;
  double *dm;
  const char *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &_X,
			&PyArray_Type, &_dm)) {
    return 0;
  }
  else {
    X = (const char*)_X->data;
    dm = (double*)_dm->data;
    m = _X->dimensions[0];
    n = _X->dimensions[1];
    
    pdist_dice_bool(X, dm, m, n);
  }
  return Py_BuildValue("");
}

extern PyObject *pdist_rogerstanimoto_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *_X, *_dm;
  int m, n;
  double *dm;
  const char *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &_X,
			&PyArray_Type, &_dm)) {
    return 0;
  }
  else {
    X = (const char*)_X->data;
    dm = (double*)_dm->data;
    m = _X->dimensions[0];
    n = _X->dimensions[1];
    
    pdist_rogerstanimoto_bool(X, dm, m, n);
  }
  return Py_BuildValue("");
}

extern PyObject *pdist_russellrao_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *_X, *_dm;
  int m, n;
  double *dm;
  const char *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &_X,
			&PyArray_Type, &_dm)) {
    return 0;
  }
  else {
    X = (const char*)_X->data;
    dm = (double*)_dm->data;
    m = _X->dimensions[0];
    n = _X->dimensions[1];
    
    pdist_russellrao_bool(X, dm, m, n);
  }
  return Py_BuildValue("");
}

extern PyObject *pdist_kulsinski_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *_X, *_dm;
  int m, n;
  double *dm;
  const char *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &_X,
			&PyArray_Type, &_dm)) {
    return 0;
  }
  else {
    X = (const char*)_X->data;
    dm = (double*)_dm->data;
    m = _X->dimensions[0];
    n = _X->dimensions[1];
    
    pdist_kulsinski_bool(X, dm, m, n);
  }
  return Py_BuildValue("");
}

extern PyObject *pdist_sokalmichener_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *_X, *_dm;
  int m, n;
  double *dm;
  const char *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &_X,
			&PyArray_Type, &_dm)) {
    return 0;
  }
  else {
    X = (const char*)_X->data;
    dm = (double*)_dm->data;
    m = _X->dimensions[0];
    n = _X->dimensions[1];
    
    pdist_sokalmichener_bool(X, dm, m, n);
  }
  return Py_BuildValue("");
}

extern PyObject *pdist_sokalsneath_bool_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *_X, *_dm;
  int m, n;
  double *dm;
  const char *X;
  if (!PyArg_ParseTuple(args, "O!O!",
			&PyArray_Type, &_X,
			&PyArray_Type, &_dm)) {
    return 0;
  }
  else {
    X = (const char*)_X->data;
    dm = (double*)_dm->data;
    m = _X->dimensions[0];
    n = _X->dimensions[1];
    
    pdist_sokalsneath_bool(X, dm, m, n);
  }
  return Py_BuildValue("");
}

extern PyObject *leaders_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *_Z, *_T, *_L, *_M;
  int kk, n, res;
  if (!PyArg_ParseTuple(args, "O!O!O!O!ii",
			&PyArray_Type, &_Z,
			&PyArray_Type, &_T,
			&PyArray_Type, &_L,
			&PyArray_Type, &_M,
			&kk, &n)) {
    return 0;
  }
  else {
    res = leaders((double*)_Z->data, (int*)_T->data,
		  (int*)_L->data, (int*)_M->data, kk, n);
  }
  return Py_BuildValue("i", res);
}

static PyMethodDef _hierarchyWrapMethods[] = {
  {"calculate_cluster_sizes_wrap", calculate_cluster_sizes_wrap, METH_VARARGS},
  {"chopmins", chopmins_wrap, METH_VARARGS},
  {"chopmins_ns_i", chopmin_ns_i_wrap, METH_VARARGS},
  {"chopmins_ns_ij", chopmin_ns_ij_wrap, METH_VARARGS},
  {"cluster_in_wrap", cluster_in_wrap, METH_VARARGS},
  {"cluster_dist_wrap", cluster_dist_wrap, METH_VARARGS},
  {"cluster_maxclust_dist_wrap", cluster_maxclust_dist_wrap, METH_VARARGS},
  {"cluster_maxclust_monocrit_wrap", cluster_maxclust_monocrit_wrap, METH_VARARGS},
  {"cluster_monocrit_wrap", cluster_monocrit_wrap, METH_VARARGS},
  {"cophenetic_distances_wrap", cophenetic_distances_wrap, METH_VARARGS},
  {"dot_product_wrap", dot_product_wrap, METH_VARARGS},
  {"get_max_dist_for_each_cluster_wrap",
   get_max_dist_for_each_cluster_wrap, METH_VARARGS},
  {"get_max_Rfield_for_each_cluster_wrap",
   get_max_Rfield_for_each_cluster_wrap, METH_VARARGS},
  {"inconsistent_wrap", inconsistent_wrap, METH_VARARGS},
  {"leaders_wrap", leaders_wrap, METH_VARARGS},
  {"linkage_euclid_wrap", linkage_euclid_wrap, METH_VARARGS},
  {"linkage_wrap", linkage_wrap, METH_VARARGS},
  {"pdist_bray_curtis_wrap", pdist_bray_curtis_wrap, METH_VARARGS},
  {"pdist_canberra_wrap", pdist_canberra_wrap, METH_VARARGS},
  {"pdist_chebyshev_wrap", pdist_chebyshev_wrap, METH_VARARGS},
  {"pdist_city_block_wrap", pdist_city_block_wrap, METH_VARARGS},
  {"pdist_cosine_wrap", pdist_cosine_wrap, METH_VARARGS},
  {"pdist_dice_bool_wrap", pdist_dice_bool_wrap, METH_VARARGS},
  {"pdist_euclidean_wrap", pdist_euclidean_wrap, METH_VARARGS},
  {"pdist_hamming_wrap", pdist_hamming_wrap, METH_VARARGS},
  {"pdist_hamming_bool_wrap", pdist_hamming_bool_wrap, METH_VARARGS},
  {"pdist_jaccard_wrap", pdist_jaccard_wrap, METH_VARARGS},
  {"pdist_jaccard_bool_wrap", pdist_jaccard_bool_wrap, METH_VARARGS},
  {"pdist_kulsinski_bool_wrap", pdist_kulsinski_bool_wrap, METH_VARARGS},
  {"pdist_mahalanobis_wrap", pdist_mahalanobis_wrap, METH_VARARGS},
  {"pdist_matching_bool_wrap", pdist_matching_bool_wrap, METH_VARARGS},
  {"pdist_minkowski_wrap", pdist_minkowski_wrap, METH_VARARGS},
  {"pdist_rogerstanimoto_bool_wrap", pdist_rogerstanimoto_bool_wrap, METH_VARARGS},
  {"pdist_russellrao_bool_wrap", pdist_russellrao_bool_wrap, METH_VARARGS},
  {"pdist_seuclidean_wrap", pdist_seuclidean_wrap, METH_VARARGS},
  {"pdist_sokalmichener_bool_wrap", pdist_sokalmichener_bool_wrap, METH_VARARGS},
  {"pdist_sokalsneath_bool_wrap", pdist_sokalsneath_bool_wrap, METH_VARARGS},
  {"pdist_yule_bool_wrap", pdist_yule_bool_wrap, METH_VARARGS},
  {"prelist_wrap", prelist_wrap, METH_VARARGS},
  {"to_squareform_from_vector_wrap",
   to_squareform_from_vector_wrap, METH_VARARGS},
  {"to_vector_from_squareform_wrap",
   to_vector_from_squareform_wrap, METH_VARARGS},
  {NULL, NULL}     /* Sentinel - marks the end of this structure */
};

void init_hierarchy_wrap(void)  {
  (void) Py_InitModule("_hierarchy_wrap", _hierarchyWrapMethods);
  import_array();  // Must be present for NumPy.  Called first after above line.
}
