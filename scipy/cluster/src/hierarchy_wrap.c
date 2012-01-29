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
    return NULL;
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
    if (linkage((double*)dm->data, (double*)Z->data,
                0, 0, n, 0, 0, df, method) == -1) {
      PyErr_SetString(PyExc_MemoryError,
                      "out of memory while computing linkage");
      return NULL;
    }
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
    return NULL;
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
    if (linkage((double*)dm->data, (double*)Z->data, (double*)X->data,
	              m, n, 1, 1, df, method) == -1) {
      PyErr_SetString(PyExc_MemoryError,
                      "out of memory while computing linkage");
      return NULL;
    }
  }
  return Py_BuildValue("d", 0.0);
}

extern PyObject *calculate_cluster_sizes_wrap(PyObject *self, PyObject *args) {
  int n;
  PyArrayObject *Z, *cs_;
  if (!PyArg_ParseTuple(args, "O!O!i",
			&PyArray_Type, &Z,
			&PyArray_Type, &cs_,
			&n)) {
    return 0;
  }
  calculate_cluster_sizes((const double*)Z->data, (double*)cs_->data, n);
  return Py_BuildValue("");
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


extern PyObject *leaders_wrap(PyObject *self, PyObject *args) {
  PyArrayObject *Z_, *T_, *L_, *M_;
  int kk, n, res;
  if (!PyArg_ParseTuple(args, "O!O!O!O!ii",
			&PyArray_Type, &Z_,
			&PyArray_Type, &T_,
			&PyArray_Type, &L_,
			&PyArray_Type, &M_,
			&kk, &n)) {
    return 0;
  }
  else {
    res = leaders((double*)Z_->data, (int*)T_->data,
		  (int*)L_->data, (int*)M_->data, kk, n);
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
  {"get_max_dist_for_each_cluster_wrap",
   get_max_dist_for_each_cluster_wrap, METH_VARARGS},
  {"get_max_Rfield_for_each_cluster_wrap",
   get_max_Rfield_for_each_cluster_wrap, METH_VARARGS},
  {"inconsistent_wrap", inconsistent_wrap, METH_VARARGS},
  {"leaders_wrap", leaders_wrap, METH_VARARGS},
  {"linkage_euclid_wrap", linkage_euclid_wrap, METH_VARARGS},
  {"linkage_wrap", linkage_wrap, METH_VARARGS},
  {"prelist_wrap", prelist_wrap, METH_VARARGS},
  {NULL, NULL}     /* Sentinel - marks the end of this structure */
};

#if defined(SCIPY_PY3K)
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_vq",
    NULL,
    -1,
    _hierarchyWrapMethods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyObject *PyInit__hierarchy_wrap(void)
{
    PyObject *m;

    m = PyModule_Create(&moduledef);
    import_array();

    return m;
}
#else
PyMODINIT_FUNC init_hierarchy_wrap(void)  {
  (void) Py_InitModule("_hierarchy_wrap", _hierarchyWrapMethods);
  import_array();  // Must be present for NumPy.  Called first after above line.
}
#endif
