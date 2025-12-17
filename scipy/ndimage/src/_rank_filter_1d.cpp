// Copyright (c) 2011 ashelly.myopenid.com under
// <http://www.opensource.org/licenses/mit-license>
// Modified in 2024 by Gideon Genadi Kogan

#include "Python.h"
#include "numpy/arrayobject.h"

#include <memory>
#include <new>

//
// This is used for rank keeping.
//
class Mediator
{
public:
  int *mem;  // Single array of ints to hold the memory for `pos` and `heap`.
  int *pos;  // index into `heap` for each value
  int *heap; // max/rank/min heap holding indexes into `data`.
  int N;     // allocated size.
  int idx;   // position in circular queue
  int minCt; // count of items in min heap
  int maxCt; // count of items in max heap

  //
  // This constructor will throw an exception if memory allocation
  // fails.  The caller *must* handle this exception so that it
  // does not propagate out to Python-land.
  //
  Mediator(int nItems, int rank) {
    mem = new int[2 * nItems];  // Might throw std::bad_alloc
    pos = mem;                  // `pos` uses the first `nItems` elements of `mem`.
    heap = mem + nItems + rank; // `heap` uses the second `nItems`; it actually
                                // points to `rank` elements into the second
                                // block of `nItems` elements.
    N = nItems;
    idx = 0;
    minCt = nItems - rank - 1;
    maxCt = rank;
    while (nItems--) {
      pos[nItems] = nItems - rank;
      heap[pos[nItems]] = nItems;
    }
  }

  ~Mediator() {
    delete [] mem;
  }
};

typedef enum {
  NEAREST = 0,
  WRAP = 1,
  REFLECT = 2,
  MIRROR = 3,
  CONSTANT = 4,
} Mode;

/*--- Helper Functions ---*/

// returns 1 if heap[i] < heap[j]
template <typename T> inline int mmless(T *data, Mediator *m, int i, int j) {
  return (data[m->heap[i]] < data[m->heap[j]]);
}

// swaps items i&j in heap, maintains indexes
int mmexchange(Mediator *m, int i, int j) {
  int t = m->heap[i];
  m->heap[i] = m->heap[j];
  m->heap[j] = t;
  m->pos[m->heap[i]] = i;
  m->pos[m->heap[j]] = j;
  return 1;
}

// swaps items i & j if i < j;  returns true if swapped
template <typename T> inline int mmCmpExch(T *data, Mediator *m, int i, int j) {
  return (mmless(data, m, i, j) && mmexchange(m, i, j));
}

// maintains minheap property for all items below i.
template <typename T> void minSortDown(T *data, Mediator *m, int i) {
  for (i *= 2; i <= m->minCt; i *= 2) {
    if (i < m->minCt && mmless(data, m, i + 1, i)) {
      ++i;
    }
    if (!mmCmpExch(data, m, i, i / 2)) {
      break;
    }
  }
}

// maintains maxheap property for all items below i. (negative indexes)
template <typename T> void maxSortDown(T *data, Mediator *m, int i) {
  for (i *= 2; i >= -m->maxCt; i *= 2) {
    if (i > -m->maxCt && mmless(data, m, i, i - 1)) {
      --i;
    }
    if (!mmCmpExch(data, m, i / 2, i)) {
      break;
    }
  }
}

// maintains minheap property for all items above i, including the rank
// returns true if rank changed
template <typename T> inline int minSortUp(T *data, Mediator *m, int i) {
  while (i > 0 && mmCmpExch(data, m, i, i / 2))
    i /= 2;
  return (i == 0);
}

// maintains maxheap property for all items above i, including rank
// returns true if rank changed
template <typename T> inline int maxSortUp(T *data, Mediator *m, int i) {
  while (i < 0 && mmCmpExch(data, m, i / 2, i))
    i /= 2;
  return (i == 0);
}

/*--- Public Interface ---*/

// Inserts item, maintains rank in O(lg nItems)
template <typename T> void MediatorInsert(T *data, Mediator *m, T v) {
  int p = m->pos[m->idx];
  T old = data[m->idx];
  data[m->idx] = v;
  m->idx++;
  if (m->idx == m->N) {
    m->idx = 0;
  }

  if (p > 0) // new item is in minHeap
  {
    if (v > old) {
      minSortDown(data, m, p);
      return;
    }
    if (minSortUp(data, m, p) && mmCmpExch(data, m, 0, -1)) {
      maxSortDown(data, m, -1);
    }
  } else if (p < 0) // new item is in maxheap
  {
    if (v < old) {
      maxSortDown(data, m, p);
      return;
    }
    if (maxSortUp(data, m, p) && mmCmpExch(data, m, 1, 0)) {
      minSortDown(data, m, 1);
    }
  } else // new item is at rank
  {
    if (maxSortUp(data, m, -1)) {
      maxSortDown(data, m, -1);
    }
    if (minSortUp(data, m, 1)) {
      minSortDown(data, m, 1);
    }
  }
}

//
// _rank_filter() requires the allocation of memory.  If the allocation fails,
// the function returns -1.  Otherwise it returns 0.
//
template <typename T>
int _rank_filter(T *in_arr, int rank, int arr_len, int win_len, T *out_arr,
                  int mode, T cval, int origin) {
  int i, arr_len_thresh, lim = (win_len - 1) / 2 - origin;
  int lim2 = arr_len - lim;
  /* Note: `arr_len == 1` is the only case implemented here for `lim2 < 0`; the calling code */
  /* in _filters.py ensures that this function isn't called otherwise. xref gh-23293 for details. */
  if (lim2 < 0 && arr_len == 1) {
      switch (mode) {
          case REFLECT:
          case NEAREST:
          case WRAP:
          case MIRROR:
              out_arr[0] = in_arr[0];
              return 0;
          case CONSTANT:
              if (win_len == 1) {
                  out_arr[0] = in_arr[0];
              }
              else {
                  out_arr[0] = cval;
              }
              return 0;
      }
  }
  int offset;
  std::unique_ptr<Mediator> m;
  std::unique_ptr<T[]> data;
  try {
    m = std::make_unique<Mediator>(win_len, rank);
    data = std::unique_ptr<T[]>(new T[win_len]());
  } catch (std::bad_alloc&) {
    return -1;
  }

  switch (mode) {
  case REFLECT:
    for (i = win_len - lim - 1; i > -1; i--) {
      MediatorInsert(data.get(), m.get(), in_arr[i]);
    }
    break;
  case CONSTANT:
    for (i = win_len - lim; i > 0; i--) {
      MediatorInsert(data.get(), m.get(), cval);
    }
    break;
  case NEAREST:
    for (i = win_len - lim; i > 0; i--) {
      MediatorInsert(data.get(), m.get(), in_arr[0]);
    }
    break;
  case MIRROR:
    for (i = win_len - lim; i > 0; i--) {
      MediatorInsert(data.get(), m.get(), in_arr[i]);
    }
    break;
  case WRAP:
    if (win_len % 2 == 0) {
        offset = 2;
    }
    else {
        offset = 0;
    }
    for (i = arr_len - lim - offset - 2 * origin; i < arr_len; i++) {
      MediatorInsert(data.get(), m.get(), in_arr[i]);
    }
    break;
  }

  for (i = 0; i < lim; i++) {
    MediatorInsert(data.get(), m.get(), in_arr[i]);
  }
  for (i = lim; i < arr_len; i++) {
    MediatorInsert(data.get(), m.get(), in_arr[i]);
    out_arr[i - lim] = data[m->heap[0]];
  }
  switch (mode) {
  case REFLECT:
    arr_len_thresh = arr_len - 1;
    for (i = 0; i < lim; i++) {
      MediatorInsert(data.get(), m.get(), in_arr[arr_len_thresh - i]);
      out_arr[lim2 + i] = data[m->heap[0]];
    }
    break;
  case CONSTANT:
    for (i = 0; i < lim; i++) {
      MediatorInsert(data.get(), m.get(), cval);
      out_arr[lim2 + i] = data[m->heap[0]];
    }
    break;
  case NEAREST:
    arr_len_thresh = arr_len - 1;
    for (i = 0; i < lim; i++) {
      MediatorInsert(data.get(), m.get(), in_arr[arr_len_thresh]);
      out_arr[lim2 + i] = data[m->heap[0]];
    }
    break;
  case MIRROR:
    arr_len_thresh = arr_len - 2;
    for (i = 0; i < lim; i++) {
      MediatorInsert(data.get(), m.get(), in_arr[arr_len_thresh - i]);
      out_arr[lim2 + i] = data[m->heap[0]];
    }
    break;
  case WRAP:
    for (i = 0; i < lim; i++) {
      MediatorInsert(data.get(), m.get(), in_arr[i]);
      out_arr[lim2 + i] = data[m->heap[0]];
    }
    break;
  }
  return 0;
}


// Python wrapper for rank_filter
static PyObject *rank_filter(PyObject *self, PyObject *args) {
  PyObject *in_arr_obj, *out_arr_obj, *cval_obj;
  int rank, arr_len, win_len, mode, origin, type;
  int rank_filter_status = 0;

  if (!PyArg_ParseTuple(args, "OiiOiOi", &in_arr_obj, &rank, &win_len,
                        &out_arr_obj, &mode, &cval_obj, &origin)) {
    return NULL;
  }

  PyArrayObject *in_arr = (PyArrayObject *)PyArray_FROM_OTF(
      in_arr_obj, NPY_NOTYPE, NPY_ARRAY_IN_ARRAY);
  if (in_arr == NULL) {
    return NULL;
  }
  PyArrayObject *out_arr = (PyArrayObject *)PyArray_FROM_OTF(
      out_arr_obj, NPY_NOTYPE, NPY_ARRAY_INOUT_ARRAY2);
  if (out_arr == NULL) {
    Py_DECREF(in_arr);
    return NULL;
  }

  arr_len = PyArray_SIZE(in_arr);
  type = PyArray_TYPE(in_arr);

  switch (type) { // the considered types are float, double, int64
  case NPY_FLOAT: {
    float *c_in_arr = (float *)PyArray_DATA(in_arr);
    float *c_out_arr = (float *)PyArray_DATA(out_arr);
    float cval = (float)PyFloat_AsDouble(cval_obj);
    rank_filter_status = _rank_filter(c_in_arr, rank, arr_len, win_len, c_out_arr,
                                      mode, cval, origin);
    break;
  }
  case NPY_DOUBLE: {
    double *c_in_arr = (double *)PyArray_DATA(in_arr);
    double *c_out_arr = (double *)PyArray_DATA(out_arr);
    double cval = PyFloat_AsDouble(cval_obj);
    rank_filter_status = _rank_filter(c_in_arr, rank, arr_len, win_len, c_out_arr,
                                      mode, cval, origin);
    break;
  }
  case NPY_INT64: {
    int64_t *c_in_arr = (int64_t *)PyArray_DATA(in_arr);
    int64_t *c_out_arr = (int64_t *)PyArray_DATA(out_arr);
    int64_t cval = PyLong_AsLongLong(cval_obj);
    rank_filter_status = _rank_filter(c_in_arr, rank, arr_len, win_len, c_out_arr,
                                      mode, cval, origin);
    break;
  }
  default:
    PyErr_SetString(PyExc_TypeError, "Unsupported array type");
    break;
  }
  if (rank_filter_status == -1) {
    PyErr_SetString(PyExc_MemoryError, "failed to allocate memory for rank filter");
  }
  Py_DECREF(in_arr);
  Py_DECREF(out_arr);
  if (PyErr_Occurred()) {
    return NULL;
  }
  Py_RETURN_NONE;
}

// define the module methods
static PyMethodDef myMethods[] = {
    {"rank_filter", rank_filter, METH_VARARGS, "1D rank filter"},
    {NULL, NULL, 0, NULL}};

// define the module
static struct PyModuleDef _rank_filter_1d = {
    PyModuleDef_HEAD_INIT, "_rank_filter_1d", "1D rank filter module", -1,
    myMethods};

// init the module
PyMODINIT_FUNC
PyInit__rank_filter_1d(void)
{
    PyObject *module;

    import_array();
    module = PyModule_Create(&_rank_filter_1d);
    if (module == NULL) {
        return module;
    }

#if Py_GIL_DISABLED
    PyUnstable_Module_SetGIL(module, Py_MOD_GIL_NOT_USED);
#endif

  return module;
}
