#include "_hampel_filter.h"
#include "_rank_filter_1d.h"

// Python wrapper for rank_filter
PyObject *rank_filter(PyObject *self, PyObject *args) {
  PyObject *in_arr_obj, *out_arr_obj, *cval_obj;
  int rank, arr_len, win_len, mode, origin, type;
  if (!PyArg_ParseTuple(args, "OiiOiOi", &in_arr_obj, &rank, &win_len,
                        &out_arr_obj, &mode, &cval_obj, &origin)) {
    return NULL;
  }
  PyArrayObject *in_arr = (PyArrayObject *)PyArray_FROM_OTF(
      in_arr_obj, NPY_NOTYPE, NPY_ARRAY_INOUT_ARRAY2);
  PyArrayObject *out_arr = (PyArrayObject *)PyArray_FROM_OTF(
      out_arr_obj, NPY_NOTYPE, NPY_ARRAY_INOUT_ARRAY2);

  if (in_arr == NULL || out_arr == NULL) {
    return NULL;
  }
  arr_len = PyArray_SIZE(in_arr);
  type = PyArray_TYPE(in_arr);
  switch (type) { // the considered types are float, double, int64
  case NPY_FLOAT: {
    float *c_in_arr = (float *)PyArray_DATA(in_arr);
    float *c_out_arr = (float *)PyArray_DATA(out_arr);
    float cval = (float)PyFloat_AsDouble(cval_obj);
    _rank_filter(c_in_arr, rank, arr_len, win_len, c_out_arr, mode, cval,
                 origin);
    break;
  }
  case NPY_DOUBLE: {
    double *c_in_arr = (double *)PyArray_DATA(in_arr);
    double *c_out_arr = (double *)PyArray_DATA(out_arr);
    double cval = PyFloat_AsDouble(cval_obj);
    _rank_filter(c_in_arr, rank, arr_len, win_len, c_out_arr, mode, cval,
                 origin);
    break;
  }
  case NPY_INT64: {
    int64_t *c_in_arr = (int64_t *)PyArray_DATA(in_arr);
    int64_t *c_out_arr = (int64_t *)PyArray_DATA(out_arr);
    int64_t cval = PyLong_AsLongLong(cval_obj);
    _rank_filter(c_in_arr, rank, arr_len, win_len, c_out_arr, mode, cval,
                 origin);
    break;
  }
  default:
    PyErr_SetString(PyExc_TypeError, "Unsupported array type");
    break;
  }
  Py_DECREF(in_arr);
  Py_DECREF(out_arr);
  Py_RETURN_NONE;
}

// Python wrapper for rank_filter
PyObject* hampel_filter(PyObject* self, PyObject* args)
{
    PyObject *in_arr_obj, *out_arr_obj, *median_obj, *mad_obj, *scale_obj;
    int arr_len, win_len, type;
    if (!PyArg_ParseTuple(args, "OiOOOO", &in_arr_obj, &win_len, &median_obj, &mad_obj, &out_arr_obj, &scale_obj))
    {
        return NULL;
    }
    PyArrayObject *in_arr = (PyArrayObject *)PyArray_FROM_OTF(in_arr_obj, NPY_NOTYPE, NPY_ARRAY_INOUT_ARRAY2);
    PyArrayObject *out_arr = (PyArrayObject *)PyArray_FROM_OTF(out_arr_obj, NPY_NOTYPE, NPY_ARRAY_INOUT_ARRAY2);
    PyArrayObject *median = (PyArrayObject *)PyArray_FROM_OTF(median_obj, NPY_NOTYPE, NPY_ARRAY_INOUT_ARRAY2);
    PyArrayObject *mad = (PyArrayObject *)PyArray_FROM_OTF(mad_obj, NPY_NOTYPE, NPY_ARRAY_INOUT_ARRAY2);
    if (in_arr == NULL || out_arr == NULL || median == NULL || mad == NULL)
    {
        return NULL;
    }
    arr_len = PyArray_SIZE(in_arr);
    type = PyArray_TYPE(in_arr);
    switch (type) { // the considered types are float, double, int64
        case NPY_FLOAT:
        {
            float *c_in_arr = (float *)PyArray_DATA(in_arr);
            float *c_out_arr = (float *)PyArray_DATA(out_arr);
            float *c_median = (float *)PyArray_DATA(median);
            float *c_mad = (float *)PyArray_DATA(mad);
            float c_scale = (float)PyFloat_AsDouble(scale_obj);
            _hampel_filter(c_in_arr, arr_len, win_len, c_median, c_mad, c_out_arr, c_scale);
            break;
        }
        case NPY_DOUBLE:
        {
            double *c_in_arr = (double *)PyArray_DATA(in_arr);
            double *c_out_arr = (double *)PyArray_DATA(out_arr);
            double *c_median = (double *)PyArray_DATA(median);
            double *c_mad = (double *)PyArray_DATA(mad);
            double c_scale = PyFloat_AsDouble(scale_obj);
            _hampel_filter(c_in_arr, arr_len, win_len, c_median, c_mad, c_out_arr, c_scale);
            break;
        }
        case NPY_INT64:
        {
            int64_t *c_in_arr = (int64_t *)PyArray_DATA(in_arr);
            int64_t *c_out_arr = (int64_t *)PyArray_DATA(out_arr);
            int64_t *c_median = (int64_t *)PyArray_DATA(median);
            int64_t *c_mad = (int64_t *)PyArray_DATA(mad);
            int64_t c_scale = PyLong_AsLong(scale_obj);
            _hampel_filter(c_in_arr, arr_len, win_len, c_median, c_mad, c_out_arr, c_scale);
            break;
        }
        default:
            PyErr_SetString(PyExc_TypeError, "Unsupported array type");
            break;
    }
    Py_DECREF(in_arr);
    Py_DECREF(out_arr);
    Py_RETURN_NONE;
}

// define the module methods
static PyMethodDef myMethods[] = {
    {"rank_filter", rank_filter, METH_VARARGS, "1D rank filter"},
    {"hampel_filter", hampel_filter, METH_VARARGS, "Hampel filter"},
    {NULL, NULL, 0, NULL}};

//define the module
static struct PyModuleDef _robust_filters_1d = {
    PyModuleDef_HEAD_INIT,
    "_robust_filters_1d",
    "1D robust filters",
    -1,
    myMethods};

//init the module
PyMODINIT_FUNC PyInit__robust_filters_1d(void)
{
    import_array();
    return PyModule_Create(&_robust_filters_1d);
}