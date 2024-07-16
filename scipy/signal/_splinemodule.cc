#include "Python.h"
#include "numpy/arrayobject.h"
#include <complex>
#include <cmath>
#include <algorithm>
#include "_splinemodule.h"


static void convert_strides(npy_intp *instrides, npy_intp *convstrides, int size, int N) {
    int n;
    npy_intp bitshift;

  bitshift = -1;

  while (size != 0) {
    size >>= 1;
    bitshift++;
  }
  for (n = 0; n < N; n++) {
    convstrides[n] = instrides[n] >> bitshift;
  }
}


static char doc_FIRsepsym2d[] = "out = sepfir2d(input, hrow, hcol)\n"
"\n"
"    Convolve with a 2-D separable FIR filter.\n"
"\n"
"    Convolve the rank-2 input array with the separable filter defined by the\n"
"    rank-1 arrays hrow, and hcol. Mirror symmetric boundary conditions are\n"
"    assumed. This function can be used to find an image given its B-spline\n"
"    representation.\n"
"\n"
"    Parameters\n"
"    ----------\n"
"    input : ndarray\n"
"        The input signal. Must be a rank-2 array.\n"
"    hrow : ndarray\n"
"        A rank-1 array defining the row direction of the filter.\n"
"        Must be odd-length\n"
"    hcol : ndarray\n"
"        A rank-1 array defining the column direction of the filter.\n"
"        Must be odd-length\n"
"\n"
"    Returns\n"
"    -------\n"
"    output : ndarray\n"
"        The filtered signal.\n"
"\n"
"    Examples\n"
"    --------\n"
"    Examples are given :ref:`in the tutorial <tutorial-signal-bsplines>`.\n"
"\n";

static PyObject *FIRsepsym2d(PyObject *NPY_UNUSED(dummy), PyObject *args) {
  PyObject *image=NULL, *hrow=NULL, *hcol=NULL;
  PyArrayObject *a_image=NULL, *a_hrow=NULL, *a_hcol=NULL, *out=NULL;
  int thetype, M, N, ret;
  npy_intp outstrides[2], instrides[2];

  if (!PyArg_ParseTuple(args, "OOO", &image, &hrow, &hcol))
      return NULL;

  thetype = PyArray_ObjectType(image, NPY_FLOAT);
  thetype = std::min(thetype, static_cast<int>(NPY_CDOUBLE));
  a_image = (PyArrayObject *)PyArray_FromObject(image, thetype, 2, 2);
  if (a_image == NULL)
      goto fail;

  a_hrow = (PyArrayObject *)PyArray_ContiguousFromObject(hrow, thetype, 1, 1);
  if (a_hrow == NULL)
      goto fail;

  a_hcol = (PyArrayObject *)PyArray_ContiguousFromObject(hcol, thetype, 1, 1);
  if (a_hcol == NULL)
      goto fail;

  out = (PyArrayObject *)PyArray_SimpleNew(2, PyArray_DIMS(a_image), thetype);
  if (out == NULL)
      goto fail;

  M = PyArray_DIMS(a_image)[0];
  N = PyArray_DIMS(a_image)[1];

  convert_strides(PyArray_STRIDES(a_image), instrides, PyArray_ITEMSIZE(a_image), 2);
  outstrides[0] = N;
  outstrides[1] = 1;

  if (PyArray_DIMS(a_hrow)[0] % 2 != 1 || PyArray_DIMS(a_hcol)[0] % 2 != 1) {
    PYERR("hrow and hcol must be odd length");
  }

  switch (thetype) {
  case NPY_FLOAT:
    ret = _separable_2Dconvolve_mirror((float *) PyArray_DATA(a_image), (float *) PyArray_DATA(out), M, N,
                                       (float *) PyArray_DATA(a_hrow), (float *) PyArray_DATA(a_hcol),
                                       PyArray_DIMS(a_hrow)[0], PyArray_DIMS(a_hcol)[0], instrides, outstrides);
    break;
  case NPY_DOUBLE:
    ret = _separable_2Dconvolve_mirror((double *) PyArray_DATA(a_image), (double *) PyArray_DATA(out), M, N,
                                       (double *) PyArray_DATA(a_hrow), (double *) PyArray_DATA(a_hcol),
                                       PyArray_DIMS(a_hrow)[0], PyArray_DIMS(a_hcol)[0], instrides, outstrides);
    break;
  case NPY_CFLOAT:
    ret = _separable_2Dconvolve_mirror(reinterpret_cast<std::complex<float> *>(PyArray_DATA(a_image)),
                                       reinterpret_cast<std::complex<float> *>(PyArray_DATA(out)), M, N,
                                       reinterpret_cast<std::complex<float> *>(PyArray_DATA(a_hrow)),
                                       reinterpret_cast<std::complex<float> *>(PyArray_DATA(a_hcol)),
                                       PyArray_DIMS(a_hrow)[0], PyArray_DIMS(a_hcol)[0], instrides, outstrides);
    break;
  case NPY_CDOUBLE:
    ret = _separable_2Dconvolve_mirror(reinterpret_cast<std::complex<double> *>(PyArray_DATA(a_image)),
                                       reinterpret_cast<std::complex<double> *>(PyArray_DATA(out)), M, N,
                                       reinterpret_cast<std::complex<double> *>(PyArray_DATA(a_hrow)),
                                       reinterpret_cast<std::complex<double> *>(PyArray_DATA(a_hcol)),
                                       PyArray_DIMS(a_hrow)[0], PyArray_DIMS(a_hcol)[0], instrides, outstrides);
    break;
  default:
    PYERR("Incorrect type.");
  }

  if (ret < 0)
      PYERR("Problem occurred inside routine.");

  Py_DECREF(a_image);
  Py_DECREF(a_hrow);
  Py_DECREF(a_hcol);
  return PyArray_Return(out);

 fail:
  Py_XDECREF(a_image);
  Py_XDECREF(a_hrow);
  Py_XDECREF(a_hcol);
  Py_XDECREF(out);
  return NULL;

}


static char doc_IIRsymorder1_ic[] = "out = symiirorder1_ic(input, z1, precision=-1.0)\n"
"\n"
"    Compute the (forward) mirror-symmetric boundary conditions for a smoothing\n"
"    IIR filter that is composed of cascaded first-order sections.\n"
"\n"
"    The starting condition returned by this function is computed based on\n"
"    the following transfer function::\n"
"\n"
"                       1         \n"
"           H(z) = ------------   \n"
"                   (1 - z1/z)    \n"
"\n"
"\n"
"    Parameters\n"
"    ----------\n"
"    input : ndarray\n"
"        The input signal. If 2D, then it will find the initial conditions \n"
"        for each of the elements on the last axis.\n"
"    z1 : scalar\n"
"        Parameter in the transfer function.\n"
"    precision :\n"
"        Specifies the precision for calculating initial conditions\n"
"        of the recursive filter based on mirror-symmetric input.\n"
"\n"
"    Returns\n"
"    -------\n"
"    z_0 : ndarray\n"
"        The mirror-symmetric initial condition for the forward IIR filter.";

static PyObject *IIRsymorder1_ic(PyObject *NPY_UNUSED(dummy), PyObject *args) {
  PyObject *sig=NULL;
  PyArrayObject *a_sig=NULL, *out=NULL;
  npy_intp* in_size;

  Py_complex z1;
  double precision = -1.0;
  int thetype, ret;
  npy_intp M, N;
  PyArray_Descr* dtype;

  if (!PyArg_ParseTuple(args, "OD|d", &sig, &z1, &precision))
    return NULL;

  thetype = PyArray_ObjectType(sig, NPY_FLOAT);
  thetype = std::min(thetype, static_cast<int>(NPY_CDOUBLE));
  a_sig = (PyArrayObject *)PyArray_FromObject(sig, thetype, 1, 2);

  if (a_sig == NULL) {
    // Duplicate the `goto fail` block to avoid "crosses initialization" build errors.
    Py_XDECREF(a_sig);
    Py_XDECREF(out);
    return NULL;
  }

  in_size = PyArray_DIMS(a_sig);
  M = 1;
  N = in_size[0];

  if(PyArray_NDIM(a_sig) > 1) {
    M = in_size[0];
    N = in_size[1];
  }

  const npy_intp sz[2] = {M, 1};
  dtype = PyArray_DescrFromType(thetype);
  out = (PyArrayObject *)PyArray_Empty(2, sz, dtype, 0);
  if (out == NULL) {
    // Duplicate the `goto fail` block to avoid "crosses initialization" build errors.
    Py_XDECREF(a_sig);
    Py_XDECREF(out);
    return NULL;
  }

  switch (thetype) {
    case NPY_FLOAT: {
      float rz1 = z1.real;

      if ((precision <= 0.0) || (precision > 1.0))
          precision = 1e-6;
      ret = _sym_iir1_initial(rz1, static_cast<float *>(PyArray_DATA(a_sig)), static_cast<float *>(PyArray_DATA(out)),
                              M, N, static_cast<float>(precision));
    } break;
    case NPY_DOUBLE: {
      double rz1 = z1.real;

      if ((precision <= 0.0) || (precision > 1.0))
          precision = 1e-11;
      ret = _sym_iir1_initial(rz1, static_cast<double *>(PyArray_DATA(a_sig)),
                            static_cast<double *>(PyArray_DATA(out)), M, N, precision);

    } break;
    case NPY_CFLOAT: {
      std::complex<float> zz1(z1.real, z1.imag);
      if ((precision <= 0.0) || (precision > 1.0))
          precision = 1e-6;
      ret = _sym_iir1_initial(zz1, reinterpret_cast<std::complex<float> *>(PyArray_DATA(a_sig)),
                              reinterpret_cast<std::complex<float> *>(PyArray_DATA(out)), M, N,
                              static_cast<float>(precision));

    } break;
    case NPY_CDOUBLE: {
      std::complex<double> zz1(z1.real, z1.imag);
      if ((precision <= 0.0) || (precision > 1.0))
          precision = 1e-11;
      ret = _sym_iir1_initial(zz1, reinterpret_cast<std::complex<double> *>(PyArray_DATA(a_sig)),
                              reinterpret_cast<std::complex<double> *>(PyArray_DATA(out)), M, N, precision);

    } break;
    default:
      PYERR("Incorrect type.");
  }

  if (ret == 0) {
      Py_DECREF(a_sig);
      return PyArray_Return(out);
  }

  if (ret == -1)
      PYERR("Could not allocate enough memory.");
  if (ret == -2)
      PYERR("|z1| must be less than 1.0");
  if (ret == -3)
      PYERR("Sum to find symmetric boundary conditions did not converge.");

  PYERR("Unknown error.");

 fail:
  Py_XDECREF(a_sig);
  Py_XDECREF(out);
  return NULL;

}

static char doc_IIRsymorder2_ic_fwd[] = "out = symiirorder2_ic_fwd(input, r, omega, precision=-1.0)\n"
"\n"
"    Compute the (forward) mirror-symmetric boundary conditions for a smoothing\n"
"    IIR filter that is composed of cascaded second-order sections.\n"
"\n"
"    The starting condition returned by this function is computed based on\n"
"    the following transfer function::\n"
"\n"
"                         cs\n"
"         H(z) = -------------------\n"
"                (1 - a2/z - a3/z^2)\n"
"\n"
"    where::\n"
"\n"
"          a2 = (2 r cos omega)\n"
"          a3 = - r^2\n"
"          cs = 1 - 2 r cos omega + r^2\n"
"\n"
"    Parameters\n"
"    ----------\n"
"    input : ndarray\n"
"        The input signal.\n"
"    r, omega : float\n"
"        Parameters in the transfer function.\n"
"    precision : float\n"
"        Specifies the precision for calculating initial conditions\n"
"        of the recursive filter based on mirror-symmetric input.\n"
"\n"
"    Returns\n"
"    -------\n"
"    zi : ndarray\n"
"        The mirror-symmetric initial condition for the forward IIR filter.";

static PyObject *IIRsymorder2_ic_fwd(PyObject *NPY_UNUSED(dummy), PyObject *args) {
  PyObject *sig=NULL;
  PyArrayObject *a_sig=NULL, *out=NULL;
  npy_intp* in_size;
  double r, omega;
  double precision = -1.0;
  int thetype, ret;
  npy_intp N, M;
  PyArray_Descr* dtype;

  if (!PyArg_ParseTuple(args, "Odd|d", &sig, &r, &omega, &precision))
    return NULL;

  thetype = PyArray_ObjectType(sig, NPY_FLOAT);
  thetype = std::min(thetype, static_cast<int>(NPY_DOUBLE));
  a_sig = (PyArrayObject *)PyArray_FromObject(sig, thetype, 1, 2);

  if (a_sig == NULL) {
    // Duplicate the `goto fail` block to avoid "crosses initialization" build errors.
    Py_XDECREF(a_sig);
    Py_XDECREF(out);
    return NULL;
  }

  in_size = PyArray_DIMS(a_sig);
  M = 1;
  N = in_size[0];

  if(PyArray_NDIM(a_sig) > 1) {
    M = in_size[0];
    N = in_size[1];
  }

  dtype = PyArray_DescrFromType(thetype);
  const npy_intp sz[2] = {M, 2};
  out = (PyArrayObject *)PyArray_Empty(2, sz, dtype, 0);
  if (out == NULL) {
    // Duplicate the `goto fail` block to avoid "crosses initialization" build errors.
    Py_XDECREF(a_sig);
    Py_XDECREF(out);
    return NULL;
  }

  switch (thetype) {
    case NPY_FLOAT: {
      if ((precision <= 0.0) || (precision > 1.0))
          precision = 1e-6;
      ret = _sym_iir2_initial_fwd<float>(r, omega, (float *) (PyArray_DATA(a_sig)), (float *) (PyArray_DATA(out)), M,
                                         N, (float) (precision));

    } break;
    case NPY_DOUBLE: {
      if ((precision <= 0.0) || (precision > 1.0))
          precision = 1e-11;
      ret = _sym_iir2_initial_fwd<double>(r, omega, static_cast<double *>(PyArray_DATA(a_sig)),
                                          static_cast<double *>(PyArray_DATA(out)), M, N, precision);
    } break;
    default:
      PYERR("Incorrect type.");
  }

  if (ret == 0) {
    Py_DECREF(a_sig);
    return PyArray_Return(out);
  }

  if (ret == -1)
      PYERR("Could not allocate enough memory.");
  if (ret == -2)
      PYERR("|z1| must be less than 1.0");
  if (ret == -3)
      PYERR("Sum to find symmetric boundary conditions did not converge.");

  PYERR("Unknown error.");


 fail:
  Py_XDECREF(a_sig);
  Py_XDECREF(out);
  return NULL;

}

static char doc_IIRsymorder2_ic_bwd[] = "out = symiirorder2_ic_bwd(input, r, omega, precision=-1.0)\n"
"\n"
"    Compute the (backward) mirror-symmetric boundary conditions for a smoothing\n"
"    IIR filter that is composed of cascaded second-order sections.\n"
"\n"
"    The starting condition returned by this function is computed based on\n"
"    the following transfer function::\n"
"\n"
"                         cs\n"
"         H(z) = -------------------\n"
"                (1 - a2 z - a3 z^2)\n"
"\n"
"    where::\n"
"\n"
"          a2 = (2 r cos omega)\n"
"          a3 = - r^2\n"
"          cs = 1 - 2 r cos omega + r^2\n"
"\n"
"    Parameters\n"
"    ----------\n"
"    input : ndarray\n"
"        The input signal.\n"
"    r, omega : float\n"
"        Parameters in the transfer function.\n"
"    precision : float\n"
"        Specifies the precision for calculating initial conditions\n"
"        of the recursive filter based on mirror-symmetric input.\n"
"\n"
"    Returns\n"
"    -------\n"
"    zi : ndarray\n"
"        The mirror-symmetric initial condition for the forward IIR filter.";

static PyObject *IIRsymorder2_ic_bwd(PyObject *NPY_UNUSED(dummy), PyObject *args) {
  PyObject *sig=NULL;
  PyArrayObject *a_sig=NULL, *out=NULL;
  npy_intp* in_size;
  double r, omega;
  double precision = -1.0;
  int thetype, ret;
  npy_intp M, N;
  PyArray_Descr* dtype;

  if (!PyArg_ParseTuple(args, "Odd|d", &sig, &r, &omega, &precision))
    return NULL;

  thetype = PyArray_ObjectType(sig, NPY_FLOAT);
  thetype = std::min(thetype, static_cast<int>(NPY_DOUBLE));
  a_sig = (PyArrayObject *)PyArray_FromObject(sig, thetype, 1, 2);

  if (a_sig == NULL) {
    // Duplicate the `goto fail` block to avoid "crosses initialization" build errors.
    Py_XDECREF(a_sig);
    Py_XDECREF(out);
    return NULL;
  }

  in_size = PyArray_DIMS(a_sig);
  M = 1;
  N = in_size[0];

  if(PyArray_NDIM(a_sig) > 1) {
    M = in_size[0];
    N = in_size[1];
  }

  dtype = PyArray_DescrFromType(thetype);
  const npy_intp sz[2] = {M, 2};
  out = (PyArrayObject *)PyArray_Zeros(2, sz, dtype, 0);
  if (out == NULL) {
    // Duplicate the `goto fail` block to avoid "crosses initialization" build errors.
    Py_XDECREF(a_sig);
    Py_XDECREF(out);
    return NULL;
  }

  switch (thetype) {
    case NPY_FLOAT: {
      if ((precision <= 0.0) || (precision > 1.0))
          precision = 1e-6;
      ret = _sym_iir2_initial_bwd(r, omega, static_cast<float *>(PyArray_DATA(a_sig)),
                                  static_cast<float *>(PyArray_DATA(out)), M, N, static_cast<float>(precision));
    } break;
    case NPY_DOUBLE: {
      if ((precision <= 0.0) || (precision > 1.0))
          precision = 1e-11;
      ret = _sym_iir2_initial_bwd(r, omega, static_cast<double *>(PyArray_DATA(a_sig)),
                                  static_cast<double *>(PyArray_DATA(out)), M, N, precision);
    } break;
    default:
      PYERR("Incorrect type.");
  }

  if (ret == 0) {
    Py_DECREF(a_sig);
    return PyArray_Return(out);
  }

  if (ret == -1)
      PYERR("Could not allocate enough memory.");
  if (ret == -2)
      PYERR("|z1| must be less than 1.0");
  if (ret == -3)
      PYERR("Sum to find symmetric boundary conditions did not converge.");

  PYERR("Unknown error.");


 fail:
  Py_XDECREF(a_sig);
  Py_XDECREF(out);
  return NULL;

}


static struct PyMethodDef toolbox_module_methods[] = {
    {"sepfir2d", FIRsepsym2d, METH_VARARGS, doc_FIRsepsym2d},
    {"symiirorder1_ic", IIRsymorder1_ic, METH_VARARGS, doc_IIRsymorder1_ic},
    {"symiirorder2_ic_fwd", IIRsymorder2_ic_fwd, METH_VARARGS, doc_IIRsymorder2_ic_fwd},
    {"symiirorder2_ic_bwd", IIRsymorder2_ic_bwd, METH_VARARGS,doc_IIRsymorder2_ic_bwd },
    {NULL, NULL, 0, NULL}       /* sentinel */
};

/* Initialization function for the module */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_spline",
    NULL,
    -1,
    toolbox_module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC PyInit__spline(void) {
    import_array();
    return PyModule_Create(&moduledef);
}
