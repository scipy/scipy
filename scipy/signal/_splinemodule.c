#include "_splinemodule.h"

static void
convert_strides(npy_intp* instrides,npy_intp* convstrides,int size,int N)
{
  int n; npy_intp bitshift;

  bitshift = -1;

  while (size != 0) {
    size >>= 1;
    bitshift++;
  }
  for (n = 0; n < N; n++) {
    convstrides[n] = instrides[n] >> bitshift;
  }
}


static char doc_cspline2d[] = "out = cspline2d(input, lambda=0.0, precision=-1.0)\n"
"\n"
"    Coefficients for 2-D cubic (3rd order) B-spline.\n"
"\n"
"    Return the third-order B-spline coefficients over a regularly spaced\n"
"    input grid for the two-dimensional input image.\n"
"\n"
"    Parameters\n"
"    ----------\n"
"    input : ndarray\n"
"        The input signal.\n"
"    lambda : float\n"
"        Specifies the amount of smoothing in the transfer function.\n"
"    precision : float\n"
"        Specifies the precision for computing the infinite sum needed to apply mirror-\n"
"        symmetric boundary conditions.\n"
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


static PyObject *cspline2d(PyObject *NPY_UNUSED(dummy), PyObject *args)
{
  PyObject *image=NULL;
  PyArrayObject *a_image=NULL, *ck=NULL;
  double lambda = 0.0;
  double precision = -1.0;
  int thetype, M, N, retval=0;
  npy_intp outstrides[2], instrides[2];

  if (!PyArg_ParseTuple(args, "O|dd", &image, &lambda, &precision)) return NULL;

  thetype = PyArray_ObjectType(image, NPY_FLOAT);
  thetype = PyArray_MIN(thetype, NPY_DOUBLE);
  a_image = (PyArrayObject *)PyArray_FromObject(image, thetype, 2, 2);
  if (a_image == NULL) goto fail;

  ck = (PyArrayObject *)PyArray_SimpleNew(2, PyArray_DIMS(a_image), thetype);
  if (ck == NULL) goto fail;
  M = PyArray_DIMS(a_image)[0];
  N = PyArray_DIMS(a_image)[1];

  convert_strides(PyArray_STRIDES(a_image), instrides, PyArray_ITEMSIZE(a_image), 2);
  outstrides[0] = N;
  outstrides[1] = 1;

  if (thetype == NPY_FLOAT) {
    if ((precision <= 0.0) || (precision > 1.0)) precision = 1e-3;
    retval = S_cubic_spline2D((float *)PyArray_DATA(a_image),
                              (float *)PyArray_DATA(ck),
                              M, N, lambda, instrides, outstrides, precision);
  }
  else if (thetype == NPY_DOUBLE) {
    if ((precision <= 0.0) || (precision > 1.0)) precision = 1e-6;
    retval = D_cubic_spline2D((double *)PyArray_DATA(a_image),
                              (double *)PyArray_DATA(ck),
                              M, N, lambda, instrides, outstrides, precision);
  }

  if (retval == -3) PYERR("Precision too high.  Error did not converge.");
  if (retval < 0) PYERR("Problem occurred inside routine");

  Py_DECREF(a_image);
  return PyArray_Return(ck);

 fail:
  Py_XDECREF(a_image);
  Py_XDECREF(ck);
  return NULL;

}

static char doc_qspline2d[] = "out = qspline2d(input, lambda=0.0, precision=-1.0)\n"
"\n"
"    Coefficients for 2-D quadratic (2nd order) B-spline:\n"
"\n"
"    Return the second-order B-spline coefficients over a regularly spaced\n"
"    input grid for the two-dimensional input image.\n"
"\n"
"    Parameters\n"
"    ----------\n"
"    input : ndarray\n"
"        The input signal.\n"
"    lambda : float\n"
"        Specifies the amount of smoothing in the transfer function.\n"
"    precision : float\n"
"        Specifies the precision for computing the infinite sum needed to apply mirror-\n"
"        symmetric boundary conditions.\n"
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

static PyObject *qspline2d(PyObject *NPY_UNUSED(dummy), PyObject *args)
{
  PyObject *image=NULL;
  PyArrayObject *a_image=NULL, *ck=NULL;
  double lambda = 0.0;
  double precision = -1.0;
  int thetype, M, N, retval=0;
  npy_intp outstrides[2], instrides[2];

  if (!PyArg_ParseTuple(args, "O|dd", &image, &lambda, &precision)) return NULL;

  if (lambda != 0.0) PYERR("Smoothing spline not yet implemented.");

  thetype = PyArray_ObjectType(image, NPY_FLOAT);
  thetype = PyArray_MIN(thetype, NPY_DOUBLE);
  a_image = (PyArrayObject *)PyArray_FromObject(image, thetype, 2, 2);
  if (a_image == NULL) goto fail;

  ck = (PyArrayObject *)PyArray_SimpleNew(2, PyArray_DIMS(a_image), thetype);
  if (ck == NULL) goto fail;
  M = PyArray_DIMS(a_image)[0];
  N = PyArray_DIMS(a_image)[1];

  convert_strides(PyArray_STRIDES(a_image), instrides, PyArray_ITEMSIZE(a_image), 2);
  outstrides[0] = N;
  outstrides[1] = 1;

  if (thetype == NPY_FLOAT) {
    if ((precision <= 0.0) || (precision > 1.0)) precision = 1e-3;
    retval = S_quadratic_spline2D((float *)PyArray_DATA(a_image),
                                  (float *)PyArray_DATA(ck),
                                  M, N, lambda, instrides, outstrides, precision);
  }
  else if (thetype == NPY_DOUBLE) {
    if ((precision <= 0.0) || (precision > 1.0)) precision = 1e-6;
    retval = D_quadratic_spline2D((double *)PyArray_DATA(a_image),
                                  (double *)PyArray_DATA(ck),
                                  M, N, lambda, instrides, outstrides, precision);
  }

  if (retval == -3) PYERR("Precision too high.  Error did not converge.");
  if (retval < 0) PYERR("Problem occurred inside routine");

  Py_DECREF(a_image);
  return PyArray_Return(ck);

 fail:
  Py_XDECREF(a_image);
  Py_XDECREF(ck);
  return NULL;

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

static PyObject *FIRsepsym2d(PyObject *NPY_UNUSED(dummy), PyObject *args)
{
  PyObject *image=NULL, *hrow=NULL, *hcol=NULL;
  PyArrayObject *a_image=NULL, *a_hrow=NULL, *a_hcol=NULL, *out=NULL;
  int thetype, M, N, ret;
  npy_intp outstrides[2], instrides[2];

  if (!PyArg_ParseTuple(args, "OOO", &image, &hrow, &hcol)) return NULL;

  thetype = PyArray_ObjectType(image, NPY_FLOAT);
  thetype = PyArray_MIN(thetype, NPY_CDOUBLE);
  a_image = (PyArrayObject *)PyArray_FromObject(image, thetype, 2, 2);
  if (a_image == NULL) goto fail;

  a_hrow = (PyArrayObject *)PyArray_ContiguousFromObject(hrow, thetype, 1, 1);
  if (a_hrow == NULL) goto fail;

  a_hcol = (PyArrayObject *)PyArray_ContiguousFromObject(hcol, thetype, 1, 1);
  if (a_hcol == NULL) goto fail;

  out = (PyArrayObject *)PyArray_SimpleNew(2, PyArray_DIMS(a_image), thetype);
  if (out == NULL) goto fail;

  M = PyArray_DIMS(a_image)[0];
  N = PyArray_DIMS(a_image)[1];

  convert_strides(PyArray_STRIDES(a_image), instrides, PyArray_ITEMSIZE(a_image), 2);
  outstrides[0] = N;
  outstrides[1] = 1;

  if (PyArray_DIMS(a_hrow)[0] % 2 != 1 ||
      PyArray_DIMS(a_hcol)[0] % 2 != 1) {
    PYERR("hrow and hcol must be odd length");
  }

  switch (thetype) {
  case NPY_FLOAT:
    ret = S_separable_2Dconvolve_mirror((float *)PyArray_DATA(a_image),
					(float *)PyArray_DATA(out), M, N,
					(float *)PyArray_DATA(a_hrow),
					(float *)PyArray_DATA(a_hcol),
					PyArray_DIMS(a_hrow)[0], PyArray_DIMS(a_hcol)[0],
					instrides, outstrides);
    break;
  case NPY_DOUBLE:
    ret = D_separable_2Dconvolve_mirror((double *)PyArray_DATA(a_image),
					(double *)PyArray_DATA(out), M, N,
					(double *)PyArray_DATA(a_hrow),
					(double *)PyArray_DATA(a_hcol),
					PyArray_DIMS(a_hrow)[0], PyArray_DIMS(a_hcol)[0],
					instrides, outstrides);
    break;
#ifdef __GNUC__
  case NPY_CFLOAT:
    ret = C_separable_2Dconvolve_mirror((__complex__ float *)PyArray_DATA(a_image),
					(__complex__ float *)PyArray_DATA(out), M, N,
					(__complex__ float *)PyArray_DATA(a_hrow),
					(__complex__ float *)PyArray_DATA(a_hcol),
					PyArray_DIMS(a_hrow)[0], PyArray_DIMS(a_hcol)[0],
					instrides, outstrides);
    break;
  case NPY_CDOUBLE:
    ret = Z_separable_2Dconvolve_mirror((__complex__ double *)PyArray_DATA(a_image),
					(__complex__ double *)PyArray_DATA(out), M, N,
					(__complex__ double *)PyArray_DATA(a_hrow),
					(__complex__ double *)PyArray_DATA(a_hcol),
					PyArray_DIMS(a_hrow)[0], PyArray_DIMS(a_hcol)[0],
					instrides, outstrides);
    break;
#endif
  default:
    PYERR("Incorrect type.");
  }

  if (ret < 0) PYERR("Problem occurred inside routine.");

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

static char doc_IIRsymorder1[] = "out = symiirorder1(input, c0, z1, precision=-1.0)\n"
"\n"
"    Implement a smoothing IIR filter with mirror-symmetric boundary conditions\n"
"    using a cascade of first-order sections.  The second section uses a\n"
"    reversed sequence.  This implements a system with the following\n"
"    transfer function and mirror-symmetric boundary conditions::\n"
"\n"
"                           c0              \n"
"           H(z) = ---------------------    \n"
"                   (1-z1/z) (1 - z1 z)     \n"
"\n"
"    The resulting signal will have mirror symmetric boundary conditions as well.\n"
"\n"
"    Parameters\n"
"    ----------\n"
"    input : ndarray\n"
"        The input signal.\n"
"    c0, z1 : scalar\n"
"        Parameters in the transfer function.\n"
"    precision :\n"
"        Specifies the precision for calculating initial conditions\n"
"        of the recursive filter based on mirror-symmetric input.\n"
"\n"
"    Returns\n"
"    -------\n"
"    output : ndarray\n"
"        The filtered signal.";

static PyObject *IIRsymorder1(PyObject *NPY_UNUSED(dummy), PyObject *args)
{
  PyObject *sig=NULL;
  PyArrayObject *a_sig=NULL, *out=NULL;
  Py_complex c0, z1;
  double precision = -1.0;
  int thetype, N, ret;
  npy_intp outstrides, instrides;

  if (!PyArg_ParseTuple(args, "ODD|d", &sig, &c0, &z1, &precision))
    return NULL;

  thetype = PyArray_ObjectType(sig, NPY_FLOAT);
  thetype = PyArray_MIN(thetype, NPY_CDOUBLE);
  a_sig = (PyArrayObject *)PyArray_FromObject(sig, thetype, 1, 1);

  if (a_sig == NULL) goto fail;

  out = (PyArrayObject *)PyArray_SimpleNew(1, PyArray_DIMS(a_sig), thetype);
  if (out == NULL) goto fail;
  N = PyArray_DIMS(a_sig)[0];

  convert_strides(PyArray_STRIDES(a_sig), &instrides, PyArray_ITEMSIZE(a_sig), 1);
  outstrides = 1;

  switch (thetype) {
  case NPY_FLOAT:
    {
      float rc0 = c0.real;
      float rz1 = z1.real;

      if ((precision <= 0.0) || (precision > 1.0)) precision = 1e-6;
      ret = S_IIR_forback1 (rc0, rz1, (float *)PyArray_DATA(a_sig),
			    (float *)PyArray_DATA(out), N,
			    instrides, outstrides, (float )precision);
    }
    break;
  case NPY_DOUBLE:
    {
      double rc0 = c0.real;
      double rz1 = z1.real;

      if ((precision <= 0.0) || (precision > 1.0)) precision = 1e-11;
      ret = D_IIR_forback1 (rc0, rz1, (double *)PyArray_DATA(a_sig),
			    (double *)PyArray_DATA(out), N,
			    instrides, outstrides, precision);
    }
    break;
#ifdef __GNUC__
  case NPY_CFLOAT:
    {
      __complex__ float zc0 = c0.real + 1.0i*c0.imag;
      __complex__ float zz1 = z1.real + 1.0i*z1.imag;
      if ((precision <= 0.0) || (precision > 1.0)) precision = 1e-6;
      ret = C_IIR_forback1 (zc0, zz1, (__complex__ float *)PyArray_DATA(a_sig),
			    (__complex__ float *)PyArray_DATA(out), N,
			    instrides, outstrides, (float )precision);
    }
    break;
  case NPY_CDOUBLE:
    {
      __complex__ double zc0 = c0.real + 1.0i*c0.imag;
      __complex__ double zz1 = z1.real + 1.0i*z1.imag;
      if ((precision <= 0.0) || (precision > 1.0)) precision = 1e-11;
      ret = Z_IIR_forback1 (zc0, zz1, (__complex__ double *)PyArray_DATA(a_sig),
			    (__complex__ double *)PyArray_DATA(out), N,
			    instrides, outstrides, precision);
    }
    break;
#endif
  default:
    PYERR("Incorrect type.");
  }

  if (ret == 0) {
    Py_DECREF(a_sig);
    return PyArray_Return(out);
  }

  if (ret == -1) PYERR("Could not allocate enough memory.");
  if (ret == -2) PYERR("|z1| must be less than 1.0");
  if (ret == -3) PYERR("Sum to find symmetric boundary conditions did not converge.");

  PYERR("Unknown error.");


 fail:
  Py_XDECREF(a_sig);
  Py_XDECREF(out);
  return NULL;

}

static char doc_IIRsymorder2[] = "out = symiirorder2(input, r, omega, precision=-1.0)\n"
"\n"
"    Implement a smoothing IIR filter with mirror-symmetric boundary conditions\n"
"    using a cascade of second-order sections.  The second section uses a\n"
"    reversed sequence.  This implements the following transfer function::\n"
"\n"
"                                  cs^2\n"
"         H(z) = ---------------------------------------\n"
"                (1 - a2/z - a3/z^2) (1 - a2 z - a3 z^2 )\n"
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
"    output : ndarray\n"
"        The filtered signal.";

static PyObject *IIRsymorder2(PyObject *NPY_UNUSED(dummy), PyObject *args)
{
  PyObject *sig=NULL;
  PyArrayObject *a_sig=NULL, *out=NULL;
  double r, omega;
  double precision = -1.0;
  int thetype, N, ret;
  npy_intp outstrides, instrides;

  if (!PyArg_ParseTuple(args, "Odd|d", &sig, &r, &omega, &precision))
    return NULL;

  thetype = PyArray_ObjectType(sig, NPY_FLOAT);
  thetype = PyArray_MIN(thetype, NPY_DOUBLE);
  a_sig = (PyArrayObject *)PyArray_FromObject(sig, thetype, 1, 1);

  if (a_sig == NULL) goto fail;

  out = (PyArrayObject *)PyArray_SimpleNew(1, PyArray_DIMS(a_sig), thetype);
  if (out == NULL) goto fail;
  N = PyArray_DIMS(a_sig)[0];

  convert_strides(PyArray_STRIDES(a_sig), &instrides, PyArray_ITEMSIZE(a_sig), 1);
  outstrides = 1;

  switch (thetype) {
  case NPY_FLOAT:
    if ((precision <= 0.0) || (precision > 1.0)) precision = 1e-6;
    ret = S_IIR_forback2 (r, omega, (float *)PyArray_DATA(a_sig),
			  (float *)PyArray_DATA(out), N,
			  instrides, outstrides, precision);
    break;
  case NPY_DOUBLE:
    if ((precision <= 0.0) || (precision > 1.0)) precision = 1e-11;
    ret = D_IIR_forback2 (r, omega, (double *)PyArray_DATA(a_sig),
			  (double *)PyArray_DATA(out), N,
			  instrides, outstrides, precision);
    break;
  default:
    PYERR("Incorrect type.");
  }

  if (ret < 0) PYERR("Problem occurred inside routine.");

  Py_DECREF(a_sig);
  return PyArray_Return(out);

 fail:
  Py_XDECREF(a_sig);
  Py_XDECREF(out);
  return NULL;

}


static struct PyMethodDef toolbox_module_methods[] = {
    {"cspline2d", cspline2d, METH_VARARGS, doc_cspline2d},
    {"qspline2d", qspline2d, METH_VARARGS, doc_qspline2d},
    {"sepfir2d", FIRsepsym2d, METH_VARARGS, doc_FIRsepsym2d},
    {"symiirorder1", IIRsymorder1, METH_VARARGS, doc_IIRsymorder1},
    {"symiirorder2", IIRsymorder2, METH_VARARGS, doc_IIRsymorder2},
    {NULL, NULL, 0, NULL}		/* sentinel */
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

PyMODINIT_FUNC
PyInit__spline(void)
{
    import_array();
    return PyModule_Create(&moduledef);
}
