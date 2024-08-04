#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "numpy/arrayobject.h"

#ifdef HAVE_BLAS_ILP64

#define F_INT npy_int64
#define F_INT_NPY NPY_INT64

#if NPY_BITSOF_SHORT == 64
#define F_INT_PYFMT   "h"
#elif NPY_BITSOF_INT == 64
#define F_INT_PYFMT   "i"
#elif NPY_BITSOF_LONG == 64
#define F_INT_PYFMT   "l"
#elif NPY_BITSOF_LONGLONG == 64
#define F_INT_PYFMT   "L"
#else
#error No compatible 64-bit integer size. \
       Please contact NumPy maintainers and give detailed information about your \
       compiler and platform, or set NPY_USE_BLAS64_=0
#endif

#else

#define F_INT int
#define F_INT_NPY NPY_INT
#define F_INT_PYFMT   "i"

#endif

#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F
#else
#define F_FUNC(f,F) f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F##_
#else
#define F_FUNC(f,F) f##_
#endif
#endif

#define PYERR(errobj,message) {PyErr_SetString(errobj,message); goto fail;}
#define PYERR2(errobj,message) {PyErr_Print(); PyErr_SetString(errobj, message); goto fail;}
#define ISCONTIGUOUS(m) ((m)->flags & CONTIGUOUS)

#define MAX(n1,n2) ((n1) > (n2))?(n1):(n2);
#define MIN(n1,n2) ((n1) > (n2))?(n2):(n1);

struct ODR_info_ {
  PyObject* fcn;
  PyObject* fjacb;
  PyObject* fjacd;
  PyObject* pyBeta;
  PyObject* extra_args;
};

typedef struct ODR_info_ ODR_info;

static ODR_info odr_global;

static PyObject *odr_error=NULL;
static PyObject *odr_stop=NULL;

void fcn_callback(F_INT *n, F_INT *m, F_INT *np, F_INT *nq, F_INT *ldn, F_INT *ldm,
		  F_INT *ldnp, double *beta, double *xplusd, F_INT *ifixb,
		  F_INT *ifixx, F_INT *ldfix, F_INT *ideval, double *f,
		  double *fjacb, double *fjacd, F_INT *istop);

PyObject *gen_output(F_INT n, F_INT m, F_INT np, F_INT nq, F_INT ldwe, F_INT ld2we,
		     PyArrayObject *beta, PyArrayObject *work, PyArrayObject *iwork,
		     F_INT isodr, F_INT info, int full_output);

PyObject *odr(PyObject *self, PyObject *args, PyObject *kwds);

#define PyArray_CONTIGUOUS(m) (ISCONTIGUOUS(m) ? Py_INCREF(m), m : \
(PyArrayObject *)(PyArray_ContiguousFromObject((PyObject *)(m), \
(m)->descr->type_num, 0,0)))
#define D(dbg) printf("we're here: %i\n", dbg)
#define EXIST(name,obj) if (obj==NULL){printf("%s\n",name);}
