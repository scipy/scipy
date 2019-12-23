#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include "numpy/arrayobject.h"

#include "numpy/npy_3kcompat.h"

#ifdef HAVE_BLAS_ILP64

#define INT_TYPE npy_int64
#define INT_TYPE_NPY NPY_INT64

#if NPY_BITSOF_SHORT == 64
#define INT_TYPE_PYFMT   "h"
#define INT_TYPE_FMT "%h"
#elif NPY_BITSOF_INT == 64
#define INT_TYPE_PYFMT   "i"
#define INT_TYPE_FMT "%d"
#elif NPY_BITSOF_LONG == 64
#define INT_TYPE_PYFMT   "l"
#define INT_TYPE_FMT "%ld"
#elif NPY_BITSOF_LONGLONG == 64
#define INT_TYPE_PYFMT   "L"
#define INT_TYPE_FMT "%lld"
#else
#error No compatible 64-bit integer size. \
       Please contact NumPy maintainers and give detailed information about your \
       compiler and platform, or set NPY_USE_BLAS64_=0
#endif

#else

#define INT_TYPE int
#define INT_TYPE_NPY NPY_INT
#define INT_TYPE_FMT NPY_INT_FMT
#define INT_TYPE_PYFMT   "i"

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

void fcn_callback(INT_TYPE *n, INT_TYPE *m, INT_TYPE *np, INT_TYPE *nq, INT_TYPE *ldn, INT_TYPE *ldm,
		  INT_TYPE *ldnp, double *beta, double *xplusd, INT_TYPE *ifixb,
		  INT_TYPE *ifixx, INT_TYPE *ldfix, INT_TYPE *ideval, double *f,
		  double *fjacb, double *fjacd, INT_TYPE *istop);

PyObject *gen_output(INT_TYPE n, INT_TYPE m, INT_TYPE np, INT_TYPE nq, INT_TYPE ldwe, INT_TYPE ld2we,
		     PyArrayObject *beta, PyArrayObject *work, PyArrayObject *iwork,
		     INT_TYPE isodr, INT_TYPE info, int full_output);

PyObject *odr(PyObject *self, PyObject *args, PyObject *kwds);

#define PyArray_CONTIGUOUS(m) (ISCONTIGUOUS(m) ? Py_INCREF(m), m : \
(PyArrayObject *)(PyArray_ContiguousFromObject((PyObject *)(m), \
(m)->descr->type_num, 0,0)))
#define D(dbg) printf("we're here: %i\n", dbg)
#define EXIST(name,obj) if (obj==NULL){printf("%s\n",name);}
