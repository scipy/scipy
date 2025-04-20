/* -*- c -*- */

/*
 *****************************************************************************
 **                            INCLUDES                                     **
 *****************************************************************************
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "numpy/arrayobject.h"
#include "numpy/ufuncobject.h"
#include "numpy/npy_math.h"

//#include "npy_config.h"

#include "../_build_utils/src/npy_cblas.h"

#include <cstddef>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <type_traits>
#include <utility>

static const char* umath_linalg_version_string = "0.1.5";

/*
 ****************************************************************************
 *                        Debugging support                                 *
 ****************************************************************************
 */
#define _UMATH_LINALG_DEBUG 0

#define TRACE_TXT(...) do { fprintf (stderr, __VA_ARGS__); } while (0)
#define STACK_TRACE do {} while (0)
#define TRACE\
    do {                                        \
        fprintf (stderr,                        \
                 "%s:%d:%s\n",                  \
                 __FILE__,                      \
                 __LINE__,                      \
                 __FUNCTION__);                 \
        STACK_TRACE;                            \
    } while (0)

#if _UMATH_LINALG_DEBUG
#if defined HAVE_EXECINFO_H
#include <execinfo.h>
#elif defined HAVE_LIBUNWIND_H
#include <libunwind.h>
#endif
void
dbg_stack_trace()
{
    void *trace[32];
    size_t size;

    size = backtrace(trace, sizeof(trace)/sizeof(trace[0]));
    backtrace_symbols_fd(trace, size, 1);
}

#undef STACK_TRACE
#define STACK_TRACE do { dbg_stack_trace(); } while (0)
#endif

/*
 *****************************************************************************
 *                    BLAS/LAPACK calling macros                             *
 *****************************************************************************
 */

#define FNAME(x) BLAS_FUNC(x)

typedef CBLAS_INT         fortran_int;

typedef struct { float r, i; } f2c_complex;
typedef struct { double r, i; } f2c_doublecomplex;
/* typedef long int (*L_fp)(); */

typedef float             fortran_real;
typedef double            fortran_doublereal;
typedef f2c_complex       fortran_complex;
typedef f2c_doublecomplex fortran_doublecomplex;


extern "C" fortran_int
BLAS_FUNC(scopy)(fortran_int *n,
        float *sx, fortran_int *incx,
        float *sy, fortran_int *incy);
extern "C" fortran_int
BLAS_FUNC(dcopy)(fortran_int *n,
        double *sx, fortran_int *incx,
        double *sy, fortran_int *incy);
extern "C" fortran_int
BLAS_FUNC(ccopy)(fortran_int *n,
        f2c_complex *sx, fortran_int *incx,
        f2c_complex *sy, fortran_int *incy);
extern "C" fortran_int
BLAS_FUNC(zcopy)(fortran_int *n,
        f2c_doublecomplex *sx, fortran_int *incx,
        f2c_doublecomplex *sy, fortran_int *incy);

extern "C" fortran_int
FNAME(sgesv)(fortran_int *n, fortran_int *nrhs,
             float a[], fortran_int *lda,
             fortran_int ipiv[],
             float b[], fortran_int *ldb,
             fortran_int *info);
extern "C" fortran_int
FNAME(dgesv)(fortran_int *n, fortran_int *nrhs,
             double a[], fortran_int *lda,
             fortran_int ipiv[],
             double b[], fortran_int *ldb,
             fortran_int *info);
extern "C" fortran_int
FNAME(cgesv)(fortran_int *n, fortran_int *nrhs,
             f2c_complex a[], fortran_int *lda,
             fortran_int ipiv[],
             f2c_complex b[], fortran_int *ldb,
             fortran_int *info);
extern "C" fortran_int
FNAME(zgesv)(fortran_int *n, fortran_int *nrhs,
             f2c_doublecomplex a[], fortran_int *lda,
             fortran_int ipiv[],
             f2c_doublecomplex b[], fortran_int *ldb,
             fortran_int *info);

#define LAPACK_T(FUNC)                                          \
    TRACE_TXT("Calling LAPACK ( " # FUNC " )\n");               \
    FNAME(FUNC)

#define BLAS(FUNC)                              \
    FNAME(FUNC)

#define LAPACK(FUNC)                            \
    FNAME(FUNC)


/*
 *****************************************************************************
 **                      Some handy functions                               **
 *****************************************************************************
 */

static inline int
get_fp_invalid_and_clear(void)
{
    int status;
    status = npy_clear_floatstatus_barrier((char*)&status);
    return !!(status & NPY_FPE_INVALID);
  return 42;
}

static inline void
set_fp_invalid_or_clear(int error_occurred)
{
    if (error_occurred) {
        npy_set_floatstatus_invalid();
    }
    else {
        npy_clear_floatstatus_barrier((char*)&error_occurred);
    }
}



/*
 *****************************************************************************
 **                      Some handy constants                               **
 *****************************************************************************
 */

#define UMATH_LINALG_MODULE_NAME "_umath_linalg"

template<typename T>
struct numeric_limits;

template<>
struct numeric_limits<float> {
static constexpr float one = 1.0f;
static constexpr float zero = 0.0f;
static constexpr float minus_one = -1.0f;
static const float ninf;
static const float nan;
};
constexpr float numeric_limits<float>::one;
constexpr float numeric_limits<float>::zero;
constexpr float numeric_limits<float>::minus_one;
const float numeric_limits<float>::ninf = -NPY_INFINITYF;
const float numeric_limits<float>::nan = NPY_NANF;

template<>
struct numeric_limits<double> {
static constexpr double one = 1.0;
static constexpr double zero = 0.0;
static constexpr double minus_one = -1.0;
static const double ninf;
static const double nan;
};
constexpr double numeric_limits<double>::one;
constexpr double numeric_limits<double>::zero;
constexpr double numeric_limits<double>::minus_one;
const double numeric_limits<double>::ninf = -NPY_INFINITY;
const double numeric_limits<double>::nan = NPY_NAN;

template<>
struct numeric_limits<npy_cfloat> {
static constexpr npy_cfloat one = {1.0f};
static constexpr npy_cfloat zero = {0.0f};
static constexpr npy_cfloat minus_one = {-1.0f};
static const npy_cfloat ninf;
static const npy_cfloat nan;
};
constexpr npy_cfloat numeric_limits<npy_cfloat>::one;
constexpr npy_cfloat numeric_limits<npy_cfloat>::zero;
constexpr npy_cfloat numeric_limits<npy_cfloat>::minus_one;
const npy_cfloat numeric_limits<npy_cfloat>::ninf = {-NPY_INFINITYF};
const npy_cfloat numeric_limits<npy_cfloat>::nan = {NPY_NANF, NPY_NANF};

template<>
struct numeric_limits<f2c_complex> {
static constexpr f2c_complex one = {1.0f, 0.0f};
static constexpr f2c_complex zero = {0.0f, 0.0f};
static constexpr f2c_complex minus_one = {-1.0f, 0.0f};
static const f2c_complex ninf;
static const f2c_complex nan;
};
constexpr f2c_complex numeric_limits<f2c_complex>::one;
constexpr f2c_complex numeric_limits<f2c_complex>::zero;
constexpr f2c_complex numeric_limits<f2c_complex>::minus_one;
const f2c_complex numeric_limits<f2c_complex>::ninf = {-NPY_INFINITYF, 0.0f};
const f2c_complex numeric_limits<f2c_complex>::nan = {NPY_NANF, NPY_NANF};

template<>
struct numeric_limits<npy_cdouble> {
static constexpr npy_cdouble one = {1.0};
static constexpr npy_cdouble zero = {0.0};
static constexpr npy_cdouble minus_one = {-1.0};
static const npy_cdouble ninf;
static const npy_cdouble nan;
};
constexpr npy_cdouble numeric_limits<npy_cdouble>::one;
constexpr npy_cdouble numeric_limits<npy_cdouble>::zero;
constexpr npy_cdouble numeric_limits<npy_cdouble>::minus_one;
const npy_cdouble numeric_limits<npy_cdouble>::ninf = {-NPY_INFINITY};
const npy_cdouble numeric_limits<npy_cdouble>::nan = {NPY_NAN, NPY_NAN};

template<>
struct numeric_limits<f2c_doublecomplex> {
static constexpr f2c_doublecomplex one = {1.0};
static constexpr f2c_doublecomplex zero = {0.0};
static constexpr f2c_doublecomplex minus_one = {-1.0};
static const f2c_doublecomplex ninf;
static const f2c_doublecomplex nan;
};
constexpr f2c_doublecomplex numeric_limits<f2c_doublecomplex>::one;
constexpr f2c_doublecomplex numeric_limits<f2c_doublecomplex>::zero;
constexpr f2c_doublecomplex numeric_limits<f2c_doublecomplex>::minus_one;
const f2c_doublecomplex numeric_limits<f2c_doublecomplex>::ninf = {-NPY_INFINITY};
const f2c_doublecomplex numeric_limits<f2c_doublecomplex>::nan = {NPY_NAN, NPY_NAN};



/*
 *****************************************************************************
 **               Structs used for data rearrangement                       **
 *****************************************************************************
 */


/*
 * this struct contains information about how to linearize a matrix in a local
 * buffer so that it can be used by blas functions.  All strides are specified
 * in bytes and are converted to elements later in type specific functions.
 *
 * rows: number of rows in the matrix
 * columns: number of columns in the matrix
 * row_strides: the number bytes between consecutive rows.
 * column_strides: the number of bytes between consecutive columns.
 * output_lead_dim: BLAS/LAPACK-side leading dimension, in elements
 */
struct linearize_data
{
  npy_intp rows;
  npy_intp columns;
  npy_intp row_strides;
  npy_intp column_strides;
  npy_intp output_lead_dim;
};

static inline
linearize_data init_linearize_data_ex(npy_intp rows,
                       npy_intp columns,
                       npy_intp row_strides,
                       npy_intp column_strides,
                       npy_intp output_lead_dim)
{
    return {rows, columns, row_strides, column_strides, output_lead_dim};
}

static inline
linearize_data init_linearize_data(npy_intp rows,
                    npy_intp columns,
                    npy_intp row_strides,
                    npy_intp column_strides)
{
    return init_linearize_data_ex(
        rows, columns, row_strides, column_strides, columns);
}



/*
 *****************************************************************************
 **                            Basics                                       **
 *****************************************************************************
 */

static inline fortran_int
fortran_int_min(fortran_int x, fortran_int y) {
    return x < y ? x : y;
}

static inline fortran_int
fortran_int_max(fortran_int x, fortran_int y) {
    return x > y ? x : y;
}

#define INIT_OUTER_LOOP_1 \
    npy_intp dN = *dimensions++;\
    npy_intp N_;\
    npy_intp s0 = *steps++;

#define INIT_OUTER_LOOP_2 \
    INIT_OUTER_LOOP_1\
    npy_intp s1 = *steps++;

#define INIT_OUTER_LOOP_3 \
    INIT_OUTER_LOOP_2\
    npy_intp s2 = *steps++;

#define INIT_OUTER_LOOP_4 \
    INIT_OUTER_LOOP_3\
    npy_intp s3 = *steps++;

#define INIT_OUTER_LOOP_5 \
    INIT_OUTER_LOOP_4\
    npy_intp s4 = *steps++;

#define INIT_OUTER_LOOP_6  \
    INIT_OUTER_LOOP_5\
    npy_intp s5 = *steps++;

#define INIT_OUTER_LOOP_7  \
    INIT_OUTER_LOOP_6\
    npy_intp s6 = *steps++;

#define BEGIN_OUTER_LOOP_2 \
    for (N_ = 0;\
         N_ < dN;\
         N_++, args[0] += s0,\
             args[1] += s1) {

#define BEGIN_OUTER_LOOP_3 \
    for (N_ = 0;\
         N_ < dN;\
         N_++, args[0] += s0,\
             args[1] += s1,\
             args[2] += s2) {

#define BEGIN_OUTER_LOOP_4 \
    for (N_ = 0;\
         N_ < dN;\
         N_++, args[0] += s0,\
             args[1] += s1,\
             args[2] += s2,\
             args[3] += s3) {

#define BEGIN_OUTER_LOOP_5 \
    for (N_ = 0;\
         N_ < dN;\
         N_++, args[0] += s0,\
             args[1] += s1,\
             args[2] += s2,\
             args[3] += s3,\
             args[4] += s4) {

#define BEGIN_OUTER_LOOP_6 \
    for (N_ = 0;\
         N_ < dN;\
         N_++, args[0] += s0,\
             args[1] += s1,\
             args[2] += s2,\
             args[3] += s3,\
             args[4] += s4,\
             args[5] += s5) {

#define BEGIN_OUTER_LOOP_7 \
    for (N_ = 0;\
         N_ < dN;\
         N_++, args[0] += s0,\
             args[1] += s1,\
             args[2] += s2,\
             args[3] += s3,\
             args[4] += s4,\
             args[5] += s5,\
             args[6] += s6) {

#define END_OUTER_LOOP  }

static inline void
update_pointers(npy_uint8** bases, ptrdiff_t* offsets, size_t count)
{
    size_t i;
    for (i = 0; i < count; ++i) {
        bases[i] += offsets[i];
    }
}


/*
 *****************************************************************************
 **                             DISPATCHER FUNCS                            **
 *****************************************************************************
 */

static fortran_int copy(fortran_int *n,
        float *sx, fortran_int *incx,
        float *sy, fortran_int *incy) { return FNAME(scopy)(n, sx, incx,
            sy, incy);
}
static fortran_int copy(fortran_int *n,
        double *sx, fortran_int *incx,
        double *sy, fortran_int *incy) { return FNAME(dcopy)(n, sx, incx,
            sy, incy);
}
static fortran_int copy(fortran_int *n,
        f2c_complex *sx, fortran_int *incx,
        f2c_complex *sy, fortran_int *incy) { return FNAME(ccopy)(n, sx, incx,
            sy, incy);
}
static fortran_int copy(fortran_int *n,
        f2c_doublecomplex *sx, fortran_int *incx,
        f2c_doublecomplex *sy, fortran_int *incy) { return FNAME(zcopy)(n, sx, incx,
            sy, incy);
}


/*
 *****************************************************************************
 **                             HELPER FUNCS                                **
 *****************************************************************************
 */
template<typename T>
struct fortran_type {
using type = T;
};

template<> struct fortran_type<npy_cfloat> { using type = f2c_complex;};
template<> struct fortran_type<npy_cdouble> { using type = f2c_doublecomplex;};
template<typename T>
using fortran_type_t = typename fortran_type<T>::type;

template<typename T>
struct basetype {
using type = T;
};
template<> struct basetype<npy_cfloat> { using type = npy_float;};
template<> struct basetype<npy_cdouble> { using type = npy_double;};
template<> struct basetype<f2c_complex> { using type = fortran_real;};
template<> struct basetype<f2c_doublecomplex> { using type = fortran_doublereal;};
template<typename T>
using basetype_t = typename basetype<T>::type;

struct scalar_trait {};
struct complex_trait {};
template<typename typ>
using dispatch_scalar = typename std::conditional<sizeof(basetype_t<typ>) == sizeof(typ), scalar_trait, complex_trait>::type;


             /* rearranging of 2D matrices using blas */

template<typename typ>
static inline void *
linearize_matrix(typ *dst,
                        typ *src,
                        const linearize_data* data)
{
    using ftyp = fortran_type_t<typ>;
    if (dst) {
        int i, j;
        typ* rv = dst;
        fortran_int columns = (fortran_int)data->columns;
        fortran_int column_strides =
            (fortran_int)(data->column_strides/sizeof(typ));
        fortran_int one = 1;
        for (i = 0; i < data->rows; i++) {
            if (column_strides > 0) {
                copy(&columns,
                              (ftyp*)src, &column_strides,
                              (ftyp*)dst, &one);
            }
            else if (column_strides < 0) {
                copy(&columns,
                              ((ftyp*)src + (columns-1)*column_strides),
                              &column_strides,
                              (ftyp*)dst, &one);
            }
            else {
                /*
                 * Zero stride has undefined behavior in some BLAS
                 * implementations (e.g. OSX Accelerate), so do it
                 * manually
                 */
                for (j = 0; j < columns; ++j) {
                    memcpy(dst + j, src, sizeof(typ));
                }
            }
            src += data->row_strides/sizeof(typ);
            dst += data->output_lead_dim;
        }
        return rv;
    } else {
        return src;
    }
}

template<typename typ>
static inline void *
delinearize_matrix(typ *dst,
                          typ *src,
                          const linearize_data* data)
{
using ftyp = fortran_type_t<typ>;

    if (src) {
        int i;
        typ *rv = src;
        fortran_int columns = (fortran_int)data->columns;
        fortran_int column_strides =
            (fortran_int)(data->column_strides/sizeof(typ));
        fortran_int one = 1;
        for (i = 0; i < data->rows; i++) {
            if (column_strides > 0) {
                copy(&columns,
                              (ftyp*)src, &one,
                              (ftyp*)dst, &column_strides);
            }
            else if (column_strides < 0) {
                copy(&columns,
                              (ftyp*)src, &one,
                              ((ftyp*)dst + (columns-1)*column_strides),
                              &column_strides);
            }
            else {
                /*
                 * Zero stride has undefined behavior in some BLAS
                 * implementations (e.g. OSX Accelerate), so do it
                 * manually
                 */
                if (columns > 0) {
                    memcpy(dst,
                           src + (columns-1),
                           sizeof(typ));
                }
            }
            src += data->output_lead_dim;
            dst += data->row_strides/sizeof(typ);
        }

        return rv;
    } else {
        return src;
    }
}

template<typename typ>
static inline void
nan_matrix(typ *dst, const linearize_data* data)
{
    int i, j;
    for (i = 0; i < data->rows; i++) {
        typ *cp = dst;
        ptrdiff_t cs = data->column_strides/sizeof(typ);
        for (j = 0; j < data->columns; ++j) {
            *cp = numeric_limits<typ>::nan;
            cp += cs;
        }
        dst += data->row_strides/sizeof(typ);
    }
}

template<typename typ>
static inline void
zero_matrix(typ *dst, const linearize_data* data)
{
    int i, j;
    for (i = 0; i < data->rows; i++) {
        typ *cp = dst;
        ptrdiff_t cs = data->column_strides/sizeof(typ);
        for (j = 0; j < data->columns; ++j) {
            *cp = numeric_limits<typ>::zero;
            cp += cs;
        }
        dst += data->row_strides/sizeof(typ);
    }
}

               /* identity square matrix generation */
template<typename typ>
static inline void
identity_matrix(typ *matrix, size_t n)
{
    size_t i;
    /* in IEEE floating point, zeroes are represented as bitwise 0 */
    memset((void *)matrix, 0, n*n*sizeof(typ));

    for (i = 0; i < n; ++i)
    {
        *matrix = numeric_limits<typ>::one;
        matrix += n+1;
    }
}



/* -------------------------------------------------------------------------- */
                  /* Solve family (includes inv) */

template<typename typ>
struct GESV_PARAMS_t
{
    typ *A; /* A is (N, N) of base type */
    typ *B; /* B is (N, NRHS) of base type */
    fortran_int * IPIV; /* IPIV is (N) */

    fortran_int N;
    fortran_int NRHS;
    fortran_int LDA;
    fortran_int LDB;
};

static inline fortran_int
call_gesv(GESV_PARAMS_t<fortran_real> *params)
{
    fortran_int rv;
    LAPACK(sgesv)(&params->N, &params->NRHS,
                          params->A, &params->LDA,
                          params->IPIV,
                          params->B, &params->LDB,
                          &rv);
    return rv;
}

static inline fortran_int
call_gesv(GESV_PARAMS_t<fortran_doublereal> *params)
{
    fortran_int rv;
    LAPACK(dgesv)(&params->N, &params->NRHS,
                          params->A, &params->LDA,
                          params->IPIV,
                          params->B, &params->LDB,
                          &rv);
    return rv;
}

static inline fortran_int
call_gesv(GESV_PARAMS_t<fortran_complex> *params)
{
    fortran_int rv;
    LAPACK(cgesv)(&params->N, &params->NRHS,
                          params->A, &params->LDA,
                          params->IPIV,
                          params->B, &params->LDB,
                          &rv);
    return rv;
}

static inline fortran_int
call_gesv(GESV_PARAMS_t<fortran_doublecomplex> *params)
{
    fortran_int rv;
    LAPACK(zgesv)(&params->N, &params->NRHS,
                          params->A, &params->LDA,
                          params->IPIV,
                          params->B, &params->LDB,
                          &rv);
    return rv;
}


/*
 * Initialize the parameters to use in for the lapack function _heev
 * Handles buffer allocation
 */
template<typename ftyp>
static inline int
init_gesv(GESV_PARAMS_t<ftyp> *params, fortran_int N, fortran_int NRHS)
{
    npy_uint8 *mem_buff = NULL;
    npy_uint8 *a, *b, *ipiv;
    size_t safe_N = N;
    size_t safe_NRHS = NRHS;
    fortran_int ld = fortran_int_max(N, 1);
    mem_buff = (npy_uint8 *)malloc(safe_N * safe_N * sizeof(ftyp) +
                      safe_N * safe_NRHS*sizeof(ftyp) +
                      safe_N * sizeof(fortran_int));
    if (!mem_buff) {
        goto error;
    }
    a = mem_buff;
    b = a + safe_N * safe_N * sizeof(ftyp);
    ipiv = b + safe_N * safe_NRHS * sizeof(ftyp);

    params->A = (ftyp*)a;
    params->B = (ftyp*)b;
    params->IPIV = (fortran_int*)ipiv;
    params->N = N;
    params->NRHS = NRHS;
    params->LDA = ld;
    params->LDB = ld;

    return 1;
 error:
    free(mem_buff);
    memset(params, 0, sizeof(*params));

    return 0;
}


template<typename ftyp>
static inline void
release_gesv(GESV_PARAMS_t<ftyp> *params)
{
    /* memory block base is in A */
    free(params->A);
    memset(params, 0, sizeof(*params));
}


template<typename typ>
static void
inv(char **args, npy_intp const *dimensions, npy_intp const *steps,
           void *NPY_UNUSED(func))
{
using ftyp = fortran_type_t<typ>;
    GESV_PARAMS_t<ftyp> params;
    fortran_int n;
    int error_occurred = get_fp_invalid_and_clear();
    INIT_OUTER_LOOP_2

    n = (fortran_int)dimensions[0];
    if (init_gesv(&params, n, n)) {
        linearize_data a_in = init_linearize_data(n, n, steps[1], steps[0]);
        linearize_data r_out = init_linearize_data(n, n, steps[3], steps[2]);

        BEGIN_OUTER_LOOP_2
            int not_ok;
            linearize_matrix((typ*)params.A, (typ*)args[0], &a_in);
            identity_matrix((typ*)params.B, n);
            not_ok = call_gesv(&params);
            if (!not_ok) {
                delinearize_matrix((typ*)args[1], (typ*)params.B, &r_out);
            } else {
                error_occurred = 1;
                nan_matrix((typ*)args[1], &r_out);
            }
        END_OUTER_LOOP

        release_gesv(&params);
    }

    set_fp_invalid_or_clear(error_occurred);
}


/* -------------------------------------------------------------------------- */
              /* gufunc registration  */

static void *array_of_nulls[] = {
    (void *)NULL,
    (void *)NULL,
    (void *)NULL,
    (void *)NULL,

    (void *)NULL,
    (void *)NULL,
    (void *)NULL,
    (void *)NULL,

    (void *)NULL,
    (void *)NULL,
    (void *)NULL,
    (void *)NULL,

    (void *)NULL,
    (void *)NULL,
    (void *)NULL,
    (void *)NULL
};

#define FUNC_ARRAY_NAME(NAME) NAME ## _funcs

#define GUFUNC_FUNC_ARRAY_REAL(NAME)                    \
    static PyUFuncGenericFunction                       \
    FUNC_ARRAY_NAME(NAME)[] = {                         \
        FLOAT_ ## NAME,                                 \
        DOUBLE_ ## NAME                                 \
    }

#define GUFUNC_FUNC_ARRAY_REAL_COMPLEX(NAME)            \
    static PyUFuncGenericFunction                       \
    FUNC_ARRAY_NAME(NAME)[] = {                         \
        FLOAT_ ## NAME,                                 \
        DOUBLE_ ## NAME,                                \
        CFLOAT_ ## NAME,                                \
        CDOUBLE_ ## NAME                                \
    }
#define GUFUNC_FUNC_ARRAY_REAL_COMPLEX_(NAME)            \
    static PyUFuncGenericFunction                       \
    FUNC_ARRAY_NAME(NAME)[] = {                         \
        NAME<npy_float, npy_float>,                                 \
        NAME<npy_double, npy_double>,                                \
        NAME<npy_cfloat, npy_float>,                                \
        NAME<npy_cdouble, npy_double>                                \
    }
#define GUFUNC_FUNC_ARRAY_REAL_COMPLEX__(NAME)            \
    static PyUFuncGenericFunction                       \
    FUNC_ARRAY_NAME(NAME)[] = {                         \
        NAME<npy_float>,                                 \
        NAME<npy_double>,                                \
        NAME<npy_cfloat>,                                \
        NAME<npy_cdouble>                                \
    }

GUFUNC_FUNC_ARRAY_REAL_COMPLEX__(inv);


static const char equal_2_types[] = {
    NPY_FLOAT, NPY_FLOAT,
    NPY_DOUBLE, NPY_DOUBLE,
    NPY_CFLOAT, NPY_CFLOAT,
    NPY_CDOUBLE, NPY_CDOUBLE
};

static const char equal_3_types[] = {
    NPY_FLOAT, NPY_FLOAT, NPY_FLOAT,
    NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
    NPY_CFLOAT, NPY_CFLOAT, NPY_CFLOAT,
    NPY_CDOUBLE, NPY_CDOUBLE, NPY_CDOUBLE
};



typedef struct gufunc_descriptor_struct {
    const char *name;
    const char *signature;
    const char *doc;
    int ntypes;
    int nin;
    int nout;
    PyUFuncGenericFunction *funcs;
    const char *types;
    // XXX: process_core_dims_func does not exist in numpy 1.26 ?
    //PyUFunc_ProcessCoreDimsFunc *process_core_dims_func;
    void *process_core_dims_func;
} GUFUNC_DESCRIPTOR_t;

GUFUNC_DESCRIPTOR_t gufunc_descriptors [] = {
    {
        "inv",
        "(m, m)->(m, m)",
        "compute the inverse of the last two dimensions and broadcast"\
        " to the rest. \n"\
        "Results in the inverse matrices. \n"\
        "    \"(m,m)->(m,m)\" \n",
        4, 1, 1,
        FUNC_ARRAY_NAME(inv),
        equal_2_types,
        nullptr
    }
};

static int
addUfuncs(PyObject *dictionary) {
    PyUFuncObject *f;
    int i;
    const int gufunc_count = sizeof(gufunc_descriptors)/
        sizeof(gufunc_descriptors[0]);
    for (i = 0; i < gufunc_count; i++) {
        GUFUNC_DESCRIPTOR_t* d = &gufunc_descriptors[i];
        f = (PyUFuncObject *) PyUFunc_FromFuncAndDataAndSignature(
                                                d->funcs,
                                                array_of_nulls,
                                                (char*)d->types, // XXX: char* cast for np 1.26
                                                d->ntypes,
                                                d->nin,
                                                d->nout,
                                                PyUFunc_None,
                                                d->name,
                                                d->doc,
                                                0,
                                                d->signature);
        if (f == NULL) {
            return -1;
        }
        // f->process_core_dims_func = d->process_core_dims_func; // XXX numpy 1.26 vs 2.x
#if _UMATH_LINALG_DEBUG
        dump_ufunc_object((PyUFuncObject*) f);
#endif
        int ret = PyDict_SetItemString(dictionary, d->name, (PyObject *)f);
        Py_DECREF(f);
        if (ret < 0) {
            return -1;
        }
    }
    return 0;
}



/* -------------------------------------------------------------------------- */
                  /* Module initialization and state  */

static PyMethodDef UMath_LinAlgMethods[] = {
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        UMATH_LINALG_MODULE_NAME,
        NULL,
        -1,
        UMath_LinAlgMethods,
        NULL,
        NULL,
        NULL,
        NULL
};

PyMODINIT_FUNC PyInit__umath_linalg(void)
{
    PyObject *m;
    PyObject *d;
    PyObject *version;

    m = PyModule_Create(&moduledef);
    if (m == NULL) {
        return NULL;
    }

    import_array();
    import_ufunc();

    d = PyModule_GetDict(m);
    if (d == NULL) {
        return NULL;
    }

    version = PyUnicode_FromString(umath_linalg_version_string);
    if (version == NULL) {
        return NULL;
    }
    int ret = PyDict_SetItemString(d, "__version__", version);
    Py_DECREF(version);
    if (ret < 0) {
        return NULL;
    }

    /* Load the ufunc operators into the module's namespace */
    if (addUfuncs(d) < 0) {
        return NULL;
    }

// XXX: ILP64 macros? GIL_DISABLED?

#ifdef HAVE_BLAS_ILP64
    PyDict_SetItemString(d, "_ilp64", Py_True);
#else
    PyDict_SetItemString(d, "_ilp64", Py_False);
#endif

#if Py_GIL_DISABLED
    // signal this module supports running with the GIL disabled
    PyUnstable_Module_SetGIL(m, Py_MOD_GIL_NOT_USED);
#endif

    return m;
}
