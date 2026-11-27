#ifndef GLB_H_GUARD
#define GLB_H_GUARD

#include "HConfig.h"

// #ifndef CUPDLP_CPU
// #include <cublas_v2.h>         // cublas
// #include <cuda_runtime_api.h>  // cudaMalloc, cudaMemcpy, etc.
// #include <cusparse.h>          // cusparseSpMV
// #endif

#ifdef __cplusplus

extern "C" {
#endif

#include <math.h>

// return code
#define RETCODE_OK (0)
#define RETCODE_FAILED (1)
#define BOOL (1)
// #define DLONG

#ifndef cupdlp
#define cupdlp(x) cupdlp_##x
#endif

/* cupdlp VERSION NUMBER ----------------------------------------------    */
#define cupdlp_VERSION \
  ("1.0.0") /* string literals automatically null-terminated */

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define cupdlp_printf mexPrintf
#define _cupdlp_free mxFree
#define _cupdlp_malloc mxMalloc
#define _cupdlp_calloc mxCalloc
#define _cupdlp_realloc mxRealloc
#elif defined PYTHON
#include <Python.h>
#include <stdlib.h>
#define cupdlp_printf(...)                           \
  {                                                  \
    PyGILState_STATE gilstate = PyGILState_Ensure(); \
    PySys_WriteStdout(__VA_ARGS__);                  \
    PyGILState_Release(gilstate);                    \
  }
#define _cupdlp_free free
#define _cupdlp_malloc malloc
#define _cupdlp_calloc calloc
#define _cupdlp_realloc realloc
#else

#include <stdio.h>
#include <stdlib.h>

#define cupdlp_printf printf
#define cupdlp_snprintf snprintf
#define _cupdlp_free free
#define _cupdlp_malloc malloc
#define _cupdlp_calloc calloc
#define _cupdlp_realloc realloc
#endif

// for cuda
#ifdef CUPDLP_CPU

#define CUPDLP_COPY_VEC(dst, src, type, size) \
  memcpy(dst, src, sizeof(type) * (size))

  //CUPDLP_INIT_VEC is not used
  //
  //#define CUPDLP_INIT_VEC(var, size)				\
//  {                                                             \
//    (var) = (typeof(var))malloc((size) * sizeof(typeof(*var))); \
//    if ((var) == cupdlp_NULL) {                                 \
//      retcode = RETCODE_FAILED;                                 \
//      goto exit_cleanup;                                        \
//    }                                                           \
//  }
  //#define CUPDLP_INIT_ZERO_VEC(var, size)		     \
//  {                                                          \
//    (var) = (typeof(var))calloc(size, sizeof(typeof(*var))); \
//    if ((var) == cupdlp_NULL) {                              \
//      retcode = RETCODE_FAILED;                              \
//      goto exit_cleanup;                                     \
//    }                                                        \
//  }
#define cupdlp_init_zero_vec_double(var, size)                      \
  {                                                          \
    (var) = (double*)calloc(size, sizeof(double));	     \
    if ((var) == cupdlp_NULL) {                              \
      retcode = RETCODE_FAILED;                              \
      goto exit_cleanup;                                     \
    }                                                        \
  }
#define cupdlp_init_zero_vec_int(var, size)                      \
  {                                                          \
    (var) = (int*)calloc(size, sizeof(int));	     \
    if ((var) == cupdlp_NULL) {                              \
      retcode = RETCODE_FAILED;                              \
      goto exit_cleanup;                                     \
    }                                                        \
  }
#define CUPDLP_FREE_VEC(x) \
  {                        \
    _cupdlp_free(x);       \
    x = cupdlp_NULL;       \
  }
#define CUPDLP_ZERO_VEC(var, type, size) memset(var, 0, sizeof(type) * (size))

#endif

#define cupdlp_free(x) \
  {                    \
    _cupdlp_free(x);   \
    x = cupdlp_NULL;   \
  }
#define cupdlp_malloc(x) _cupdlp_malloc(x)
#define cupdlp_calloc(x, y) _cupdlp_calloc(x, y)
#define cupdlp_realloc(x, y) _cupdlp_realloc(x, y)
#define cupdlp_zero(var, type, size) memset(var, 0, sizeof(type) * (size))
#define cupdlp_copy(dst, src, type, size) \
  memcpy(dst, src, sizeof(type) * (size))
//#define CUPDLP_INIT(var, size)				\
//  {                                                             \
//    (var) = (typeof(var))malloc((size) * sizeof(typeof(*var))); \
//    if ((var) == cupdlp_NULL) {                                 \
//      retcode = RETCODE_FAILED;                                 \
//      goto exit_cleanup;                                        \
//    }                                                           \
//  }
#define CUPDLP_INIT_DOUBLE(var, size)                                  \
  {                                                             \
    (var) = (double*)malloc((size) * sizeof(double)); \
    if ((var) == cupdlp_NULL) {                                 \
      retcode = RETCODE_FAILED;                                 \
      goto exit_cleanup;                                        \
    }                                                           \
  }
#define CUPDLP_INIT_CUPDLP_VEC(var, size)                                  \
  {                                                             \
    (var) = (CUPDLPvec*)malloc((size) * sizeof(CUPDLPvec)); \
    if ((var) == cupdlp_NULL) {                                 \
      retcode = RETCODE_FAILED;                                 \
      goto exit_cleanup;                                        \
    }                                                           \
  }
#define CUPDLP_INIT_DENSE_MATRIX(var, size)                                  \
  {                                                             \
    (var) = (CUPDLPdense*)malloc((size) * sizeof(CUPDLPdense)); \
    if ((var) == cupdlp_NULL) {                                 \
      retcode = RETCODE_FAILED;                                 \
      goto exit_cleanup;                                        \
    }                                                           \
  }
#define CUPDLP_INIT_CSR_MATRIX(var, size)                                  \
  {                                                             \
    (var) = (CUPDLPcsr*)malloc((size) * sizeof(CUPDLPcsr)); \
    if ((var) == cupdlp_NULL) {                                 \
      retcode = RETCODE_FAILED;                                 \
      goto exit_cleanup;                                        \
    }                                                           \
  }
#define CUPDLP_INIT_CSC_MATRIX(var, size)                                  \
  {                                                             \
    (var) = (CUPDLPcsc*)malloc((size) * sizeof(CUPDLPcsc));	\
    if ((var) == cupdlp_NULL) {                                 \
      retcode = RETCODE_FAILED;                                 \
      goto exit_cleanup;                                        \
    }                                                           \
  }
#define CUPDLP_INIT_SETTINGS(var, size)                                  \
  {                                                             \
    (var) = (CUPDLPsettings*)malloc((size) * sizeof(CUPDLPsettings));	\
    if ((var) == cupdlp_NULL) {                                 \
      retcode = RETCODE_FAILED;                                 \
      goto exit_cleanup;                                        \
    }                                                           \
  }
#define CUPDLP_INIT_RESOBJ(var, size)                                  \
  {                                                             \
    (var) = (CUPDLPresobj*)malloc((size) * sizeof(CUPDLPresobj));	\
    if ((var) == cupdlp_NULL) {                                 \
      retcode = RETCODE_FAILED;                                 \
      goto exit_cleanup;                                        \
    }                                                           \
  }
#define CUPDLP_INIT_ITERATES(var, size)                                  \
  {                                                             \
    (var) = (CUPDLPiterates*)malloc((size) * sizeof(CUPDLPiterates));	\
    if ((var) == cupdlp_NULL) {                                 \
      retcode = RETCODE_FAILED;                                 \
      goto exit_cleanup;                                        \
    }                                                           \
  }
#define CUPDLP_INIT_STEPSIZE(var, size)                                  \
  {                                                             \
    (var) = (CUPDLPstepsize*)malloc((size) * sizeof(CUPDLPstepsize));	\
    if ((var) == cupdlp_NULL) {                                 \
      retcode = RETCODE_FAILED;                                 \
      goto exit_cleanup;                                        \
    }                                                           \
  }
#define CUPDLP_INIT_TIMERS(var, size)                                  \
  {                                                             \
    (var) = (CUPDLPtimers*)malloc((size) * sizeof(CUPDLPtimers));	\
    if ((var) == cupdlp_NULL) {                                 \
      retcode = RETCODE_FAILED;                                 \
      goto exit_cleanup;                                        \
    }                                                           \
  }
  //#define CUPDLP_INIT_ZERO(var, size)			     \
//  {                                                          \
//    (var) = (typeof(var))calloc(size, sizeof(typeof(*var))); \
//    if ((var) == cupdlp_NULL) {                              \
//      retcode = RETCODE_FAILED;                              \
//      goto exit_cleanup;                                     \
//    }                                                        \
//  }
#define CUPDLP_INIT_ZERO_DOUBLE(var, size)                          \
  {                                                          \
    (var) = (double*)calloc(size, sizeof(double)); \
    if ((var) == cupdlp_NULL) {                              \
      retcode = RETCODE_FAILED;                              \
      goto exit_cleanup;                                     \
    }                                                        \
  }
#define CUPDLP_INIT_ZERO_CUPDLP_WORK(var, size)                          \
  {                                                          \
    (var) = (CUPDLPwork*)calloc(size, sizeof(CUPDLPwork)); \
    if ((var) == cupdlp_NULL) {                              \
      retcode = RETCODE_FAILED;                              \
      goto exit_cleanup;                                     \
    }                                                        \
  }
#define CUPDLP_FREE(var) cupdlp_free(var)
#define CUPDLP_CALL(func)       \
  {                             \
    if ((func) != RETCODE_OK) { \
      retcode = RETCODE_FAILED; \
      goto exit_cleanup;        \
    }                           \
  }

#ifndef SFLOAT
#ifdef DLONG
typedef long long cupdlp_int;
#else
typedef int cupdlp_int;
#endif
typedef double cupdlp_float;
#ifndef NAN
#define NAN ((cupdlp_float)0x7ff8000000000000)
#endif
#ifndef INFINITY
#define INFINITY NAN
#endif
#else
typedef float cupdlp_float;
#ifndef NAN
#define NAN ((float)0x7fc00000)
#endif
#ifndef INFINITY
#define INFINITY NAN
#endif
#endif

#ifdef BOOL

#include <stdbool.h>

typedef bool cupdlp_bool;
#endif

#define cupdlp_NULL 0

#ifndef MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#endif

#ifndef ABS
#define ABS(x) (((x) < 0) ? -(x) : (x))
#endif

#ifndef POWF
#ifdef SFLOAT
#define POWF powf
#else
#define POWF pow
#endif
#endif

#ifndef SQRTF
#ifdef SFLOAT
#define SQRTF sqrtf
#else
#define SQRTF sqrt
#endif
#endif

#if EXTRA_VERBOSE > 1
#if (defined _WIN32 || defined _WIN64 || defined _WINDLL)
#define __func__ __FUNCTION__
#endif
#define DEBUG_FUNC                                                     \
  cupdlp_printf("IN function: %s, time: %4f ms, file: %s, line: %i\n", \
                __func__, cupdlp(tocq)(&global_timer), __FILE__, __LINE__);
#define RETURN
cupdlp_printf("EXIT function: %s, time: %4f ms, file: %s, line: %i\n", __func__,
              cupdlp(tocq)(&global_timer), __FILE__, __LINE__);
return
#else
#define DEBUG_FUNC
#define RETURN return
#endif

#define EPS_TOL (1E-18)
#define EPS (1E-8)  // for condition number in subnp
#define SAFEDIV_POS(X, Y) ((Y) < EPS_TOL ? ((X) / EPS_TOL) : (X) / (Y))

#define CONVERGED_INTERVAL (1)
#define INDETERMINATE_TOL (1e-9)

#define OUR_DBL_MAX 1E+20

#ifndef CUPDLP_ASSERT_H

#include <assert.h>

#define CUPDLP_ASSERT assert
#endif

#ifdef __cplusplus
}
#endif
#endif
