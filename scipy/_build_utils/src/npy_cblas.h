/*
 * This header provides numpy a consistent interface to CBLAS code. It is needed
 * because not all providers of cblas provide cblas.h. For instance, MKL provides
 * mkl_cblas.h and also typedefs the CBLAS_XXX enums.
 */
#ifndef _NPY_CBLAS_H_
#define _NPY_CBLAS_H_

#include <numpy/npy_common.h>
#include <stddef.h>

/* Allow the use in C++ code.  */
#ifdef __cplusplus
extern "C"
{
#endif

/*
 * Enumerated and derived types
 */
enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
enum CBLAS_UPLO {CblasUpper=121, CblasLower=122};
enum CBLAS_DIAG {CblasNonUnit=131, CblasUnit=132};
enum CBLAS_SIDE {CblasLeft=141, CblasRight=142};

#define CBLAS_INDEX size_t  /* this may vary between platforms */

#include "scipy_blas_defines.h"


#define BLASNAME(name) CBLAS_FUNC(name)
#define BLASINT CBLAS_INT

#include "npy_cblas_base.h"

#undef BLASINT
#undef BLASNAME


/*
 * Convert NumPy stride to BLAS stride. Returns 0 if conversion cannot be done
 * (BLAS won't handle negative or zero strides the way we want).
 */
static inline CBLAS_INT
blas_stride(npy_intp stride, unsigned itemsize)
{
    /*
     * Should probably check pointer alignment also, but this may cause
     * problems if we require complex to be 16 byte aligned.
     */
    if (stride > 0 && (stride % itemsize) == 0) {
        stride /= itemsize;
        if (stride <= CBLAS_INT_MAX) {
            return stride;
        }
    }
    return 0;
}

/*
 * Define a chunksize for CBLAS.
 *
 * The chunksize is the greatest power of two less than CBLAS_INT_MAX.
 */
#if NPY_MAX_INTP > CBLAS_INT_MAX
# define NPY_CBLAS_CHUNK  (CBLAS_INT_MAX / 2 + 1)
#else
# define NPY_CBLAS_CHUNK  NPY_MAX_INTP
#endif


#ifdef __cplusplus
}
#endif

#endif  /* _NPY_CBLAS_H_ */
