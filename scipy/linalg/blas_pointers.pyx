"""
Function pointers to BLAS functions.

BLAS Functions
==============

- caxpy
- ccopy
- cdotc
- cdotu
- cgemm
- cgemv
- cgerc
- cgeru
- chemm
- chemv
- cher
- cher2
- cher2k
- cherk
- crotg
- cscal
- csrot
- csscal
- cswap
- csymm
- csyr
- csyr2k
- csyrk
- ctrmm
- ctrmv
- dasum
- daxpy
- dcopy
- ddot
- dgemm
- dgemv
- dger
- dnrm2
- drot
- drotg
- drotm
- drotmg
- dscal
- dswap
- dsymm
- dsymv
- dsyr
- dsyr2
- dsyr2k
- dsyrk
- dtrmm
- dtrmv
- dzasum
- dznrm2
- icamax
- idamax
- isamax
- izamax
- sasum
- saxpy
- scasum
- scnrm2
- scopy
- sdot
- sgemm
- sgemv
- sger
- snrm2
- srot
- srotg
- srotm
- srotmg
- sscal
- sswap
- ssymm
- ssymv
- ssyr
- ssyr2
- ssyr2k
- ssyrk
- strmm
- strmv
- zaxpy
- zcopy
- zdotc
- zdotu
- zdrot
- zdscal
- zgemm
- zgemv
- zgerc
- zgeru
- zhemm
- zhemv
- zher
- zher2
- zher2k
- zherk
- zrotg
- zscal
- zswap
- zsymm
- zsyr
- zsyr2k
- zsyrk
- ztrmm
- ztrmv

"""

from . import  blas
from . import _blas_subroutine_wrappers as wblas

cdef extern from "f2pyptr.h":
    void *f2py_ptr(object) except NULL

# Use the subroutine wrappers for the
# functions with specific return values.
ctypedef int wcdotc_t(c *out, int *n, c *cx, int *incx, c *cy, int *incy) nogil
ctypedef int wcdotu_t(c *out, int *n, c *cx, int *incx, c *cy, int *incy) nogil
ctypedef int wdasum_t(d *out, int *n, d *dx, int *incx) nogil
ctypedef int wddot_t(d *out, int *n, d *dx, int *incx, d *dy, int *incy) nogil
ctypedef int wdnrm2_t(d *out, int *n, d *x, int *incx) nogil
ctypedef int wdzasum_t(d *out, int *n, z *zx, int *incx) nogil
ctypedef int wdznrm2_t(d *out, int *n, z *x, int *incx) nogil
ctypedef int wicamax_t(int *out, int *n, c *cx, int *incx) nogil
ctypedef int widamax_t(int *out, int *n, d *dx, int *incx) nogil
ctypedef int wisamax_t(int *out, int *n, s *sx, int *incx) nogil
ctypedef int wizamax_t(int *out, int *n, z *zx, int *incx) nogil
ctypedef int wsasum_t(s *out, int *n, s *sx, int *incx) nogil
ctypedef int wscasum_t(s *out, int *n, c *cx, int *incx) nogil
ctypedef int wscnrm2_t(s *out, int *n, c *x, int *incx) nogil
ctypedef int wsdot_t(s *out, int *n, s *sx, int *incx, s *sy, int *incy) nogil
ctypedef int wsnrm2_t(s *out, int *n, s *x, int *incx) nogil
ctypedef int wzdotc_t(z *out, int *n, z *zx, int *incx, z *zy, int *incy) nogil
ctypedef int wzdotu_t(z *out, int *n, z *zx, int *incx, z *zy, int *incy) nogil

# Get function pointers to the wrapper subroutines.
cdef:
    wcdotc_t *wcdotc = <wcdotc_t*>f2py_ptr(wblas.cdotcwrapper._cpointer)
    wcdotu_t *wcdotu = <wcdotu_t*>f2py_ptr(wblas.cdotuwrapper._cpointer)
    wdasum_t *wdasum = <wdasum_t*>f2py_ptr(wblas.dasumwrapper._cpointer)
    wddot_t *wddot = <wddot_t*>f2py_ptr(wblas.ddotwrapper._cpointer)
    wdnrm2_t *wdnrm2 = <wdnrm2_t*>f2py_ptr(wblas.dnrm2wrapper._cpointer)
    wdzasum_t *wdzasum = <wdzasum_t*>f2py_ptr(wblas.dzasumwrapper._cpointer)
    wdznrm2_t *wdznrm2 = <wdznrm2_t*>f2py_ptr(wblas.dznrm2wrapper._cpointer)
    wicamax_t *wicamax = <wicamax_t*>f2py_ptr(wblas.icamaxwrapper._cpointer)
    widamax_t *widamax = <widamax_t*>f2py_ptr(wblas.idamaxwrapper._cpointer)
    wisamax_t *wisamax = <wisamax_t*>f2py_ptr(wblas.isamaxwrapper._cpointer)
    wizamax_t *wizamax = <wizamax_t*>f2py_ptr(wblas.izamaxwrapper._cpointer)
    wsasum_t *wsasum = <wsasum_t*>f2py_ptr(wblas.sasumwrapper._cpointer)
    wscasum_t *wscasum = <wscasum_t*>f2py_ptr(wblas.scasumwrapper._cpointer)
    wscnrm2_t *wscnrm2 = <wscnrm2_t*>f2py_ptr(wblas.scnrm2wrapper._cpointer)
    wsdot_t *wsdot = <wsdot_t*>f2py_ptr(wblas.sdotwrapper._cpointer)
    wsnrm2_t *wsnrm2 = <wsnrm2_t*>f2py_ptr(wblas.snrm2wrapper._cpointer)
    wzdotc_t *wzdotc = <wzdotc_t*>f2py_ptr(wblas.zdotcwrapper._cpointer)
    wzdotu_t *wzdotu = <wzdotu_t*>f2py_ptr(wblas.zdotuwrapper._cpointer)

# Define functions to get the return values from the wrapper subroutines.
cdef c _cdotc(int *n, c *cx, int *incx, c *cy, int *incy) nogil:
    cdef c out
    wcdotc(&out, n, cx, incx, cy, incy)
    return out

cdef c _cdotu(int *n, c *cx, int *incx, c *cy, int *incy) nogil:
    cdef c out
    wcdotu(&out, n, cx, incx, cy, incy)
    return out

cdef d _dasum(int *n, d *dx, int *incx) nogil:
    cdef d out
    wdasum(&out, n, dx, incx)
    return out

cdef d _ddot(int *n, d *dx, int *incx, d *dy, int *incy) nogil:
    cdef d out
    wddot(&out, n, dx, incx, dy, incy)
    return out

cdef d _dnrm2(int *n, d *x, int *incx) nogil:
    cdef d out
    wdnrm2(&out, n, x, incx)
    return out

cdef d _dzasum(int *n, z *zx, int *incx) nogil:
    cdef d out
    wdzasum(&out, n, zx, incx)
    return out

cdef d _dznrm2(int *n, z *x, int *incx) nogil:
    cdef d out
    wdznrm2(&out, n, x, incx)
    return out

# Offset indices by 1 in the i*amax routines to account for the
# difference between Fortran and C indexing.
cdef int _icamax(int *n, c *cx, int *incx) nogil:
    cdef int out
    wicamax(&out, n, cx, incx)
    return out - 1

cdef int _idamax(int *n, d *dx, int *incx) nogil:
    cdef int out
    widamax(&out, n, dx, incx)
    return out - 1

cdef int _isamax(int *n, s *sx, int *incx) nogil:
    cdef int out
    wisamax(&out, n, sx, incx)
    return out - 1

cdef int _izamax(int *n, z *zx, int *incx) nogil:
    cdef int out
    wizamax(&out, n, zx, incx)
    return out - 1

cdef s _sasum(int *n, s *sx, int *incx) nogil:
    cdef s out
    wsasum(&out, n, sx, incx)
    return out

cdef s _scasum(int *n, c *cx, int *incx) nogil:
    cdef s out
    wscasum(&out, n , cx, incx)
    return out

cdef s _scnrm2(int *n, c *x, int *incx) nogil:
    cdef s out
    wscnrm2(&out, n, x, incx)
    return out

cdef s _sdot(int *n, s *sx, int *incx, s *sy, int *incy) nogil:
    cdef s out
    wsdot(&out, n, sx, incx, sy, incy)
    return out

cdef s _snrm2(int *n, s *x, int *incx) nogil:
    cdef s out
    wsnrm2(&out, n, x, incx)
    return out

cdef z _zdotc(int *n, z *zx, int *incx, z *zy, int *incy) nogil:
    cdef z out
    wzdotc(&out, n, zx, incx, zy, incy)
    return out

cdef z _zdotu(int *n, z *zx, int *incx, z *zy, int *incy) nogil:
    cdef z out
    wzdotu(&out, n, zx, incx, zy, incy)
    return out

# When assigning the pointers for routines expected to return a value,
# take the address of the wrapper function instead of getting the
# cpointer attribute of the f2py wrapper.
cdef:
    caxpy_t *caxpy = <caxpy_t*>f2py_ptr(blas.caxpy._cpointer)
    ccopy_t *ccopy = <ccopy_t*>f2py_ptr(blas.ccopy._cpointer)
    cdotc_t *cdotc = &_cdotc
    cdotu_t *cdotu = &_cdotu
    cgemm_t *cgemm = <cgemm_t*>f2py_ptr(blas.cgemm._cpointer)
    cgemv_t *cgemv = <cgemv_t*>f2py_ptr(blas.cgemv._cpointer)
    cgerc_t *cgerc = <cgerc_t*>f2py_ptr(blas.cgerc._cpointer)
    cgeru_t *cgeru = <cgeru_t*>f2py_ptr(blas.cgeru._cpointer)
    chemm_t *chemm = <chemm_t*>f2py_ptr(blas.chemm._cpointer)
    chemv_t *chemv = <chemv_t*>f2py_ptr(blas.chemv._cpointer)
    cher_t *cher = <cher_t*>f2py_ptr(blas.cher._cpointer)
    cher2_t *cher2 = <cher2_t*>f2py_ptr(blas.cher2._cpointer)
    cher2k_t *cher2k = <cher2k_t*>f2py_ptr(blas.cher2k._cpointer)
    cherk_t *cherk = <cherk_t*>f2py_ptr(blas.cherk._cpointer)
    crotg_t *crotg = <crotg_t*>f2py_ptr(blas.crotg._cpointer)
    cscal_t *cscal = <cscal_t*>f2py_ptr(blas.cscal._cpointer)
    csrot_t *csrot = <csrot_t*>f2py_ptr(blas.csrot._cpointer)
    csscal_t *csscal = <csscal_t*>f2py_ptr(blas.csscal._cpointer)
    cswap_t *cswap = <cswap_t*>f2py_ptr(blas.cswap._cpointer)
    csymm_t *csymm = <csymm_t*>f2py_ptr(blas.csymm._cpointer)
    csyr_t *csyr = <csyr_t*>f2py_ptr(blas.csyr._cpointer)
    csyr2k_t *csyr2k = <csyr2k_t*>f2py_ptr(blas.csyr2k._cpointer)
    csyrk_t *csyrk = <csyrk_t*>f2py_ptr(blas.csyrk._cpointer)
    ctrmm_t *ctrmm = <ctrmm_t*>f2py_ptr(blas.ctrmm._cpointer)
    ctrmv_t *ctrmv = <ctrmv_t*>f2py_ptr(blas.ctrmv._cpointer)
    dasum_t *dasum = &_dasum
    daxpy_t *daxpy = <daxpy_t*>f2py_ptr(blas.daxpy._cpointer)
    dcopy_t *dcopy = <dcopy_t*>f2py_ptr(blas.dcopy._cpointer)
    ddot_t *ddot = &_ddot
    dgemm_t *dgemm = <dgemm_t*>f2py_ptr(blas.dgemm._cpointer)
    dgemv_t *dgemv = <dgemv_t*>f2py_ptr(blas.dgemv._cpointer)
    dger_t *dger = <dger_t*>f2py_ptr(blas.dger._cpointer)
    dnrm2_t *dnrm2 = &_dnrm2
    drot_t *drot = <drot_t*>f2py_ptr(blas.drot._cpointer)
    drotg_t *drotg = <drotg_t*>f2py_ptr(blas.drotg._cpointer)
    drotm_t *drotm = <drotm_t*>f2py_ptr(blas.drotm._cpointer)
    drotmg_t *drotmg = <drotmg_t*>f2py_ptr(blas.drotmg._cpointer)
    dscal_t *dscal = <dscal_t*>f2py_ptr(blas.dscal._cpointer)
    dswap_t *dswap = <dswap_t*>f2py_ptr(blas.dswap._cpointer)
    dsymm_t *dsymm = <dsymm_t*>f2py_ptr(blas.dsymm._cpointer)
    dsymv_t *dsymv = <dsymv_t*>f2py_ptr(blas.dsymv._cpointer)
    dsyr_t *dsyr = <dsyr_t*>f2py_ptr(blas.dsyr._cpointer)
    dsyr2_t *dsyr2 = <dsyr2_t*>f2py_ptr(blas.dsyr2._cpointer)
    dsyr2k_t *dsyr2k = <dsyr2k_t*>f2py_ptr(blas.dsyr2k._cpointer)
    dsyrk_t *dsyrk = <dsyrk_t*>f2py_ptr(blas.dsyrk._cpointer)
    dtrmm_t *dtrmm = <dtrmm_t*>f2py_ptr(blas.dtrmm._cpointer)
    dtrmv_t *dtrmv = <dtrmv_t*>f2py_ptr(blas.dtrmv._cpointer)
    dzasum_t *dzasum = &_dzasum
    dznrm2_t *dznrm2 = &_dznrm2
    icamax_t *icamax = &_icamax
    idamax_t *idamax = &_idamax
    isamax_t *isamax = &_isamax
    izamax_t *izamax = &_izamax
    sasum_t *sasum = &_sasum
    saxpy_t *saxpy = <saxpy_t*>f2py_ptr(blas.saxpy._cpointer)
    scasum_t *scasum = &_scasum
    scnrm2_t *scnrm2 = &_scnrm2
    scopy_t *scopy = <scopy_t*>f2py_ptr(blas.scopy._cpointer)
    sdot_t *sdot = &_sdot
    sgemm_t *sgemm = <sgemm_t*>f2py_ptr(blas.sgemm._cpointer)
    sgemv_t *sgemv = <sgemv_t*>f2py_ptr(blas.sgemv._cpointer)
    sger_t *sger = <sger_t*>f2py_ptr(blas.sger._cpointer)
    snrm2_t *snrm2 = &_snrm2
    srot_t *srot = <srot_t*>f2py_ptr(blas.srot._cpointer)
    srotg_t *srotg = <srotg_t*>f2py_ptr(blas.srotg._cpointer)
    srotm_t *srotm = <srotm_t*>f2py_ptr(blas.srotm._cpointer)
    srotmg_t *srotmg = <srotmg_t*>f2py_ptr(blas.srotmg._cpointer)
    sscal_t *sscal = <sscal_t*>f2py_ptr(blas.sscal._cpointer)
    sswap_t *sswap = <sswap_t*>f2py_ptr(blas.sswap._cpointer)
    ssymm_t *ssymm = <ssymm_t*>f2py_ptr(blas.ssymm._cpointer)
    ssymv_t *ssymv = <ssymv_t*>f2py_ptr(blas.ssymv._cpointer)
    ssyr_t *ssyr = <ssyr_t*>f2py_ptr(blas.ssyr._cpointer)
    ssyr2_t *ssyr2 = <ssyr2_t*>f2py_ptr(blas.ssyr2._cpointer)
    ssyr2k_t *ssyr2k = <ssyr2k_t*>f2py_ptr(blas.ssyr2k._cpointer)
    ssyrk_t *ssyrk = <ssyrk_t*>f2py_ptr(blas.ssyrk._cpointer)
    strmm_t *strmm = <strmm_t*>f2py_ptr(blas.strmm._cpointer)
    strmv_t *strmv = <strmv_t*>f2py_ptr(blas.strmv._cpointer)
    zaxpy_t *zaxpy = <zaxpy_t*>f2py_ptr(blas.zaxpy._cpointer)
    zcopy_t *zcopy = <zcopy_t*>f2py_ptr(blas.zcopy._cpointer)
    zdotc_t *zdotc = &_zdotc
    zdotu_t *zdotu = &_zdotu
    zdrot_t *zdrot = <zdrot_t*>f2py_ptr(blas.zdrot._cpointer)
    zdscal_t *zdscal = <zdscal_t*>f2py_ptr(blas.zdscal._cpointer)
    zgemm_t *zgemm = <zgemm_t*>f2py_ptr(blas.zgemm._cpointer)
    zgemv_t *zgemv = <zgemv_t*>f2py_ptr(blas.zgemv._cpointer)
    zgerc_t *zgerc = <zgerc_t*>f2py_ptr(blas.zgerc._cpointer)
    zgeru_t *zgeru = <zgeru_t*>f2py_ptr(blas.zgeru._cpointer)
    zhemm_t *zhemm = <zhemm_t*>f2py_ptr(blas.zhemm._cpointer)
    zhemv_t *zhemv = <zhemv_t*>f2py_ptr(blas.zhemv._cpointer)
    zher_t *zher = <zher_t*>f2py_ptr(blas.zher._cpointer)
    zher2_t *zher2 = <zher2_t*>f2py_ptr(blas.zher2._cpointer)
    zher2k_t *zher2k = <zher2k_t*>f2py_ptr(blas.zher2k._cpointer)
    zherk_t *zherk = <zherk_t*>f2py_ptr(blas.zherk._cpointer)
    zrotg_t *zrotg = <zrotg_t*>f2py_ptr(blas.zrotg._cpointer)
    zscal_t *zscal = <zscal_t*>f2py_ptr(blas.zscal._cpointer)
    zswap_t *zswap = <zswap_t*>f2py_ptr(blas.zswap._cpointer)
    zsymm_t *zsymm = <zsymm_t*>f2py_ptr(blas.zsymm._cpointer)
    zsyr_t *zsyr = <zsyr_t*>f2py_ptr(blas.zsyr._cpointer)
    zsyr2k_t *zsyr2k = <zsyr2k_t*>f2py_ptr(blas.zsyr2k._cpointer)
    zsyrk_t *zsyrk = <zsyrk_t*>f2py_ptr(blas.zsyrk._cpointer)
    ztrmm_t *ztrmm = <ztrmm_t*>f2py_ptr(blas.ztrmm._cpointer)
    ztrmv_t *ztrmv = <ztrmv_t*>f2py_ptr(blas.ztrmv._cpointer)
