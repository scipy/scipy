"""
Function pointers to LAPACK functions.

LAPACK Functions
================

- cgbsv
- cgbtrf
- cgbtrs
- cgebal
- cgees
- cgeev
- cgegv
- cgehrd
- cgelss
- cgeqp3
- cgeqrf
- cgerqf
- cgesdd
- cgesv
- cgetrf
- cgetri
- cgetrs
- cgges
- cggev
- chbevd
- chbevx
- cheev
- cheevd
- cheevr
- chegv
- chegvd
- chegvx
- claswp
- clauum
- cpbsv
- cpbtrf
- cpbtrs
- cposv
- cpotrf
- cpotri
- cpotrs
- ctrsyl
- ctrtri
- ctrtrs
- cungqr
- cungrq
- cunmqr
- dgbsv
- dgbtrf
- dgbtrs
- dgebal
- dgees
- dgeev
- dgegv
- dgehrd
- dgelss
- dgeqp3
- dgeqrf
- dgerqf
- dgesdd
- dgesv
- dgetrf
- dgetri
- dgetrs
- dgges
- dggev
- dlamch
- dlaswp
- dlauum
- dorgqr
- dorgrq
- dormqr
- dpbsv
- dpbtrf
- dpbtrs
- dposv
- dpotrf
- dpotri
- dpotrs
- dsbev
- dsbevd
- dsbevx
- dsyev
- dsyevd
- dsyevr
- dsygv
- dsygvd
- dsygvx
- dtrsyl
- dtrtri
- dtrtrs
- sgbsv
- sgbtrf
- sgbtrs
- sgebal
- sgees
- sgeev
- sgegv
- sgehrd
- sgelss
- sgeqp3
- sgeqrf
- sgerqf
- sgesdd
- sgesv
- sgetrf
- sgetri
- sgetrs
- sgges
- sggev
- slamch
- slaswp
- slauum
- sorgqr
- sorgrq
- sormqr
- spbsv
- spbtrf
- spbtrs
- sposv
- spotrf
- spotri
- spotrs
- ssbev
- ssbevd
- ssbevx
- ssyev
- ssyevd
- ssyevr
- ssygv
- ssygvd
- ssygvx
- strsyl
- strtri
- strtrs
- zgbsv
- zgbtrf
- zgbtrs
- zgebal
- zgees
- zgeev
- zgegv
- zgehrd
- zgelss
- zgeqp3
- zgeqrf
- zgerqf
- zgesdd
- zgesv
- zgetrf
- zgetri
- zgetrs
- zgges
- zggev
- zhbevd
- zhbevx
- zheev
- zheevd
- zheevr
- zhegv
- zhegvd
- zhegvx
- zlaswp
- zlauum
- zpbsv
- zpbtrf
- zpbtrs
- zposv
- zpotrf
- zpotri
- zpotrs
- ztrsyl
- ztrtri
- ztrtrs
- zungqr
- zungrq
- zunmqr

"""

from . import lapack

cdef extern from "f2pyptr.h":
    void *f2py_ptr(object) except NULL

cdef:
    cgbsv_t *cgbsv = <cgbsv_t*>f2py_ptr(lapack.cgbsv._cpointer)
    cgbtrf_t *cgbtrf = <cgbtrf_t*>f2py_ptr(lapack.cgbtrf._cpointer)
    cgbtrs_t *cgbtrs = <cgbtrs_t*>f2py_ptr(lapack.cgbtrs._cpointer)
    cgebal_t *cgebal = <cgebal_t*>f2py_ptr(lapack.cgebal._cpointer)
    cgees_t *cgees = <cgees_t*>f2py_ptr(lapack.cgees._cpointer)
    cgeev_t *cgeev = <cgeev_t*>f2py_ptr(lapack.cgeev._cpointer)
    cgegv_t *cgegv = <cgegv_t*>f2py_ptr(lapack.cgegv._cpointer)
    cgehrd_t *cgehrd = <cgehrd_t*>f2py_ptr(lapack.cgehrd._cpointer)
    cgelss_t *cgelss = <cgelss_t*>f2py_ptr(lapack.cgelss._cpointer)
    cgeqp3_t *cgeqp3 = <cgeqp3_t*>f2py_ptr(lapack.cgeqp3._cpointer)
    cgeqrf_t *cgeqrf = <cgeqrf_t*>f2py_ptr(lapack.cgeqrf._cpointer)
    cgerqf_t *cgerqf = <cgerqf_t*>f2py_ptr(lapack.cgerqf._cpointer)
    cgesdd_t *cgesdd = <cgesdd_t*>f2py_ptr(lapack.cgesdd._cpointer)
    cgesv_t *cgesv = <cgesv_t*>f2py_ptr(lapack.cgesv._cpointer)
    cgetrf_t *cgetrf = <cgetrf_t*>f2py_ptr(lapack.cgetrf._cpointer)
    cgetri_t *cgetri = <cgetri_t*>f2py_ptr(lapack.cgetri._cpointer)
    cgetrs_t *cgetrs = <cgetrs_t*>f2py_ptr(lapack.cgetrs._cpointer)
    cgges_t *cgges = <cgges_t*>f2py_ptr(lapack.cgges._cpointer)
    cggev_t *cggev = <cggev_t*>f2py_ptr(lapack.cggev._cpointer)
    chbevd_t *chbevd = <chbevd_t*>f2py_ptr(lapack.chbevd._cpointer)
    chbevx_t *chbevx = <chbevx_t*>f2py_ptr(lapack.chbevx._cpointer)
    cheev_t *cheev = <cheev_t*>f2py_ptr(lapack.cheev._cpointer)
    cheevd_t *cheevd = <cheevd_t*>f2py_ptr(lapack.cheevd._cpointer)
    cheevr_t *cheevr = <cheevr_t*>f2py_ptr(lapack.cheevr._cpointer)
    chegv_t *chegv = <chegv_t*>f2py_ptr(lapack.chegv._cpointer)
    chegvd_t *chegvd = <chegvd_t*>f2py_ptr(lapack.chegvd._cpointer)
    chegvx_t *chegvx = <chegvx_t*>f2py_ptr(lapack.chegvx._cpointer)
    claswp_t *claswp = <claswp_t*>f2py_ptr(lapack.claswp._cpointer)
    clauum_t *clauum = <clauum_t*>f2py_ptr(lapack.clauum._cpointer)
    cpbsv_t *cpbsv = <cpbsv_t*>f2py_ptr(lapack.cpbsv._cpointer)
    cpbtrf_t *cpbtrf = <cpbtrf_t*>f2py_ptr(lapack.cpbtrf._cpointer)
    cpbtrs_t *cpbtrs = <cpbtrs_t*>f2py_ptr(lapack.cpbtrs._cpointer)
    cposv_t *cposv = <cposv_t*>f2py_ptr(lapack.cposv._cpointer)
    cpotrf_t *cpotrf = <cpotrf_t*>f2py_ptr(lapack.cpotrf._cpointer)
    cpotri_t *cpotri = <cpotri_t*>f2py_ptr(lapack.cpotri._cpointer)
    cpotrs_t *cpotrs = <cpotrs_t*>f2py_ptr(lapack.cpotrs._cpointer)
    ctrsyl_t *ctrsyl = <ctrsyl_t*>f2py_ptr(lapack.ctrsyl._cpointer)
    ctrtri_t *ctrtri = <ctrtri_t*>f2py_ptr(lapack.ctrtri._cpointer)
    ctrtrs_t *ctrtrs = <ctrtrs_t*>f2py_ptr(lapack.ctrtrs._cpointer)
    cungqr_t *cungqr = <cungqr_t*>f2py_ptr(lapack.cungqr._cpointer)
    cungrq_t *cungrq = <cungrq_t*>f2py_ptr(lapack.cungrq._cpointer)
    cunmqr_t *cunmqr = <cunmqr_t*>f2py_ptr(lapack.cunmqr._cpointer)
    dgbsv_t *dgbsv = <dgbsv_t*>f2py_ptr(lapack.dgbsv._cpointer)
    dgbtrf_t *dgbtrf = <dgbtrf_t*>f2py_ptr(lapack.dgbtrf._cpointer)
    dgbtrs_t *dgbtrs = <dgbtrs_t*>f2py_ptr(lapack.dgbtrs._cpointer)
    dgebal_t *dgebal = <dgebal_t*>f2py_ptr(lapack.dgebal._cpointer)
    dgees_t *dgees = <dgees_t*>f2py_ptr(lapack.dgees._cpointer)
    dgeev_t *dgeev = <dgeev_t*>f2py_ptr(lapack.dgeev._cpointer)
    dgegv_t *dgegv = <dgegv_t*>f2py_ptr(lapack.dgegv._cpointer)
    dgehrd_t *dgehrd = <dgehrd_t*>f2py_ptr(lapack.dgehrd._cpointer)
    dgelss_t *dgelss = <dgelss_t*>f2py_ptr(lapack.dgelss._cpointer)
    dgeqp3_t *dgeqp3 = <dgeqp3_t*>f2py_ptr(lapack.dgeqp3._cpointer)
    dgeqrf_t *dgeqrf = <dgeqrf_t*>f2py_ptr(lapack.dgeqrf._cpointer)
    dgerqf_t *dgerqf = <dgerqf_t*>f2py_ptr(lapack.dgerqf._cpointer)
    dgesdd_t *dgesdd = <dgesdd_t*>f2py_ptr(lapack.dgesdd._cpointer)
    dgesv_t *dgesv = <dgesv_t*>f2py_ptr(lapack.dgesv._cpointer)
    dgetrf_t *dgetrf = <dgetrf_t*>f2py_ptr(lapack.dgetrf._cpointer)
    dgetri_t *dgetri = <dgetri_t*>f2py_ptr(lapack.dgetri._cpointer)
    dgetrs_t *dgetrs = <dgetrs_t*>f2py_ptr(lapack.dgetrs._cpointer)
    dgges_t *dgges = <dgges_t*>f2py_ptr(lapack.dgges._cpointer)
    dggev_t *dggev = <dggev_t*>f2py_ptr(lapack.dggev._cpointer)
    dlamch_t *dlamch = <dlamch_t*>f2py_ptr(lapack.dlamch._cpointer)
    dlaswp_t *dlaswp = <dlaswp_t*>f2py_ptr(lapack.dlaswp._cpointer)
    dlauum_t *dlauum = <dlauum_t*>f2py_ptr(lapack.dlauum._cpointer)
    dorgqr_t *dorgqr = <dorgqr_t*>f2py_ptr(lapack.dorgqr._cpointer)
    dorgrq_t *dorgrq = <dorgrq_t*>f2py_ptr(lapack.dorgrq._cpointer)
    dormqr_t *dormqr = <dormqr_t*>f2py_ptr(lapack.dormqr._cpointer)
    dpbsv_t *dpbsv = <dpbsv_t*>f2py_ptr(lapack.dpbsv._cpointer)
    dpbtrf_t *dpbtrf = <dpbtrf_t*>f2py_ptr(lapack.dpbtrf._cpointer)
    dpbtrs_t *dpbtrs = <dpbtrs_t*>f2py_ptr(lapack.dpbtrs._cpointer)
    dposv_t *dposv = <dposv_t*>f2py_ptr(lapack.dposv._cpointer)
    dpotrf_t *dpotrf = <dpotrf_t*>f2py_ptr(lapack.dpotrf._cpointer)
    dpotri_t *dpotri = <dpotri_t*>f2py_ptr(lapack.dpotri._cpointer)
    dpotrs_t *dpotrs = <dpotrs_t*>f2py_ptr(lapack.dpotrs._cpointer)
    dsbev_t *dsbev = <dsbev_t*>f2py_ptr(lapack.dsbev._cpointer)
    dsbevd_t *dsbevd = <dsbevd_t*>f2py_ptr(lapack.dsbevd._cpointer)
    dsbevx_t *dsbevx = <dsbevx_t*>f2py_ptr(lapack.dsbevx._cpointer)
    dsyev_t *dsyev = <dsyev_t*>f2py_ptr(lapack.dsyev._cpointer)
    dsyevd_t *dsyevd = <dsyevd_t*>f2py_ptr(lapack.dsyevd._cpointer)
    dsyevr_t *dsyevr = <dsyevr_t*>f2py_ptr(lapack.dsyevr._cpointer)
    dsygv_t *dsygv = <dsygv_t*>f2py_ptr(lapack.dsygv._cpointer)
    dsygvd_t *dsygvd = <dsygvd_t*>f2py_ptr(lapack.dsygvd._cpointer)
    dsygvx_t *dsygvx = <dsygvx_t*>f2py_ptr(lapack.dsygvx._cpointer)
    dtrsyl_t *dtrsyl = <dtrsyl_t*>f2py_ptr(lapack.dtrsyl._cpointer)
    dtrtri_t *dtrtri = <dtrtri_t*>f2py_ptr(lapack.dtrtri._cpointer)
    dtrtrs_t *dtrtrs = <dtrtrs_t*>f2py_ptr(lapack.dtrtrs._cpointer)
    sgbsv_t *sgbsv = <sgbsv_t*>f2py_ptr(lapack.sgbsv._cpointer)
    sgbtrf_t *sgbtrf = <sgbtrf_t*>f2py_ptr(lapack.sgbtrf._cpointer)
    sgbtrs_t *sgbtrs = <sgbtrs_t*>f2py_ptr(lapack.sgbtrs._cpointer)
    sgebal_t *sgebal = <sgebal_t*>f2py_ptr(lapack.sgebal._cpointer)
    sgees_t *sgees = <sgees_t*>f2py_ptr(lapack.sgees._cpointer)
    sgeev_t *sgeev = <sgeev_t*>f2py_ptr(lapack.sgeev._cpointer)
    sgegv_t *sgegv = <sgegv_t*>f2py_ptr(lapack.sgegv._cpointer)
    sgehrd_t *sgehrd = <sgehrd_t*>f2py_ptr(lapack.sgehrd._cpointer)
    sgelss_t *sgelss = <sgelss_t*>f2py_ptr(lapack.sgelss._cpointer)
    sgeqp3_t *sgeqp3 = <sgeqp3_t*>f2py_ptr(lapack.sgeqp3._cpointer)
    sgeqrf_t *sgeqrf = <sgeqrf_t*>f2py_ptr(lapack.sgeqrf._cpointer)
    sgerqf_t *sgerqf = <sgerqf_t*>f2py_ptr(lapack.sgerqf._cpointer)
    sgesdd_t *sgesdd = <sgesdd_t*>f2py_ptr(lapack.sgesdd._cpointer)
    sgesv_t *sgesv = <sgesv_t*>f2py_ptr(lapack.sgesv._cpointer)
    sgetrf_t *sgetrf = <sgetrf_t*>f2py_ptr(lapack.sgetrf._cpointer)
    sgetri_t *sgetri = <sgetri_t*>f2py_ptr(lapack.sgetri._cpointer)
    sgetrs_t *sgetrs = <sgetrs_t*>f2py_ptr(lapack.sgetrs._cpointer)
    sgges_t *sgges = <sgges_t*>f2py_ptr(lapack.sgges._cpointer)
    sggev_t *sggev = <sggev_t*>f2py_ptr(lapack.sggev._cpointer)
    slamch_t *slamch = <slamch_t*>f2py_ptr(lapack.wslamch._cpointer)
    slaswp_t *slaswp = <slaswp_t*>f2py_ptr(lapack.slaswp._cpointer)
    slauum_t *slauum = <slauum_t*>f2py_ptr(lapack.slauum._cpointer)
    sorgqr_t *sorgqr = <sorgqr_t*>f2py_ptr(lapack.sorgqr._cpointer)
    sorgrq_t *sorgrq = <sorgrq_t*>f2py_ptr(lapack.sorgrq._cpointer)
    sormqr_t *sormqr = <sormqr_t*>f2py_ptr(lapack.sormqr._cpointer)
    spbsv_t *spbsv = <spbsv_t*>f2py_ptr(lapack.spbsv._cpointer)
    spbtrf_t *spbtrf = <spbtrf_t*>f2py_ptr(lapack.spbtrf._cpointer)
    spbtrs_t *spbtrs = <spbtrs_t*>f2py_ptr(lapack.spbtrs._cpointer)
    sposv_t *sposv = <sposv_t*>f2py_ptr(lapack.sposv._cpointer)
    spotrf_t *spotrf = <spotrf_t*>f2py_ptr(lapack.spotrf._cpointer)
    spotri_t *spotri = <spotri_t*>f2py_ptr(lapack.spotri._cpointer)
    spotrs_t *spotrs = <spotrs_t*>f2py_ptr(lapack.spotrs._cpointer)
    ssbev_t *ssbev = <ssbev_t*>f2py_ptr(lapack.ssbev._cpointer)
    ssbevd_t *ssbevd = <ssbevd_t*>f2py_ptr(lapack.ssbevd._cpointer)
    ssbevx_t *ssbevx = <ssbevx_t*>f2py_ptr(lapack.ssbevx._cpointer)
    ssyev_t *ssyev = <ssyev_t*> f2py_ptr(lapack.ssyev._cpointer)
    ssyevd_t *ssyevd = <ssyevd_t*>f2py_ptr(lapack.ssyevd._cpointer)
    ssyevr_t *ssyevr = <ssyevr_t*>f2py_ptr(lapack.ssyevr._cpointer)
    ssygv_t *ssygv = <ssygv_t*>f2py_ptr(lapack.ssygv._cpointer)
    ssygvd_t *ssygvd = <ssygvd_t*>f2py_ptr(lapack.ssygvd._cpointer)
    ssygvx_t *ssygvx = <ssygvx_t*>f2py_ptr(lapack.ssygvx._cpointer)
    strsyl_t *strsyl = <strsyl_t*>f2py_ptr(lapack.strsyl._cpointer)
    strtri_t *strtri = <strtri_t*>f2py_ptr(lapack.strtri._cpointer)
    strtrs_t *strtrs = <strtrs_t*>f2py_ptr(lapack.strtrs._cpointer)
    zgbsv_t *zgbsv = <zgbsv_t*>f2py_ptr(lapack.zgbsv._cpointer)
    zgbtrf_t *zgbtrf = <zgbtrf_t*>f2py_ptr(lapack.zgbtrf._cpointer)
    zgbtrs_t *zgbtrs = <zgbtrs_t*>f2py_ptr(lapack.zgbtrs._cpointer)
    zgebal_t *zgebal = <zgebal_t*>f2py_ptr(lapack.zgebal._cpointer)
    zgees_t *zgees = <zgees_t*>f2py_ptr(lapack.zgees._cpointer)
    zgeev_t *zgeev = <zgeev_t*>f2py_ptr(lapack.zgeev._cpointer)
    zgegv_t *zgegv = <zgegv_t*>f2py_ptr(lapack.zgegv._cpointer)
    zgehrd_t *zgehrd = <zgehrd_t*>f2py_ptr(lapack.zgehrd._cpointer)
    zgelss_t *zgelss = <zgelss_t*>f2py_ptr(lapack.zgelss._cpointer)
    zgeqp3_t *zgeqp3 = <zgeqp3_t*>f2py_ptr(lapack.zgeqp3._cpointer)
    zgeqrf_t *zgeqrf = <zgeqrf_t*>f2py_ptr(lapack.zgeqrf._cpointer)
    zgerqf_t *zgerqf = <zgerqf_t*>f2py_ptr(lapack.zgerqf._cpointer)
    zgesdd_t *zgesdd = <zgesdd_t*>f2py_ptr(lapack.zgesdd._cpointer)
    zgesv_t *zgesv = <zgesv_t*>f2py_ptr(lapack.zgesv._cpointer)
    zgetrf_t *zgetrf = <zgetrf_t*>f2py_ptr(lapack.zgetrf._cpointer)
    zgetri_t *zgetri = <zgetri_t*>f2py_ptr(lapack.zgetri._cpointer)
    zgetrs_t *zgetrs = <zgetrs_t*>f2py_ptr(lapack.zgetrs._cpointer)
    zgges_t *zgges = <zgges_t*>f2py_ptr(lapack.zgges._cpointer)
    zggev_t *zggev = <zggev_t*>f2py_ptr(lapack.zggev._cpointer)
    zhbevd_t *zhbevd = <zhbevd_t*>f2py_ptr(lapack.zhbevd._cpointer)
    zhbevx_t *zhbevx = <zhbevx_t*>f2py_ptr(lapack.zhbevx._cpointer)
    zheev_t *zheev = <zheev_t*>f2py_ptr(lapack.zheev._cpointer)
    zheevd_t *zheevd = <zheevd_t*>f2py_ptr(lapack.zheevd._cpointer)
    zheevr_t *zheevr = <zheevr_t*>f2py_ptr(lapack.zheevr._cpointer)
    zhegv_t *zhegv = <zhegv_t*>f2py_ptr(lapack.zhegv._cpointer)
    zhegvd_t *zhegvd = <zhegvd_t*>f2py_ptr(lapack.zhegvd._cpointer)
    zhegvx_t *zhegvx = <zhegvx_t*>f2py_ptr(lapack.zhegvx._cpointer)
    zlaswp_t *zlaswp = <zlaswp_t*>f2py_ptr(lapack.zlaswp._cpointer)
    zlauum_t *zlauum = <zlauum_t*>f2py_ptr(lapack.zlauum._cpointer)
    zpbsv_t *zpbsv = <zpbsv_t*>f2py_ptr(lapack.zpbsv._cpointer)
    zpbtrf_t *zpbtrf = <zpbtrf_t*>f2py_ptr(lapack.zpbtrf._cpointer)
    zpbtrs_t *zpbtrs = <zpbtrs_t*>f2py_ptr(lapack.zpbtrs._cpointer)
    zposv_t *zposv = <zposv_t*>f2py_ptr(lapack.zposv._cpointer)
    zpotrf_t *zpotrf = <zpotrf_t*>f2py_ptr(lapack.zpotrf._cpointer)
    zpotri_t *zpotri = <zpotri_t*>f2py_ptr(lapack.zpotri._cpointer)
    zpotrs_t *zpotrs = <zpotrs_t*>f2py_ptr(lapack.zpotrs._cpointer)
    ztrsyl_t *ztrsyl = <ztrsyl_t*>f2py_ptr(lapack.ztrsyl._cpointer)
    ztrtri_t *ztrtri = <ztrtri_t*>f2py_ptr(lapack.ztrtri._cpointer)
    ztrtrs_t *ztrtrs = <ztrtrs_t*>f2py_ptr(lapack.ztrtrs._cpointer)
    zungqr_t *zungqr = <zungqr_t*>f2py_ptr(lapack.zungqr._cpointer)
    zungrq_t *zungrq = <zungrq_t*>f2py_ptr(lapack.zungrq._cpointer)
    zunmqr_t *zunmqr = <zunmqr_t*>f2py_ptr(lapack.zunmqr._cpointer)
