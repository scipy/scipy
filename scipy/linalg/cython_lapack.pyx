"""
LAPACK functions for Cython
===========================

Usable from Cython via::

    cimport scipy.linalg.cython_lapack

This module provides Cython-level wrappers for all primary routines included
in LAPACK 3.1.0 and some of the fixed-api auxiliary routines.

The signature for dcgesv changed from LAPACK 3.1.1 to LAPACK 3.2.0.
The version here is the newer of the two since it matches the signature
from later versions of LAPACK and the version in the CLAPACK included in OSX.

Raw function pointers (Fortran-style pointer arguments):

- cbdsqr
- cgbbrd
- cgbcon
- cgbequ
- cgbrfs
- cgbsv
- cgbsvx
- cgbtf2
- cgbtrf
- cgbtrs
- cgebak
- cgebal
- cgebd2
- cgebrd
- cgecon
- cgeequ
- cgees
- cgeesx
- cgeev
- cgeevx
- cgegs
- cgegv
- cgehd2
- cgehrd
- cgelq2
- cgelqf
- cgels
- cgelsd
- cgelss
- cgelsx
- cgelsy
- cgeql2
- cgeqlf
- cgeqp3
- cgeqpf
- cgeqr2
- cgeqrf
- cgerfs
- cgerq2
- cgerqf
- cgesc2
- cgesdd
- cgesv
- cgesvd
- cgesvx
- cgetc2
- cgetf2
- cgetrf
- cgetri
- cgetrs
- cggbak
- cggbal
- cgges
- cggesx
- cggev
- cggevx
- cggglm
- cgghrd
- cgglse
- cggqrf
- cggrqf
- cggsvd
- cggsvp
- cgtcon
- cgtrfs
- cgtsv
- cgtsvx
- cgttrf
- cgttrs
- cgtts2
- chbev
- chbevd
- chbevx
- chbgst
- chbgv
- chbgvd
- chbgvx
- chbtrd
- checon
- cheev
- cheevd
- cheevr
- cheevx
- chegs2
- chegst
- chegv
- chegvd
- chegvx
- cherfs
- chesv
- chesvx
- chetd2
- chetf2
- chetrd
- chetrf
- chetri
- chetrs
- chgeqz
- chpcon
- chpev
- chpevd
- chpevx
- chpgst
- chpgv
- chpgvd
- chpgvx
- chprfs
- chpsv
- chpsvx
- chptrd
- chptrf
- chptri
- chptrs
- chsein
- chseqr
- clacn2
- clacon
- clangb
- clange
- clangt
- clanhb
- clanhe
- clanhp
- clanhs
- clanht
- clansb
- clansp
- clansy
- clantb
- clantp
- clantr
- clarf
- clarz
- claswp
- clauum
- cpbcon
- cpbequ
- cpbrfs
- cpbstf
- cpbsv
- cpbsvx
- cpbtf2
- cpbtrf
- cpbtrs
- cpocon
- cpoequ
- cporfs
- cposv
- cposvx
- cpotf2
- cpotrf
- cpotri
- cpotrs
- cppcon
- cppequ
- cpprfs
- cppsv
- cppsvx
- cpptrf
- cpptri
- cpptrs
- cptcon
- cpteqr
- cptrfs
- cptsv
- cptsvx
- cpttrf
- cpttrs
- cptts2
- crot
- cspcon
- cspmv
- cspr
- csprfs
- cspsv
- cspsvx
- csptrf
- csptri
- csptrs
- csrscl
- cstedc
- cstegr
- cstein
- cstemr
- csteqr
- csycon
- csymv
- csyr
- csyrfs
- csysv
- csysvx
- csytf2
- csytrf
- csytri
- csytrs
- ctbcon
- ctbrfs
- ctbtrs
- ctgevc
- ctgex2
- ctgexc
- ctgsen
- ctgsja
- ctgsna
- ctgsy2
- ctgsyl
- ctpcon
- ctprfs
- ctptri
- ctptrs
- ctrcon
- ctrevc
- ctrexc
- ctrrfs
- ctrsen
- ctrsna
- ctrsyl
- ctrti2
- ctrtri
- ctrtrs
- ctzrqf
- ctzrzf
- cung2l
- cung2r
- cungbr
- cunghr
- cungl2
- cunglq
- cungql
- cungqr
- cungr2
- cungrq
- cungtr
- cunm2l
- cunm2r
- cunmbr
- cunmhr
- cunml2
- cunmlq
- cunmql
- cunmqr
- cunmr2
- cunmr3
- cunmrq
- cunmrz
- cunmtr
- cupgtr
- cupmtr
- dbdsdc
- dbdsqr
- ddisna
- dgbbrd
- dgbcon
- dgbequ
- dgbrfs
- dgbsv
- dgbsvx
- dgbtf2
- dgbtrf
- dgbtrs
- dgebak
- dgebal
- dgebd2
- dgebrd
- dgecon
- dgeequ
- dgees
- dgeesx
- dgeev
- dgeevx
- dgegs
- dgegv
- dgehd2
- dgehrd
- dgelq2
- dgelqf
- dgels
- dgelsd
- dgelss
- dgelsx
- dgelsy
- dgeql2
- dgeqlf
- dgeqp3
- dgeqpf
- dgeqr2
- dgeqrf
- dgerfs
- dgerq2
- dgerqf
- dgesc2
- dgesdd
- dgesv
- dgesvd
- dgesvx
- dgetc2
- dgetf2
- dgetrf
- dgetri
- dgetrs
- dggbak
- dggbal
- dgges
- dggesx
- dggev
- dggevx
- dggglm
- dgghrd
- dgglse
- dggqrf
- dggrqf
- dggsvd
- dggsvp
- dgtcon
- dgtrfs
- dgtsv
- dgtsvx
- dgttrf
- dgttrs
- dgtts2
- dhgeqz
- dhsein
- dhseqr
- disnan
- dlacn2
- dlacon
- dlamch
- dlangb
- dlange
- dlangt
- dlanhs
- dlansb
- dlansp
- dlanst
- dlansy
- dlantb
- dlantp
- dlantr
- dlanv2
- dlarf
- dlarz
- dlaswp
- dlauum
- dopgtr
- dopmtr
- dorg2l
- dorg2r
- dorgbr
- dorghr
- dorgl2
- dorglq
- dorgql
- dorgqr
- dorgr2
- dorgrq
- dorgtr
- dorm2l
- dorm2r
- dormbr
- dormhr
- dorml2
- dormlq
- dormql
- dormqr
- dormr2
- dormr3
- dormrq
- dormrz
- dormtr
- dpbcon
- dpbequ
- dpbrfs
- dpbstf
- dpbsv
- dpbsvx
- dpbtf2
- dpbtrf
- dpbtrs
- dpocon
- dpoequ
- dporfs
- dposv
- dposvx
- dpotf2
- dpotrf
- dpotri
- dpotrs
- dppcon
- dppequ
- dpprfs
- dppsv
- dppsvx
- dpptrf
- dpptri
- dpptrs
- dptcon
- dpteqr
- dptrfs
- dptsv
- dptsvx
- dpttrf
- dpttrs
- dptts2
- drscl
- dsbev
- dsbevd
- dsbevx
- dsbgst
- dsbgv
- dsbgvd
- dsbgvx
- dsbtrd
- dsgesv
- dspcon
- dspev
- dspevd
- dspevx
- dspgst
- dspgv
- dspgvd
- dspgvx
- dsprfs
- dspsv
- dspsvx
- dsptrd
- dsptrf
- dsptri
- dsptrs
- dstebz
- dstedc
- dstegr
- dstein
- dstemr
- dsteqr
- dsterf
- dstev
- dstevd
- dstevr
- dstevx
- dsycon
- dsyev
- dsyevd
- dsyevr
- dsyevx
- dsygs2
- dsygst
- dsygv
- dsygvd
- dsygvx
- dsyrfs
- dsysv
- dsysvx
- dsytd2
- dsytf2
- dsytrd
- dsytrf
- dsytri
- dsytrs
- dtbcon
- dtbrfs
- dtbtrs
- dtgevc
- dtgex2
- dtgexc
- dtgsen
- dtgsja
- dtgsna
- dtgsy2
- dtgsyl
- dtpcon
- dtprfs
- dtptri
- dtptrs
- dtrcon
- dtrevc
- dtrexc
- dtrrfs
- dtrsen
- dtrsna
- dtrsyl
- dtrti2
- dtrtri
- dtrtrs
- dtzrqf
- dtzrzf
- dzsum1
- icmax1
- ieeeck
- ilaenv
- iparmq
- izmax1
- lsamen
- sbdsdc
- sbdsqr
- scsum1
- sdisna
- sgbbrd
- sgbcon
- sgbequ
- sgbrfs
- sgbsv
- sgbsvx
- sgbtf2
- sgbtrf
- sgbtrs
- sgebak
- sgebal
- sgebd2
- sgebrd
- sgecon
- sgeequ
- sgees
- sgeesx
- sgeev
- sgeevx
- sgegs
- sgegv
- sgehd2
- sgehrd
- sgelq2
- sgelqf
- sgels
- sgelsd
- sgelss
- sgelsx
- sgelsy
- sgeql2
- sgeqlf
- sgeqp3
- sgeqpf
- sgeqr2
- sgeqrf
- sgerfs
- sgerq2
- sgerqf
- sgesc2
- sgesdd
- sgesv
- sgesvd
- sgesvx
- sgetc2
- sgetf2
- sgetrf
- sgetri
- sgetrs
- sggbak
- sggbal
- sgges
- sggesx
- sggev
- sggevx
- sggglm
- sgghrd
- sgglse
- sggqrf
- sggrqf
- sggsvd
- sggsvp
- sgtcon
- sgtrfs
- sgtsv
- sgtsvx
- sgttrf
- sgttrs
- sgtts2
- shgeqz
- shsein
- shseqr
- slacn2
- slacon
- slamch
- slangb
- slange
- slangt
- slanhs
- slansb
- slansp
- slanst
- slansy
- slantb
- slantp
- slantr
- slanv2
- slarf
- slarz
- slaswp
- slauum
- sopgtr
- sopmtr
- sorg2l
- sorg2r
- sorgbr
- sorghr
- sorgl2
- sorglq
- sorgql
- sorgqr
- sorgr2
- sorgrq
- sorgtr
- sorm2l
- sorm2r
- sormbr
- sormhr
- sorml2
- sormlq
- sormql
- sormqr
- sormr2
- sormr3
- sormrq
- sormrz
- sormtr
- spbcon
- spbequ
- spbrfs
- spbstf
- spbsv
- spbsvx
- spbtf2
- spbtrf
- spbtrs
- spocon
- spoequ
- sporfs
- sposv
- sposvx
- spotf2
- spotrf
- spotri
- spotrs
- sppcon
- sppequ
- spprfs
- sppsv
- sppsvx
- spptrf
- spptri
- spptrs
- sptcon
- spteqr
- sptrfs
- sptsv
- sptsvx
- spttrf
- spttrs
- sptts2
- srscl
- ssbev
- ssbevd
- ssbevx
- ssbgst
- ssbgv
- ssbgvd
- ssbgvx
- ssbtrd
- sspcon
- sspev
- sspevd
- sspevx
- sspgst
- sspgv
- sspgvd
- sspgvx
- ssprfs
- sspsv
- sspsvx
- ssptrd
- ssptrf
- ssptri
- ssptrs
- sstebz
- sstedc
- sstegr
- sstein
- sstemr
- ssteqr
- ssterf
- sstev
- sstevd
- sstevr
- sstevx
- ssycon
- ssyev
- ssyevd
- ssyevr
- ssyevx
- ssygs2
- ssygst
- ssygv
- ssygvd
- ssygvx
- ssyrfs
- ssysv
- ssysvx
- ssytd2
- ssytf2
- ssytrd
- ssytrf
- ssytri
- ssytrs
- stbcon
- stbrfs
- stbtrs
- stgevc
- stgex2
- stgexc
- stgsen
- stgsja
- stgsna
- stgsy2
- stgsyl
- stpcon
- stprfs
- stptri
- stptrs
- strcon
- strevc
- strexc
- strrfs
- strsen
- strsna
- strsyl
- strti2
- strtri
- strtrs
- stzrqf
- stzrzf
- xerbla
- zbdsqr
- zcgesv
- zdrscl
- zgbbrd
- zgbcon
- zgbequ
- zgbrfs
- zgbsv
- zgbsvx
- zgbtf2
- zgbtrf
- zgbtrs
- zgebak
- zgebal
- zgebd2
- zgebrd
- zgecon
- zgeequ
- zgees
- zgeesx
- zgeev
- zgeevx
- zgegs
- zgegv
- zgehd2
- zgehrd
- zgelq2
- zgelqf
- zgels
- zgelsd
- zgelss
- zgelsx
- zgelsy
- zgeql2
- zgeqlf
- zgeqp3
- zgeqpf
- zgeqr2
- zgeqrf
- zgerfs
- zgerq2
- zgerqf
- zgesc2
- zgesdd
- zgesv
- zgesvd
- zgesvx
- zgetc2
- zgetf2
- zgetrf
- zgetri
- zgetrs
- zggbak
- zggbal
- zgges
- zggesx
- zggev
- zggevx
- zggglm
- zgghrd
- zgglse
- zggqrf
- zggrqf
- zggsvd
- zggsvp
- zgtcon
- zgtrfs
- zgtsv
- zgtsvx
- zgttrf
- zgttrs
- zgtts2
- zhbev
- zhbevd
- zhbevx
- zhbgst
- zhbgv
- zhbgvd
- zhbgvx
- zhbtrd
- zhecon
- zheev
- zheevd
- zheevr
- zheevx
- zhegs2
- zhegst
- zhegv
- zhegvd
- zhegvx
- zherfs
- zhesv
- zhesvx
- zhetd2
- zhetf2
- zhetrd
- zhetrf
- zhetri
- zhetrs
- zhgeqz
- zhpcon
- zhpev
- zhpevd
- zhpevx
- zhpgst
- zhpgv
- zhpgvd
- zhpgvx
- zhprfs
- zhpsv
- zhpsvx
- zhptrd
- zhptrf
- zhptri
- zhptrs
- zhsein
- zhseqr
- zlacn2
- zlacon
- zlangb
- zlange
- zlangt
- zlanhb
- zlanhe
- zlanhp
- zlanhs
- zlanht
- zlansb
- zlansp
- zlansy
- zlantb
- zlantp
- zlantr
- zlarf
- zlarz
- zlaswp
- zlauum
- zpbcon
- zpbequ
- zpbrfs
- zpbstf
- zpbsv
- zpbsvx
- zpbtf2
- zpbtrf
- zpbtrs
- zpocon
- zpoequ
- zporfs
- zposv
- zposvx
- zpotf2
- zpotrf
- zpotri
- zpotrs
- zppcon
- zppequ
- zpprfs
- zppsv
- zppsvx
- zpptrf
- zpptri
- zpptrs
- zptcon
- zpteqr
- zptrfs
- zptsv
- zptsvx
- zpttrf
- zpttrs
- zptts2
- zrot
- zspcon
- zspmv
- zspr
- zsprfs
- zspsv
- zspsvx
- zsptrf
- zsptri
- zsptrs
- zstedc
- zstegr
- zstein
- zstemr
- zsteqr
- zsycon
- zsymv
- zsyr
- zsyrfs
- zsysv
- zsysvx
- zsytf2
- zsytrf
- zsytri
- zsytrs
- ztbcon
- ztbrfs
- ztbtrs
- ztgevc
- ztgex2
- ztgexc
- ztgsen
- ztgsja
- ztgsna
- ztgsy2
- ztgsyl
- ztpcon
- ztprfs
- ztptri
- ztptrs
- ztrcon
- ztrevc
- ztrexc
- ztrrfs
- ztrsen
- ztrsna
- ztrsyl
- ztrti2
- ztrtri
- ztrtrs
- ztzrqf
- ztzrzf
- zung2l
- zung2r
- zungbr
- zunghr
- zungl2
- zunglq
- zungql
- zungqr
- zungr2
- zungrq
- zungtr
- zunm2l
- zunm2r
- zunmbr
- zunmhr
- zunml2
- zunmlq
- zunmql
- zunmqr
- zunmr2
- zunmr3
- zunmrq
- zunmrz
- zunmtr
- zupgtr
- zupmtr


"""

cdef extern from "fortran_defs.h":
    pass

ctypedef void _wclangb_t(s *out, char *norm, int *n, int *kl, int *ku, c *ab, int *ldab, s *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_clangb "F_FUNC(clangbwrp, CLANGBWRP)"(s *out, char *norm, int *n, int *kl, int *ku, c *ab, int *ldab, s *work) nogil
cdef s _wrap_clangb(char *norm, int *n, int *kl, int *ku, c *ab, int *ldab, s *work) nogil:
    cdef s out
    _fortran_clangb(&out, norm, n, kl, ku, ab, ldab, work)
    return out
cdef clangb_t *clangb_f = &_wrap_clangb

ctypedef void _wclange_t(s *out, char *norm, int *m, int *n, c *a, int *lda, s *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_clange "F_FUNC(clangewrp, CLANGEWRP)"(s *out, char *norm, int *m, int *n, c *a, int *lda, s *work) nogil
cdef s _wrap_clange(char *norm, int *m, int *n, c *a, int *lda, s *work) nogil:
    cdef s out
    _fortran_clange(&out, norm, m, n, a, lda, work)
    return out
cdef clange_t *clange_f = &_wrap_clange

ctypedef void _wclangt_t(s *out, char *norm, int *n, c *dl, c *d, c *du) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_clangt "F_FUNC(clangtwrp, CLANGTWRP)"(s *out, char *norm, int *n, c *dl, c *d, c *du) nogil
cdef s _wrap_clangt(char *norm, int *n, c *dl, c *d, c *du) nogil:
    cdef s out
    _fortran_clangt(&out, norm, n, dl, d, du)
    return out
cdef clangt_t *clangt_f = &_wrap_clangt

ctypedef void _wclanhb_t(s *out, char *norm, char *uplo, int *n, int *k, c *ab, int *ldab, s *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_clanhb "F_FUNC(clanhbwrp, CLANHBWRP)"(s *out, char *norm, char *uplo, int *n, int *k, c *ab, int *ldab, s *work) nogil
cdef s _wrap_clanhb(char *norm, char *uplo, int *n, int *k, c *ab, int *ldab, s *work) nogil:
    cdef s out
    _fortran_clanhb(&out, norm, uplo, n, k, ab, ldab, work)
    return out
cdef clanhb_t *clanhb_f = &_wrap_clanhb

ctypedef void _wclanhe_t(s *out, char *norm, char *uplo, int *n, c *a, int *lda, s *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_clanhe "F_FUNC(clanhewrp, CLANHEWRP)"(s *out, char *norm, char *uplo, int *n, c *a, int *lda, s *work) nogil
cdef s _wrap_clanhe(char *norm, char *uplo, int *n, c *a, int *lda, s *work) nogil:
    cdef s out
    _fortran_clanhe(&out, norm, uplo, n, a, lda, work)
    return out
cdef clanhe_t *clanhe_f = &_wrap_clanhe

ctypedef void _wclanhp_t(s *out, char *norm, char *uplo, int *n, c *ap, s *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_clanhp "F_FUNC(clanhpwrp, CLANHPWRP)"(s *out, char *norm, char *uplo, int *n, c *ap, s *work) nogil
cdef s _wrap_clanhp(char *norm, char *uplo, int *n, c *ap, s *work) nogil:
    cdef s out
    _fortran_clanhp(&out, norm, uplo, n, ap, work)
    return out
cdef clanhp_t *clanhp_f = &_wrap_clanhp

ctypedef void _wclanhs_t(s *out, char *norm, int *n, c *a, int *lda, s *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_clanhs "F_FUNC(clanhswrp, CLANHSWRP)"(s *out, char *norm, int *n, c *a, int *lda, s *work) nogil
cdef s _wrap_clanhs(char *norm, int *n, c *a, int *lda, s *work) nogil:
    cdef s out
    _fortran_clanhs(&out, norm, n, a, lda, work)
    return out
cdef clanhs_t *clanhs_f = &_wrap_clanhs

ctypedef void _wclanht_t(s *out, char *norm, int *n, s *d, c *e) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_clanht "F_FUNC(clanhtwrp, CLANHTWRP)"(s *out, char *norm, int *n, s *d, c *e) nogil
cdef s _wrap_clanht(char *norm, int *n, s *d, c *e) nogil:
    cdef s out
    _fortran_clanht(&out, norm, n, d, e)
    return out
cdef clanht_t *clanht_f = &_wrap_clanht

ctypedef void _wclansb_t(s *out, char *norm, char *uplo, int *n, int *k, c *ab, int *ldab, s *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_clansb "F_FUNC(clansbwrp, CLANSBWRP)"(s *out, char *norm, char *uplo, int *n, int *k, c *ab, int *ldab, s *work) nogil
cdef s _wrap_clansb(char *norm, char *uplo, int *n, int *k, c *ab, int *ldab, s *work) nogil:
    cdef s out
    _fortran_clansb(&out, norm, uplo, n, k, ab, ldab, work)
    return out
cdef clansb_t *clansb_f = &_wrap_clansb

ctypedef void _wclansp_t(s *out, char *norm, char *uplo, int *n, c *ap, s *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_clansp "F_FUNC(clanspwrp, CLANSPWRP)"(s *out, char *norm, char *uplo, int *n, c *ap, s *work) nogil
cdef s _wrap_clansp(char *norm, char *uplo, int *n, c *ap, s *work) nogil:
    cdef s out
    _fortran_clansp(&out, norm, uplo, n, ap, work)
    return out
cdef clansp_t *clansp_f = &_wrap_clansp

ctypedef void _wclansy_t(s *out, char *norm, char *uplo, int *n, c *a, int *lda, s *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_clansy "F_FUNC(clansywrp, CLANSYWRP)"(s *out, char *norm, char *uplo, int *n, c *a, int *lda, s *work) nogil
cdef s _wrap_clansy(char *norm, char *uplo, int *n, c *a, int *lda, s *work) nogil:
    cdef s out
    _fortran_clansy(&out, norm, uplo, n, a, lda, work)
    return out
cdef clansy_t *clansy_f = &_wrap_clansy

ctypedef void _wclantb_t(s *out, char *norm, char *uplo, char *diag, int *n, int *k, c *ab, int *ldab, s *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_clantb "F_FUNC(clantbwrp, CLANTBWRP)"(s *out, char *norm, char *uplo, char *diag, int *n, int *k, c *ab, int *ldab, s *work) nogil
cdef s _wrap_clantb(char *norm, char *uplo, char *diag, int *n, int *k, c *ab, int *ldab, s *work) nogil:
    cdef s out
    _fortran_clantb(&out, norm, uplo, diag, n, k, ab, ldab, work)
    return out
cdef clantb_t *clantb_f = &_wrap_clantb

ctypedef void _wclantp_t(s *out, char *norm, char *uplo, char *diag, int *n, c *ap, s *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_clantp "F_FUNC(clantpwrp, CLANTPWRP)"(s *out, char *norm, char *uplo, char *diag, int *n, c *ap, s *work) nogil
cdef s _wrap_clantp(char *norm, char *uplo, char *diag, int *n, c *ap, s *work) nogil:
    cdef s out
    _fortran_clantp(&out, norm, uplo, diag, n, ap, work)
    return out
cdef clantp_t *clantp_f = &_wrap_clantp

ctypedef void _wclantr_t(s *out, char *norm, char *uplo, char *diag, int *m, int *n, c *a, int *lda, s *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_clantr "F_FUNC(clantrwrp, CLANTRWRP)"(s *out, char *norm, char *uplo, char *diag, int *m, int *n, c *a, int *lda, s *work) nogil
cdef s _wrap_clantr(char *norm, char *uplo, char *diag, int *m, int *n, c *a, int *lda, s *work) nogil:
    cdef s out
    _fortran_clantr(&out, norm, uplo, diag, m, n, a, lda, work)
    return out
cdef clantr_t *clantr_f = &_wrap_clantr

ctypedef void _wdisnan_t(bint *out, d *din) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_disnan "F_FUNC(disnanwrp, DISNANWRP)"(bint *out, d *din) nogil
cdef bint _wrap_disnan(d *din) nogil:
    cdef bint out
    _fortran_disnan(&out, din)
    return out
cdef disnan_t *disnan_f = &_wrap_disnan

ctypedef void _wdlamch_t(d *out, char *cmach) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dlamch "F_FUNC(dlamchwrp, DLAMCHWRP)"(d *out, char *cmach) nogil
cdef d _wrap_dlamch(char *cmach) nogil:
    cdef d out
    _fortran_dlamch(&out, cmach)
    return out
cdef dlamch_t *dlamch_f = &_wrap_dlamch

ctypedef void _wdlangb_t(d *out, char *norm, int *n, int *kl, int *ku, d *ab, int *ldab, d *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dlangb "F_FUNC(dlangbwrp, DLANGBWRP)"(d *out, char *norm, int *n, int *kl, int *ku, d *ab, int *ldab, d *work) nogil
cdef d _wrap_dlangb(char *norm, int *n, int *kl, int *ku, d *ab, int *ldab, d *work) nogil:
    cdef d out
    _fortran_dlangb(&out, norm, n, kl, ku, ab, ldab, work)
    return out
cdef dlangb_t *dlangb_f = &_wrap_dlangb

ctypedef void _wdlange_t(d *out, char *norm, int *m, int *n, d *a, int *lda, d *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dlange "F_FUNC(dlangewrp, DLANGEWRP)"(d *out, char *norm, int *m, int *n, d *a, int *lda, d *work) nogil
cdef d _wrap_dlange(char *norm, int *m, int *n, d *a, int *lda, d *work) nogil:
    cdef d out
    _fortran_dlange(&out, norm, m, n, a, lda, work)
    return out
cdef dlange_t *dlange_f = &_wrap_dlange

ctypedef void _wdlangt_t(d *out, char *norm, int *n, d *dl, d *d_, d *du) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dlangt "F_FUNC(dlangtwrp, DLANGTWRP)"(d *out, char *norm, int *n, d *dl, d *d_, d *du) nogil
cdef d _wrap_dlangt(char *norm, int *n, d *dl, d *d_, d *du) nogil:
    cdef d out
    _fortran_dlangt(&out, norm, n, dl, d_, du)
    return out
cdef dlangt_t *dlangt_f = &_wrap_dlangt

ctypedef void _wdlanhs_t(d *out, char *norm, int *n, d *a, int *lda, d *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dlanhs "F_FUNC(dlanhswrp, DLANHSWRP)"(d *out, char *norm, int *n, d *a, int *lda, d *work) nogil
cdef d _wrap_dlanhs(char *norm, int *n, d *a, int *lda, d *work) nogil:
    cdef d out
    _fortran_dlanhs(&out, norm, n, a, lda, work)
    return out
cdef dlanhs_t *dlanhs_f = &_wrap_dlanhs

ctypedef void _wdlansb_t(d *out, char *norm, char *uplo, int *n, int *k, d *ab, int *ldab, d *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dlansb "F_FUNC(dlansbwrp, DLANSBWRP)"(d *out, char *norm, char *uplo, int *n, int *k, d *ab, int *ldab, d *work) nogil
cdef d _wrap_dlansb(char *norm, char *uplo, int *n, int *k, d *ab, int *ldab, d *work) nogil:
    cdef d out
    _fortran_dlansb(&out, norm, uplo, n, k, ab, ldab, work)
    return out
cdef dlansb_t *dlansb_f = &_wrap_dlansb

ctypedef void _wdlansp_t(d *out, char *norm, char *uplo, int *n, d *ap, d *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dlansp "F_FUNC(dlanspwrp, DLANSPWRP)"(d *out, char *norm, char *uplo, int *n, d *ap, d *work) nogil
cdef d _wrap_dlansp(char *norm, char *uplo, int *n, d *ap, d *work) nogil:
    cdef d out
    _fortran_dlansp(&out, norm, uplo, n, ap, work)
    return out
cdef dlansp_t *dlansp_f = &_wrap_dlansp

ctypedef void _wdlanst_t(d *out, char *norm, int *n, d *d_, d *e) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dlanst "F_FUNC(dlanstwrp, DLANSTWRP)"(d *out, char *norm, int *n, d *d_, d *e) nogil
cdef d _wrap_dlanst(char *norm, int *n, d *d_, d *e) nogil:
    cdef d out
    _fortran_dlanst(&out, norm, n, d_, e)
    return out
cdef dlanst_t *dlanst_f = &_wrap_dlanst

ctypedef void _wdlansy_t(d *out, char *norm, char *uplo, int *n, d *a, int *lda, d *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dlansy "F_FUNC(dlansywrp, DLANSYWRP)"(d *out, char *norm, char *uplo, int *n, d *a, int *lda, d *work) nogil
cdef d _wrap_dlansy(char *norm, char *uplo, int *n, d *a, int *lda, d *work) nogil:
    cdef d out
    _fortran_dlansy(&out, norm, uplo, n, a, lda, work)
    return out
cdef dlansy_t *dlansy_f = &_wrap_dlansy

ctypedef void _wdlantb_t(d *out, char *norm, char *uplo, char *diag, int *n, int *k, d *ab, int *ldab, d *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dlantb "F_FUNC(dlantbwrp, DLANTBWRP)"(d *out, char *norm, char *uplo, char *diag, int *n, int *k, d *ab, int *ldab, d *work) nogil
cdef d _wrap_dlantb(char *norm, char *uplo, char *diag, int *n, int *k, d *ab, int *ldab, d *work) nogil:
    cdef d out
    _fortran_dlantb(&out, norm, uplo, diag, n, k, ab, ldab, work)
    return out
cdef dlantb_t *dlantb_f = &_wrap_dlantb

ctypedef void _wdlantp_t(d *out, char *norm, char *uplo, char *diag, int *n, d *ap, d *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dlantp "F_FUNC(dlantpwrp, DLANTPWRP)"(d *out, char *norm, char *uplo, char *diag, int *n, d *ap, d *work) nogil
cdef d _wrap_dlantp(char *norm, char *uplo, char *diag, int *n, d *ap, d *work) nogil:
    cdef d out
    _fortran_dlantp(&out, norm, uplo, diag, n, ap, work)
    return out
cdef dlantp_t *dlantp_f = &_wrap_dlantp

ctypedef void _wdlantr_t(d *out, char *norm, char *uplo, char *diag, int *m, int *n, d *a, int *lda, d *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dlantr "F_FUNC(dlantrwrp, DLANTRWRP)"(d *out, char *norm, char *uplo, char *diag, int *m, int *n, d *a, int *lda, d *work) nogil
cdef d _wrap_dlantr(char *norm, char *uplo, char *diag, int *m, int *n, d *a, int *lda, d *work) nogil:
    cdef d out
    _fortran_dlantr(&out, norm, uplo, diag, m, n, a, lda, work)
    return out
cdef dlantr_t *dlantr_f = &_wrap_dlantr

ctypedef void _wdzsum1_t(d *out, int *n, z *cx, int *incx) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dzsum1 "F_FUNC(dzsum1wrp, DZSUM1WRP)"(d *out, int *n, z *cx, int *incx) nogil
cdef d _wrap_dzsum1(int *n, z *cx, int *incx) nogil:
    cdef d out
    _fortran_dzsum1(&out, n, cx, incx)
    return out
cdef dzsum1_t *dzsum1_f = &_wrap_dzsum1

ctypedef void _wicmax1_t(int *out, int *n, c *cx, int *incx) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_icmax1 "F_FUNC(icmax1wrp, ICMAX1WRP)"(int *out, int *n, c *cx, int *incx) nogil
cdef int _wrap_icmax1(int *n, c *cx, int *incx) nogil:
    cdef int out
    _fortran_icmax1(&out, n, cx, incx)
    return out
cdef icmax1_t *icmax1_f = &_wrap_icmax1

ctypedef void _wieeeck_t(int *out, int *ispec, s *zero, s *one) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ieeeck "F_FUNC(ieeeckwrp, IEEECKWRP)"(int *out, int *ispec, s *zero, s *one) nogil
cdef int _wrap_ieeeck(int *ispec, s *zero, s *one) nogil:
    cdef int out
    _fortran_ieeeck(&out, ispec, zero, one)
    return out
cdef ieeeck_t *ieeeck_f = &_wrap_ieeeck

ctypedef void _wilaenv_t(int *out, int *ispec, char *name, char *opts, int *n1, int *n2, int *n3, int *n4) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ilaenv "F_FUNC(ilaenvwrp, ILAENVWRP)"(int *out, int *ispec, char *name, char *opts, int *n1, int *n2, int *n3, int *n4) nogil
cdef int _wrap_ilaenv(int *ispec, char *name, char *opts, int *n1, int *n2, int *n3, int *n4) nogil:
    cdef int out
    _fortran_ilaenv(&out, ispec, name, opts, n1, n2, n3, n4)
    return out
cdef ilaenv_t *ilaenv_f = &_wrap_ilaenv

ctypedef void _wiparmq_t(int *out, int *ispec, char *name, char *opts, int *n, int *ilo, int *ihi, int *lwork) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_iparmq "F_FUNC(iparmqwrp, IPARMQWRP)"(int *out, int *ispec, char *name, char *opts, int *n, int *ilo, int *ihi, int *lwork) nogil
cdef int _wrap_iparmq(int *ispec, char *name, char *opts, int *n, int *ilo, int *ihi, int *lwork) nogil:
    cdef int out
    _fortran_iparmq(&out, ispec, name, opts, n, ilo, ihi, lwork)
    return out
cdef iparmq_t *iparmq_f = &_wrap_iparmq

ctypedef void _wizmax1_t(int *out, int *n, z *cx, int *incx) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_izmax1 "F_FUNC(izmax1wrp, IZMAX1WRP)"(int *out, int *n, z *cx, int *incx) nogil
cdef int _wrap_izmax1(int *n, z *cx, int *incx) nogil:
    cdef int out
    _fortran_izmax1(&out, n, cx, incx)
    return out
cdef izmax1_t *izmax1_f = &_wrap_izmax1

ctypedef void _wlsamen_t(bint *out, int *n, char *ca, char *cb) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_lsamen "F_FUNC(lsamenwrp, LSAMENWRP)"(bint *out, int *n, char *ca, char *cb) nogil
cdef bint _wrap_lsamen(int *n, char *ca, char *cb) nogil:
    cdef bint out
    _fortran_lsamen(&out, n, ca, cb)
    return out
cdef lsamen_t *lsamen_f = &_wrap_lsamen

ctypedef void _wscsum1_t(s *out, int *n, c *cx, int *incx) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_scsum1 "F_FUNC(scsum1wrp, SCSUM1WRP)"(s *out, int *n, c *cx, int *incx) nogil
cdef s _wrap_scsum1(int *n, c *cx, int *incx) nogil:
    cdef s out
    _fortran_scsum1(&out, n, cx, incx)
    return out
cdef scsum1_t *scsum1_f = &_wrap_scsum1

ctypedef void _wslamch_t(s *out, char *cmach) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_slamch "F_FUNC(slamchwrp, SLAMCHWRP)"(s *out, char *cmach) nogil
cdef s _wrap_slamch(char *cmach) nogil:
    cdef s out
    _fortran_slamch(&out, cmach)
    return out
cdef slamch_t *slamch_f = &_wrap_slamch

ctypedef void _wslangb_t(s *out, char *norm, int *n, int *kl, int *ku, s *ab, int *ldab, s *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_slangb "F_FUNC(slangbwrp, SLANGBWRP)"(s *out, char *norm, int *n, int *kl, int *ku, s *ab, int *ldab, s *work) nogil
cdef s _wrap_slangb(char *norm, int *n, int *kl, int *ku, s *ab, int *ldab, s *work) nogil:
    cdef s out
    _fortran_slangb(&out, norm, n, kl, ku, ab, ldab, work)
    return out
cdef slangb_t *slangb_f = &_wrap_slangb

ctypedef void _wslange_t(s *out, char *norm, int *m, int *n, s *a, int *lda, s *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_slange "F_FUNC(slangewrp, SLANGEWRP)"(s *out, char *norm, int *m, int *n, s *a, int *lda, s *work) nogil
cdef s _wrap_slange(char *norm, int *m, int *n, s *a, int *lda, s *work) nogil:
    cdef s out
    _fortran_slange(&out, norm, m, n, a, lda, work)
    return out
cdef slange_t *slange_f = &_wrap_slange

ctypedef void _wslangt_t(s *out, char *norm, int *n, s *dl, s *d, s *du) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_slangt "F_FUNC(slangtwrp, SLANGTWRP)"(s *out, char *norm, int *n, s *dl, s *d, s *du) nogil
cdef s _wrap_slangt(char *norm, int *n, s *dl, s *d, s *du) nogil:
    cdef s out
    _fortran_slangt(&out, norm, n, dl, d, du)
    return out
cdef slangt_t *slangt_f = &_wrap_slangt

ctypedef void _wslanhs_t(s *out, char *norm, int *n, s *a, int *lda, s *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_slanhs "F_FUNC(slanhswrp, SLANHSWRP)"(s *out, char *norm, int *n, s *a, int *lda, s *work) nogil
cdef s _wrap_slanhs(char *norm, int *n, s *a, int *lda, s *work) nogil:
    cdef s out
    _fortran_slanhs(&out, norm, n, a, lda, work)
    return out
cdef slanhs_t *slanhs_f = &_wrap_slanhs

ctypedef void _wslansb_t(s *out, char *norm, char *uplo, int *n, int *k, s *ab, int *ldab, s *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_slansb "F_FUNC(slansbwrp, SLANSBWRP)"(s *out, char *norm, char *uplo, int *n, int *k, s *ab, int *ldab, s *work) nogil
cdef s _wrap_slansb(char *norm, char *uplo, int *n, int *k, s *ab, int *ldab, s *work) nogil:
    cdef s out
    _fortran_slansb(&out, norm, uplo, n, k, ab, ldab, work)
    return out
cdef slansb_t *slansb_f = &_wrap_slansb

ctypedef void _wslansp_t(s *out, char *norm, char *uplo, int *n, s *ap, s *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_slansp "F_FUNC(slanspwrp, SLANSPWRP)"(s *out, char *norm, char *uplo, int *n, s *ap, s *work) nogil
cdef s _wrap_slansp(char *norm, char *uplo, int *n, s *ap, s *work) nogil:
    cdef s out
    _fortran_slansp(&out, norm, uplo, n, ap, work)
    return out
cdef slansp_t *slansp_f = &_wrap_slansp

ctypedef void _wslanst_t(s *out, char *norm, int *n, s *d, s *e) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_slanst "F_FUNC(slanstwrp, SLANSTWRP)"(s *out, char *norm, int *n, s *d, s *e) nogil
cdef s _wrap_slanst(char *norm, int *n, s *d, s *e) nogil:
    cdef s out
    _fortran_slanst(&out, norm, n, d, e)
    return out
cdef slanst_t *slanst_f = &_wrap_slanst

ctypedef void _wslansy_t(s *out, char *norm, char *uplo, int *n, s *a, int *lda, s *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_slansy "F_FUNC(slansywrp, SLANSYWRP)"(s *out, char *norm, char *uplo, int *n, s *a, int *lda, s *work) nogil
cdef s _wrap_slansy(char *norm, char *uplo, int *n, s *a, int *lda, s *work) nogil:
    cdef s out
    _fortran_slansy(&out, norm, uplo, n, a, lda, work)
    return out
cdef slansy_t *slansy_f = &_wrap_slansy

ctypedef void _wslantb_t(s *out, char *norm, char *uplo, char *diag, int *n, int *k, s *ab, int *ldab, s *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_slantb "F_FUNC(slantbwrp, SLANTBWRP)"(s *out, char *norm, char *uplo, char *diag, int *n, int *k, s *ab, int *ldab, s *work) nogil
cdef s _wrap_slantb(char *norm, char *uplo, char *diag, int *n, int *k, s *ab, int *ldab, s *work) nogil:
    cdef s out
    _fortran_slantb(&out, norm, uplo, diag, n, k, ab, ldab, work)
    return out
cdef slantb_t *slantb_f = &_wrap_slantb

ctypedef void _wslantp_t(s *out, char *norm, char *uplo, char *diag, int *n, s *ap, s *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_slantp "F_FUNC(slantpwrp, SLANTPWRP)"(s *out, char *norm, char *uplo, char *diag, int *n, s *ap, s *work) nogil
cdef s _wrap_slantp(char *norm, char *uplo, char *diag, int *n, s *ap, s *work) nogil:
    cdef s out
    _fortran_slantp(&out, norm, uplo, diag, n, ap, work)
    return out
cdef slantp_t *slantp_f = &_wrap_slantp

ctypedef void _wslantr_t(s *out, char *norm, char *uplo, char *diag, int *m, int *n, s *a, int *lda, s *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_slantr "F_FUNC(slantrwrp, SLANTRWRP)"(s *out, char *norm, char *uplo, char *diag, int *m, int *n, s *a, int *lda, s *work) nogil
cdef s _wrap_slantr(char *norm, char *uplo, char *diag, int *m, int *n, s *a, int *lda, s *work) nogil:
    cdef s out
    _fortran_slantr(&out, norm, uplo, diag, m, n, a, lda, work)
    return out
cdef slantr_t *slantr_f = &_wrap_slantr

ctypedef void _wzlangb_t(d *out, char *norm, int *n, int *kl, int *ku, z *ab, int *ldab, d *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zlangb "F_FUNC(zlangbwrp, ZLANGBWRP)"(d *out, char *norm, int *n, int *kl, int *ku, z *ab, int *ldab, d *work) nogil
cdef d _wrap_zlangb(char *norm, int *n, int *kl, int *ku, z *ab, int *ldab, d *work) nogil:
    cdef d out
    _fortran_zlangb(&out, norm, n, kl, ku, ab, ldab, work)
    return out
cdef zlangb_t *zlangb_f = &_wrap_zlangb

ctypedef void _wzlange_t(d *out, char *norm, int *m, int *n, z *a, int *lda, d *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zlange "F_FUNC(zlangewrp, ZLANGEWRP)"(d *out, char *norm, int *m, int *n, z *a, int *lda, d *work) nogil
cdef d _wrap_zlange(char *norm, int *m, int *n, z *a, int *lda, d *work) nogil:
    cdef d out
    _fortran_zlange(&out, norm, m, n, a, lda, work)
    return out
cdef zlange_t *zlange_f = &_wrap_zlange

ctypedef void _wzlangt_t(d *out, char *norm, int *n, z *dl, z *d_, z *du) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zlangt "F_FUNC(zlangtwrp, ZLANGTWRP)"(d *out, char *norm, int *n, z *dl, z *d_, z *du) nogil
cdef d _wrap_zlangt(char *norm, int *n, z *dl, z *d_, z *du) nogil:
    cdef d out
    _fortran_zlangt(&out, norm, n, dl, d_, du)
    return out
cdef zlangt_t *zlangt_f = &_wrap_zlangt

ctypedef void _wzlanhb_t(d *out, char *norm, char *uplo, int *n, int *k, z *ab, int *ldab, d *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zlanhb "F_FUNC(zlanhbwrp, ZLANHBWRP)"(d *out, char *norm, char *uplo, int *n, int *k, z *ab, int *ldab, d *work) nogil
cdef d _wrap_zlanhb(char *norm, char *uplo, int *n, int *k, z *ab, int *ldab, d *work) nogil:
    cdef d out
    _fortran_zlanhb(&out, norm, uplo, n, k, ab, ldab, work)
    return out
cdef zlanhb_t *zlanhb_f = &_wrap_zlanhb

ctypedef void _wzlanhe_t(d *out, char *norm, char *uplo, int *n, z *a, int *lda, d *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zlanhe "F_FUNC(zlanhewrp, ZLANHEWRP)"(d *out, char *norm, char *uplo, int *n, z *a, int *lda, d *work) nogil
cdef d _wrap_zlanhe(char *norm, char *uplo, int *n, z *a, int *lda, d *work) nogil:
    cdef d out
    _fortran_zlanhe(&out, norm, uplo, n, a, lda, work)
    return out
cdef zlanhe_t *zlanhe_f = &_wrap_zlanhe

ctypedef void _wzlanhp_t(d *out, char *norm, char *uplo, int *n, z *ap, d *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zlanhp "F_FUNC(zlanhpwrp, ZLANHPWRP)"(d *out, char *norm, char *uplo, int *n, z *ap, d *work) nogil
cdef d _wrap_zlanhp(char *norm, char *uplo, int *n, z *ap, d *work) nogil:
    cdef d out
    _fortran_zlanhp(&out, norm, uplo, n, ap, work)
    return out
cdef zlanhp_t *zlanhp_f = &_wrap_zlanhp

ctypedef void _wzlanhs_t(d *out, char *norm, int *n, z *a, int *lda, d *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zlanhs "F_FUNC(zlanhswrp, ZLANHSWRP)"(d *out, char *norm, int *n, z *a, int *lda, d *work) nogil
cdef d _wrap_zlanhs(char *norm, int *n, z *a, int *lda, d *work) nogil:
    cdef d out
    _fortran_zlanhs(&out, norm, n, a, lda, work)
    return out
cdef zlanhs_t *zlanhs_f = &_wrap_zlanhs

ctypedef void _wzlanht_t(d *out, char *norm, int *n, d *d_, z *e) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zlanht "F_FUNC(zlanhtwrp, ZLANHTWRP)"(d *out, char *norm, int *n, d *d_, z *e) nogil
cdef d _wrap_zlanht(char *norm, int *n, d *d_, z *e) nogil:
    cdef d out
    _fortran_zlanht(&out, norm, n, d_, e)
    return out
cdef zlanht_t *zlanht_f = &_wrap_zlanht

ctypedef void _wzlansb_t(d *out, char *norm, char *uplo, int *n, int *k, z *ab, int *ldab, d *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zlansb "F_FUNC(zlansbwrp, ZLANSBWRP)"(d *out, char *norm, char *uplo, int *n, int *k, z *ab, int *ldab, d *work) nogil
cdef d _wrap_zlansb(char *norm, char *uplo, int *n, int *k, z *ab, int *ldab, d *work) nogil:
    cdef d out
    _fortran_zlansb(&out, norm, uplo, n, k, ab, ldab, work)
    return out
cdef zlansb_t *zlansb_f = &_wrap_zlansb

ctypedef void _wzlansp_t(d *out, char *norm, char *uplo, int *n, z *ap, d *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zlansp "F_FUNC(zlanspwrp, ZLANSPWRP)"(d *out, char *norm, char *uplo, int *n, z *ap, d *work) nogil
cdef d _wrap_zlansp(char *norm, char *uplo, int *n, z *ap, d *work) nogil:
    cdef d out
    _fortran_zlansp(&out, norm, uplo, n, ap, work)
    return out
cdef zlansp_t *zlansp_f = &_wrap_zlansp

ctypedef void _wzlansy_t(d *out, char *norm, char *uplo, int *n, z *a, int *lda, d *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zlansy "F_FUNC(zlansywrp, ZLANSYWRP)"(d *out, char *norm, char *uplo, int *n, z *a, int *lda, d *work) nogil
cdef d _wrap_zlansy(char *norm, char *uplo, int *n, z *a, int *lda, d *work) nogil:
    cdef d out
    _fortran_zlansy(&out, norm, uplo, n, a, lda, work)
    return out
cdef zlansy_t *zlansy_f = &_wrap_zlansy

ctypedef void _wzlantb_t(d *out, char *norm, char *uplo, char *diag, int *n, int *k, z *ab, int *ldab, d *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zlantb "F_FUNC(zlantbwrp, ZLANTBWRP)"(d *out, char *norm, char *uplo, char *diag, int *n, int *k, z *ab, int *ldab, d *work) nogil
cdef d _wrap_zlantb(char *norm, char *uplo, char *diag, int *n, int *k, z *ab, int *ldab, d *work) nogil:
    cdef d out
    _fortran_zlantb(&out, norm, uplo, diag, n, k, ab, ldab, work)
    return out
cdef zlantb_t *zlantb_f = &_wrap_zlantb

ctypedef void _wzlantp_t(d *out, char *norm, char *uplo, char *diag, int *n, z *ap, d *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zlantp "F_FUNC(zlantpwrp, ZLANTPWRP)"(d *out, char *norm, char *uplo, char *diag, int *n, z *ap, d *work) nogil
cdef d _wrap_zlantp(char *norm, char *uplo, char *diag, int *n, z *ap, d *work) nogil:
    cdef d out
    _fortran_zlantp(&out, norm, uplo, diag, n, ap, work)
    return out
cdef zlantp_t *zlantp_f = &_wrap_zlantp

ctypedef void _wzlantr_t(d *out, char *norm, char *uplo, char *diag, int *m, int *n, z *a, int *lda, d *work) nogil
cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zlantr "F_FUNC(zlantrwrp, ZLANTRWRP)"(d *out, char *norm, char *uplo, char *diag, int *m, int *n, z *a, int *lda, d *work) nogil
cdef d _wrap_zlantr(char *norm, char *uplo, char *diag, int *m, int *n, z *a, int *lda, d *work) nogil:
    cdef d out
    _fortran_zlantr(&out, norm, uplo, diag, m, n, a, lda, work)
    return out
cdef zlantr_t *zlantr_f = &_wrap_zlantr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cbdsqr "F_FUNC(cbdsqr,CBDSQR)"(char *uplo, int *n, int *ncvt, int *nru, int *ncc, s *d, s *e, c *vt, int *ldvt, c *u, int *ldu, c *c, int *ldc, s *rwork, int *info) nogil
cdef cbdsqr_t *cbdsqr_f = &_fortran_cbdsqr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgbbrd "F_FUNC(cgbbrd,CGBBRD)"(char *vect, int *m, int *n, int *ncc, int *kl, int *ku, c *ab, int *ldab, s *d, s *e, c *q, int *ldq, c *pt, int *ldpt, c *c, int *ldc, c *work, s *rwork, int *info) nogil
cdef cgbbrd_t *cgbbrd_f = &_fortran_cgbbrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgbcon "F_FUNC(cgbcon,CGBCON)"(char *norm, int *n, int *kl, int *ku, c *ab, int *ldab, int *ipiv, s *anorm, s *rcond, c *work, s *rwork, int *info) nogil
cdef cgbcon_t *cgbcon_f = &_fortran_cgbcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgbequ "F_FUNC(cgbequ,CGBEQU)"(int *m, int *n, int *kl, int *ku, c *ab, int *ldab, s *r, s *c, s *rowcnd, s *colcnd, s *amax, int *info) nogil
cdef cgbequ_t *cgbequ_f = &_fortran_cgbequ

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgbrfs "F_FUNC(cgbrfs,CGBRFS)"(char *trans, int *n, int *kl, int *ku, int *nrhs, c *ab, int *ldab, c *afb, int *ldafb, int *ipiv, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cgbrfs_t *cgbrfs_f = &_fortran_cgbrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgbsv "F_FUNC(cgbsv,CGBSV)"(int *n, int *kl, int *ku, int *nrhs, c *ab, int *ldab, int *ipiv, c *b, int *ldb, int *info) nogil
cdef cgbsv_t *cgbsv_f = &_fortran_cgbsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgbsvx "F_FUNC(cgbsvx,CGBSVX)"(char *fact, char *trans, int *n, int *kl, int *ku, int *nrhs, c *ab, int *ldab, c *afb, int *ldafb, int *ipiv, char *equed, s *r, s *c, c *b, int *ldb, c *x, int *ldx, s *rcond, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cgbsvx_t *cgbsvx_f = &_fortran_cgbsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgbtf2 "F_FUNC(cgbtf2,CGBTF2)"(int *m, int *n, int *kl, int *ku, c *ab, int *ldab, int *ipiv, int *info) nogil
cdef cgbtf2_t *cgbtf2_f = &_fortran_cgbtf2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgbtrf "F_FUNC(cgbtrf,CGBTRF)"(int *m, int *n, int *kl, int *ku, c *ab, int *ldab, int *ipiv, int *info) nogil
cdef cgbtrf_t *cgbtrf_f = &_fortran_cgbtrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgbtrs "F_FUNC(cgbtrs,CGBTRS)"(char *trans, int *n, int *kl, int *ku, int *nrhs, c *ab, int *ldab, int *ipiv, c *b, int *ldb, int *info) nogil
cdef cgbtrs_t *cgbtrs_f = &_fortran_cgbtrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgebak "F_FUNC(cgebak,CGEBAK)"(char *job, char *side, int *n, int *ilo, int *ihi, s *scale, int *m, c *v, int *ldv, int *info) nogil
cdef cgebak_t *cgebak_f = &_fortran_cgebak

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgebal "F_FUNC(cgebal,CGEBAL)"(char *job, int *n, c *a, int *lda, int *ilo, int *ihi, s *scale, int *info) nogil
cdef cgebal_t *cgebal_f = &_fortran_cgebal

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgebd2 "F_FUNC(cgebd2,CGEBD2)"(int *m, int *n, c *a, int *lda, s *d, s *e, c *tauq, c *taup, c *work, int *info) nogil
cdef cgebd2_t *cgebd2_f = &_fortran_cgebd2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgebrd "F_FUNC(cgebrd,CGEBRD)"(int *m, int *n, c *a, int *lda, s *d, s *e, c *tauq, c *taup, c *work, int *lwork, int *info) nogil
cdef cgebrd_t *cgebrd_f = &_fortran_cgebrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgecon "F_FUNC(cgecon,CGECON)"(char *norm, int *n, c *a, int *lda, s *anorm, s *rcond, c *work, s *rwork, int *info) nogil
cdef cgecon_t *cgecon_f = &_fortran_cgecon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgeequ "F_FUNC(cgeequ,CGEEQU)"(int *m, int *n, c *a, int *lda, s *r, s *c, s *rowcnd, s *colcnd, s *amax, int *info) nogil
cdef cgeequ_t *cgeequ_f = &_fortran_cgeequ

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgees "F_FUNC(cgees,CGEES)"(char *jobvs, char *sort, cselect1 *select, int *n, c *a, int *lda, int *sdim, c *w, c *vs, int *ldvs, c *work, int *lwork, s *rwork, bint *bwork, int *info) nogil
cdef cgees_t *cgees_f = &_fortran_cgees

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgeesx "F_FUNC(cgeesx,CGEESX)"(char *jobvs, char *sort, cselect1 *select, char *sense, int *n, c *a, int *lda, int *sdim, c *w, c *vs, int *ldvs, s *rconde, s *rcondv, c *work, int *lwork, s *rwork, bint *bwork, int *info) nogil
cdef cgeesx_t *cgeesx_f = &_fortran_cgeesx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgeev "F_FUNC(cgeev,CGEEV)"(char *jobvl, char *jobvr, int *n, c *a, int *lda, c *w, c *vl, int *ldvl, c *vr, int *ldvr, c *work, int *lwork, s *rwork, int *info) nogil
cdef cgeev_t *cgeev_f = &_fortran_cgeev

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgeevx "F_FUNC(cgeevx,CGEEVX)"(char *balanc, char *jobvl, char *jobvr, char *sense, int *n, c *a, int *lda, c *w, c *vl, int *ldvl, c *vr, int *ldvr, int *ilo, int *ihi, s *scale, s *abnrm, s *rconde, s *rcondv, c *work, int *lwork, s *rwork, int *info) nogil
cdef cgeevx_t *cgeevx_f = &_fortran_cgeevx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgegs "F_FUNC(cgegs,CGEGS)"(char *jobvsl, char *jobvsr, int *n, c *a, int *lda, c *b, int *ldb, c *alpha, c *beta, c *vsl, int *ldvsl, c *vsr, int *ldvsr, c *work, int *lwork, s *rwork, int *info) nogil
cdef cgegs_t *cgegs_f = &_fortran_cgegs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgegv "F_FUNC(cgegv,CGEGV)"(char *jobvl, char *jobvr, int *n, c *a, int *lda, c *b, int *ldb, c *alpha, c *beta, c *vl, int *ldvl, c *vr, int *ldvr, c *work, int *lwork, s *rwork, int *info) nogil
cdef cgegv_t *cgegv_f = &_fortran_cgegv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgehd2 "F_FUNC(cgehd2,CGEHD2)"(int *n, int *ilo, int *ihi, c *a, int *lda, c *tau, c *work, int *info) nogil
cdef cgehd2_t *cgehd2_f = &_fortran_cgehd2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgehrd "F_FUNC(cgehrd,CGEHRD)"(int *n, int *ilo, int *ihi, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef cgehrd_t *cgehrd_f = &_fortran_cgehrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgelq2 "F_FUNC(cgelq2,CGELQ2)"(int *m, int *n, c *a, int *lda, c *tau, c *work, int *info) nogil
cdef cgelq2_t *cgelq2_f = &_fortran_cgelq2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgelqf "F_FUNC(cgelqf,CGELQF)"(int *m, int *n, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef cgelqf_t *cgelqf_f = &_fortran_cgelqf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgels "F_FUNC(cgels,CGELS)"(char *trans, int *m, int *n, int *nrhs, c *a, int *lda, c *b, int *ldb, c *work, int *lwork, int *info) nogil
cdef cgels_t *cgels_f = &_fortran_cgels

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgelsd "F_FUNC(cgelsd,CGELSD)"(int *m, int *n, int *nrhs, c *a, int *lda, c *b, int *ldb, s *s, s *rcond, int *rank, c *work, int *lwork, s *rwork, int *iwork, int *info) nogil
cdef cgelsd_t *cgelsd_f = &_fortran_cgelsd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgelss "F_FUNC(cgelss,CGELSS)"(int *m, int *n, int *nrhs, c *a, int *lda, c *b, int *ldb, s *s, s *rcond, int *rank, c *work, int *lwork, s *rwork, int *info) nogil
cdef cgelss_t *cgelss_f = &_fortran_cgelss

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgelsx "F_FUNC(cgelsx,CGELSX)"(int *m, int *n, int *nrhs, c *a, int *lda, c *b, int *ldb, int *jpvt, s *rcond, int *rank, c *work, s *rwork, int *info) nogil
cdef cgelsx_t *cgelsx_f = &_fortran_cgelsx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgelsy "F_FUNC(cgelsy,CGELSY)"(int *m, int *n, int *nrhs, c *a, int *lda, c *b, int *ldb, int *jpvt, s *rcond, int *rank, c *work, int *lwork, s *rwork, int *info) nogil
cdef cgelsy_t *cgelsy_f = &_fortran_cgelsy

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgeql2 "F_FUNC(cgeql2,CGEQL2)"(int *m, int *n, c *a, int *lda, c *tau, c *work, int *info) nogil
cdef cgeql2_t *cgeql2_f = &_fortran_cgeql2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgeqlf "F_FUNC(cgeqlf,CGEQLF)"(int *m, int *n, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef cgeqlf_t *cgeqlf_f = &_fortran_cgeqlf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgeqp3 "F_FUNC(cgeqp3,CGEQP3)"(int *m, int *n, c *a, int *lda, int *jpvt, c *tau, c *work, int *lwork, s *rwork, int *info) nogil
cdef cgeqp3_t *cgeqp3_f = &_fortran_cgeqp3

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgeqpf "F_FUNC(cgeqpf,CGEQPF)"(int *m, int *n, c *a, int *lda, int *jpvt, c *tau, c *work, s *rwork, int *info) nogil
cdef cgeqpf_t *cgeqpf_f = &_fortran_cgeqpf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgeqr2 "F_FUNC(cgeqr2,CGEQR2)"(int *m, int *n, c *a, int *lda, c *tau, c *work, int *info) nogil
cdef cgeqr2_t *cgeqr2_f = &_fortran_cgeqr2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgeqrf "F_FUNC(cgeqrf,CGEQRF)"(int *m, int *n, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef cgeqrf_t *cgeqrf_f = &_fortran_cgeqrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgerfs "F_FUNC(cgerfs,CGERFS)"(char *trans, int *n, int *nrhs, c *a, int *lda, c *af, int *ldaf, int *ipiv, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cgerfs_t *cgerfs_f = &_fortran_cgerfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgerq2 "F_FUNC(cgerq2,CGERQ2)"(int *m, int *n, c *a, int *lda, c *tau, c *work, int *info) nogil
cdef cgerq2_t *cgerq2_f = &_fortran_cgerq2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgerqf "F_FUNC(cgerqf,CGERQF)"(int *m, int *n, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef cgerqf_t *cgerqf_f = &_fortran_cgerqf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgesc2 "F_FUNC(cgesc2,CGESC2)"(int *n, c *a, int *lda, c *rhs, int *ipiv, int *jpiv, s *scale) nogil
cdef cgesc2_t *cgesc2_f = &_fortran_cgesc2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgesdd "F_FUNC(cgesdd,CGESDD)"(char *jobz, int *m, int *n, c *a, int *lda, s *s, c *u, int *ldu, c *vt, int *ldvt, c *work, int *lwork, s *rwork, int *iwork, int *info) nogil
cdef cgesdd_t *cgesdd_f = &_fortran_cgesdd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgesv "F_FUNC(cgesv,CGESV)"(int *n, int *nrhs, c *a, int *lda, int *ipiv, c *b, int *ldb, int *info) nogil
cdef cgesv_t *cgesv_f = &_fortran_cgesv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgesvd "F_FUNC(cgesvd,CGESVD)"(char *jobu, char *jobvt, int *m, int *n, c *a, int *lda, s *s, c *u, int *ldu, c *vt, int *ldvt, c *work, int *lwork, s *rwork, int *info) nogil
cdef cgesvd_t *cgesvd_f = &_fortran_cgesvd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgesvx "F_FUNC(cgesvx,CGESVX)"(char *fact, char *trans, int *n, int *nrhs, c *a, int *lda, c *af, int *ldaf, int *ipiv, char *equed, s *r, s *c, c *b, int *ldb, c *x, int *ldx, s *rcond, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cgesvx_t *cgesvx_f = &_fortran_cgesvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgetc2 "F_FUNC(cgetc2,CGETC2)"(int *n, c *a, int *lda, int *ipiv, int *jpiv, int *info) nogil
cdef cgetc2_t *cgetc2_f = &_fortran_cgetc2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgetf2 "F_FUNC(cgetf2,CGETF2)"(int *m, int *n, c *a, int *lda, int *ipiv, int *info) nogil
cdef cgetf2_t *cgetf2_f = &_fortran_cgetf2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgetrf "F_FUNC(cgetrf,CGETRF)"(int *m, int *n, c *a, int *lda, int *ipiv, int *info) nogil
cdef cgetrf_t *cgetrf_f = &_fortran_cgetrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgetri "F_FUNC(cgetri,CGETRI)"(int *n, c *a, int *lda, int *ipiv, c *work, int *lwork, int *info) nogil
cdef cgetri_t *cgetri_f = &_fortran_cgetri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgetrs "F_FUNC(cgetrs,CGETRS)"(char *trans, int *n, int *nrhs, c *a, int *lda, int *ipiv, c *b, int *ldb, int *info) nogil
cdef cgetrs_t *cgetrs_f = &_fortran_cgetrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cggbak "F_FUNC(cggbak,CGGBAK)"(char *job, char *side, int *n, int *ilo, int *ihi, s *lscale, s *rscale, int *m, c *v, int *ldv, int *info) nogil
cdef cggbak_t *cggbak_f = &_fortran_cggbak

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cggbal "F_FUNC(cggbal,CGGBAL)"(char *job, int *n, c *a, int *lda, c *b, int *ldb, int *ilo, int *ihi, s *lscale, s *rscale, s *work, int *info) nogil
cdef cggbal_t *cggbal_f = &_fortran_cggbal

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgges "F_FUNC(cgges,CGGES)"(char *jobvsl, char *jobvsr, char *sort, cselect2 *selctg, int *n, c *a, int *lda, c *b, int *ldb, int *sdim, c *alpha, c *beta, c *vsl, int *ldvsl, c *vsr, int *ldvsr, c *work, int *lwork, s *rwork, bint *bwork, int *info) nogil
cdef cgges_t *cgges_f = &_fortran_cgges

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cggesx "F_FUNC(cggesx,CGGESX)"(char *jobvsl, char *jobvsr, char *sort, cselect2 *selctg, char *sense, int *n, c *a, int *lda, c *b, int *ldb, int *sdim, c *alpha, c *beta, c *vsl, int *ldvsl, c *vsr, int *ldvsr, s *rconde, s *rcondv, c *work, int *lwork, s *rwork, int *iwork, int *liwork, bint *bwork, int *info) nogil
cdef cggesx_t *cggesx_f = &_fortran_cggesx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cggev "F_FUNC(cggev,CGGEV)"(char *jobvl, char *jobvr, int *n, c *a, int *lda, c *b, int *ldb, c *alpha, c *beta, c *vl, int *ldvl, c *vr, int *ldvr, c *work, int *lwork, s *rwork, int *info) nogil
cdef cggev_t *cggev_f = &_fortran_cggev

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cggevx "F_FUNC(cggevx,CGGEVX)"(char *balanc, char *jobvl, char *jobvr, char *sense, int *n, c *a, int *lda, c *b, int *ldb, c *alpha, c *beta, c *vl, int *ldvl, c *vr, int *ldvr, int *ilo, int *ihi, s *lscale, s *rscale, s *abnrm, s *bbnrm, s *rconde, s *rcondv, c *work, int *lwork, s *rwork, int *iwork, bint *bwork, int *info) nogil
cdef cggevx_t *cggevx_f = &_fortran_cggevx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cggglm "F_FUNC(cggglm,CGGGLM)"(int *n, int *m, int *p, c *a, int *lda, c *b, int *ldb, c *d, c *x, c *y, c *work, int *lwork, int *info) nogil
cdef cggglm_t *cggglm_f = &_fortran_cggglm

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgghrd "F_FUNC(cgghrd,CGGHRD)"(char *compq, char *compz, int *n, int *ilo, int *ihi, c *a, int *lda, c *b, int *ldb, c *q, int *ldq, c *z, int *ldz, int *info) nogil
cdef cgghrd_t *cgghrd_f = &_fortran_cgghrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgglse "F_FUNC(cgglse,CGGLSE)"(int *m, int *n, int *p, c *a, int *lda, c *b, int *ldb, c *c, c *d, c *x, c *work, int *lwork, int *info) nogil
cdef cgglse_t *cgglse_f = &_fortran_cgglse

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cggqrf "F_FUNC(cggqrf,CGGQRF)"(int *n, int *m, int *p, c *a, int *lda, c *taua, c *b, int *ldb, c *taub, c *work, int *lwork, int *info) nogil
cdef cggqrf_t *cggqrf_f = &_fortran_cggqrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cggrqf "F_FUNC(cggrqf,CGGRQF)"(int *m, int *p, int *n, c *a, int *lda, c *taua, c *b, int *ldb, c *taub, c *work, int *lwork, int *info) nogil
cdef cggrqf_t *cggrqf_f = &_fortran_cggrqf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cggsvd "F_FUNC(cggsvd,CGGSVD)"(char *jobu, char *jobv, char *jobq, int *m, int *n, int *p, int *k, int *l, c *a, int *lda, c *b, int *ldb, s *alpha, s *beta, c *u, int *ldu, c *v, int *ldv, c *q, int *ldq, c *work, s *rwork, int *iwork, int *info) nogil
cdef cggsvd_t *cggsvd_f = &_fortran_cggsvd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cggsvp "F_FUNC(cggsvp,CGGSVP)"(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, c *a, int *lda, c *b, int *ldb, s *tola, s *tolb, int *k, int *l, c *u, int *ldu, c *v, int *ldv, c *q, int *ldq, int *iwork, s *rwork, c *tau, c *work, int *info) nogil
cdef cggsvp_t *cggsvp_f = &_fortran_cggsvp

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgtcon "F_FUNC(cgtcon,CGTCON)"(char *norm, int *n, c *dl, c *d, c *du, c *du2, int *ipiv, s *anorm, s *rcond, c *work, int *info) nogil
cdef cgtcon_t *cgtcon_f = &_fortran_cgtcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgtrfs "F_FUNC(cgtrfs,CGTRFS)"(char *trans, int *n, int *nrhs, c *dl, c *d, c *du, c *dlf, c *df, c *duf, c *du2, int *ipiv, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cgtrfs_t *cgtrfs_f = &_fortran_cgtrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgtsv "F_FUNC(cgtsv,CGTSV)"(int *n, int *nrhs, c *dl, c *d, c *du, c *b, int *ldb, int *info) nogil
cdef cgtsv_t *cgtsv_f = &_fortran_cgtsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgtsvx "F_FUNC(cgtsvx,CGTSVX)"(char *fact, char *trans, int *n, int *nrhs, c *dl, c *d, c *du, c *dlf, c *df, c *duf, c *du2, int *ipiv, c *b, int *ldb, c *x, int *ldx, s *rcond, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cgtsvx_t *cgtsvx_f = &_fortran_cgtsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgttrf "F_FUNC(cgttrf,CGTTRF)"(int *n, c *dl, c *d, c *du, c *du2, int *ipiv, int *info) nogil
cdef cgttrf_t *cgttrf_f = &_fortran_cgttrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgttrs "F_FUNC(cgttrs,CGTTRS)"(char *trans, int *n, int *nrhs, c *dl, c *d, c *du, c *du2, int *ipiv, c *b, int *ldb, int *info) nogil
cdef cgttrs_t *cgttrs_f = &_fortran_cgttrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cgtts2 "F_FUNC(cgtts2,CGTTS2)"(int *itrans, int *n, int *nrhs, c *dl, c *d, c *du, c *du2, int *ipiv, c *b, int *ldb) nogil
cdef cgtts2_t *cgtts2_f = &_fortran_cgtts2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chbev "F_FUNC(chbev,CHBEV)"(char *jobz, char *uplo, int *n, int *kd, c *ab, int *ldab, s *w, c *z, int *ldz, c *work, s *rwork, int *info) nogil
cdef chbev_t *chbev_f = &_fortran_chbev

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chbevd "F_FUNC(chbevd,CHBEVD)"(char *jobz, char *uplo, int *n, int *kd, c *ab, int *ldab, s *w, c *z, int *ldz, c *work, int *lwork, s *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef chbevd_t *chbevd_f = &_fortran_chbevd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chbevx "F_FUNC(chbevx,CHBEVX)"(char *jobz, char *range, char *uplo, int *n, int *kd, c *ab, int *ldab, c *q, int *ldq, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, c *z, int *ldz, c *work, s *rwork, int *iwork, int *ifail, int *info) nogil
cdef chbevx_t *chbevx_f = &_fortran_chbevx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chbgst "F_FUNC(chbgst,CHBGST)"(char *vect, char *uplo, int *n, int *ka, int *kb, c *ab, int *ldab, c *bb, int *ldbb, c *x, int *ldx, c *work, s *rwork, int *info) nogil
cdef chbgst_t *chbgst_f = &_fortran_chbgst

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chbgv "F_FUNC(chbgv,CHBGV)"(char *jobz, char *uplo, int *n, int *ka, int *kb, c *ab, int *ldab, c *bb, int *ldbb, s *w, c *z, int *ldz, c *work, s *rwork, int *info) nogil
cdef chbgv_t *chbgv_f = &_fortran_chbgv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chbgvd "F_FUNC(chbgvd,CHBGVD)"(char *jobz, char *uplo, int *n, int *ka, int *kb, c *ab, int *ldab, c *bb, int *ldbb, s *w, c *z, int *ldz, c *work, int *lwork, s *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef chbgvd_t *chbgvd_f = &_fortran_chbgvd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chbgvx "F_FUNC(chbgvx,CHBGVX)"(char *jobz, char *range, char *uplo, int *n, int *ka, int *kb, c *ab, int *ldab, c *bb, int *ldbb, c *q, int *ldq, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, c *z, int *ldz, c *work, s *rwork, int *iwork, int *ifail, int *info) nogil
cdef chbgvx_t *chbgvx_f = &_fortran_chbgvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chbtrd "F_FUNC(chbtrd,CHBTRD)"(char *vect, char *uplo, int *n, int *kd, c *ab, int *ldab, s *d, s *e, c *q, int *ldq, c *work, int *info) nogil
cdef chbtrd_t *chbtrd_f = &_fortran_chbtrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_checon "F_FUNC(checon,CHECON)"(char *uplo, int *n, c *a, int *lda, int *ipiv, s *anorm, s *rcond, c *work, int *info) nogil
cdef checon_t *checon_f = &_fortran_checon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cheev "F_FUNC(cheev,CHEEV)"(char *jobz, char *uplo, int *n, c *a, int *lda, s *w, c *work, int *lwork, s *rwork, int *info) nogil
cdef cheev_t *cheev_f = &_fortran_cheev

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cheevd "F_FUNC(cheevd,CHEEVD)"(char *jobz, char *uplo, int *n, c *a, int *lda, s *w, c *work, int *lwork, s *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef cheevd_t *cheevd_f = &_fortran_cheevd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cheevr "F_FUNC(cheevr,CHEEVR)"(char *jobz, char *range, char *uplo, int *n, c *a, int *lda, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, c *z, int *ldz, int *isuppz, c *work, int *lwork, s *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef cheevr_t *cheevr_f = &_fortran_cheevr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cheevx "F_FUNC(cheevx,CHEEVX)"(char *jobz, char *range, char *uplo, int *n, c *a, int *lda, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, c *z, int *ldz, c *work, int *lwork, s *rwork, int *iwork, int *ifail, int *info) nogil
cdef cheevx_t *cheevx_f = &_fortran_cheevx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chegs2 "F_FUNC(chegs2,CHEGS2)"(int *itype, char *uplo, int *n, c *a, int *lda, c *b, int *ldb, int *info) nogil
cdef chegs2_t *chegs2_f = &_fortran_chegs2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chegst "F_FUNC(chegst,CHEGST)"(int *itype, char *uplo, int *n, c *a, int *lda, c *b, int *ldb, int *info) nogil
cdef chegst_t *chegst_f = &_fortran_chegst

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chegv "F_FUNC(chegv,CHEGV)"(int *itype, char *jobz, char *uplo, int *n, c *a, int *lda, c *b, int *ldb, s *w, c *work, int *lwork, s *rwork, int *info) nogil
cdef chegv_t *chegv_f = &_fortran_chegv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chegvd "F_FUNC(chegvd,CHEGVD)"(int *itype, char *jobz, char *uplo, int *n, c *a, int *lda, c *b, int *ldb, s *w, c *work, int *lwork, s *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef chegvd_t *chegvd_f = &_fortran_chegvd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chegvx "F_FUNC(chegvx,CHEGVX)"(int *itype, char *jobz, char *range, char *uplo, int *n, c *a, int *lda, c *b, int *ldb, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, c *z, int *ldz, c *work, int *lwork, s *rwork, int *iwork, int *ifail, int *info) nogil
cdef chegvx_t *chegvx_f = &_fortran_chegvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cherfs "F_FUNC(cherfs,CHERFS)"(char *uplo, int *n, int *nrhs, c *a, int *lda, c *af, int *ldaf, int *ipiv, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cherfs_t *cherfs_f = &_fortran_cherfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chesv "F_FUNC(chesv,CHESV)"(char *uplo, int *n, int *nrhs, c *a, int *lda, int *ipiv, c *b, int *ldb, c *work, int *lwork, int *info) nogil
cdef chesv_t *chesv_f = &_fortran_chesv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chesvx "F_FUNC(chesvx,CHESVX)"(char *fact, char *uplo, int *n, int *nrhs, c *a, int *lda, c *af, int *ldaf, int *ipiv, c *b, int *ldb, c *x, int *ldx, s *rcond, s *ferr, s *berr, c *work, int *lwork, s *rwork, int *info) nogil
cdef chesvx_t *chesvx_f = &_fortran_chesvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chetd2 "F_FUNC(chetd2,CHETD2)"(char *uplo, int *n, c *a, int *lda, s *d, s *e, c *tau, int *info) nogil
cdef chetd2_t *chetd2_f = &_fortran_chetd2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chetf2 "F_FUNC(chetf2,CHETF2)"(char *uplo, int *n, c *a, int *lda, int *ipiv, int *info) nogil
cdef chetf2_t *chetf2_f = &_fortran_chetf2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chetrd "F_FUNC(chetrd,CHETRD)"(char *uplo, int *n, c *a, int *lda, s *d, s *e, c *tau, c *work, int *lwork, int *info) nogil
cdef chetrd_t *chetrd_f = &_fortran_chetrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chetrf "F_FUNC(chetrf,CHETRF)"(char *uplo, int *n, c *a, int *lda, int *ipiv, c *work, int *lwork, int *info) nogil
cdef chetrf_t *chetrf_f = &_fortran_chetrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chetri "F_FUNC(chetri,CHETRI)"(char *uplo, int *n, c *a, int *lda, int *ipiv, c *work, int *info) nogil
cdef chetri_t *chetri_f = &_fortran_chetri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chetrs "F_FUNC(chetrs,CHETRS)"(char *uplo, int *n, int *nrhs, c *a, int *lda, int *ipiv, c *b, int *ldb, int *info) nogil
cdef chetrs_t *chetrs_f = &_fortran_chetrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chgeqz "F_FUNC(chgeqz,CHGEQZ)"(char *job, char *compq, char *compz, int *n, int *ilo, int *ihi, c *h, int *ldh, c *t, int *ldt, c *alpha, c *beta, c *q, int *ldq, c *z, int *ldz, c *work, int *lwork, s *rwork, int *info) nogil
cdef chgeqz_t *chgeqz_f = &_fortran_chgeqz

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chpcon "F_FUNC(chpcon,CHPCON)"(char *uplo, int *n, c *ap, int *ipiv, s *anorm, s *rcond, c *work, int *info) nogil
cdef chpcon_t *chpcon_f = &_fortran_chpcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chpev "F_FUNC(chpev,CHPEV)"(char *jobz, char *uplo, int *n, c *ap, s *w, c *z, int *ldz, c *work, s *rwork, int *info) nogil
cdef chpev_t *chpev_f = &_fortran_chpev

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chpevd "F_FUNC(chpevd,CHPEVD)"(char *jobz, char *uplo, int *n, c *ap, s *w, c *z, int *ldz, c *work, int *lwork, s *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef chpevd_t *chpevd_f = &_fortran_chpevd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chpevx "F_FUNC(chpevx,CHPEVX)"(char *jobz, char *range, char *uplo, int *n, c *ap, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, c *z, int *ldz, c *work, s *rwork, int *iwork, int *ifail, int *info) nogil
cdef chpevx_t *chpevx_f = &_fortran_chpevx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chpgst "F_FUNC(chpgst,CHPGST)"(int *itype, char *uplo, int *n, c *ap, c *bp, int *info) nogil
cdef chpgst_t *chpgst_f = &_fortran_chpgst

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chpgv "F_FUNC(chpgv,CHPGV)"(int *itype, char *jobz, char *uplo, int *n, c *ap, c *bp, s *w, c *z, int *ldz, c *work, s *rwork, int *info) nogil
cdef chpgv_t *chpgv_f = &_fortran_chpgv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chpgvd "F_FUNC(chpgvd,CHPGVD)"(int *itype, char *jobz, char *uplo, int *n, c *ap, c *bp, s *w, c *z, int *ldz, c *work, int *lwork, s *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef chpgvd_t *chpgvd_f = &_fortran_chpgvd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chpgvx "F_FUNC(chpgvx,CHPGVX)"(int *itype, char *jobz, char *range, char *uplo, int *n, c *ap, c *bp, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, c *z, int *ldz, c *work, s *rwork, int *iwork, int *ifail, int *info) nogil
cdef chpgvx_t *chpgvx_f = &_fortran_chpgvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chprfs "F_FUNC(chprfs,CHPRFS)"(char *uplo, int *n, int *nrhs, c *ap, c *afp, int *ipiv, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef chprfs_t *chprfs_f = &_fortran_chprfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chpsv "F_FUNC(chpsv,CHPSV)"(char *uplo, int *n, int *nrhs, c *ap, int *ipiv, c *b, int *ldb, int *info) nogil
cdef chpsv_t *chpsv_f = &_fortran_chpsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chpsvx "F_FUNC(chpsvx,CHPSVX)"(char *fact, char *uplo, int *n, int *nrhs, c *ap, c *afp, int *ipiv, c *b, int *ldb, c *x, int *ldx, s *rcond, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef chpsvx_t *chpsvx_f = &_fortran_chpsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chptrd "F_FUNC(chptrd,CHPTRD)"(char *uplo, int *n, c *ap, s *d, s *e, c *tau, int *info) nogil
cdef chptrd_t *chptrd_f = &_fortran_chptrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chptrf "F_FUNC(chptrf,CHPTRF)"(char *uplo, int *n, c *ap, int *ipiv, int *info) nogil
cdef chptrf_t *chptrf_f = &_fortran_chptrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chptri "F_FUNC(chptri,CHPTRI)"(char *uplo, int *n, c *ap, int *ipiv, c *work, int *info) nogil
cdef chptri_t *chptri_f = &_fortran_chptri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chptrs "F_FUNC(chptrs,CHPTRS)"(char *uplo, int *n, int *nrhs, c *ap, int *ipiv, c *b, int *ldb, int *info) nogil
cdef chptrs_t *chptrs_f = &_fortran_chptrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chsein "F_FUNC(chsein,CHSEIN)"(char *side, char *eigsrc, char *initv, bint *select, int *n, c *h, int *ldh, c *w, c *vl, int *ldvl, c *vr, int *ldvr, int *mm, int *m, c *work, s *rwork, int *ifaill, int *ifailr, int *info) nogil
cdef chsein_t *chsein_f = &_fortran_chsein

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_chseqr "F_FUNC(chseqr,CHSEQR)"(char *job, char *compz, int *n, int *ilo, int *ihi, c *h, int *ldh, c *w, c *z, int *ldz, c *work, int *lwork, int *info) nogil
cdef chseqr_t *chseqr_f = &_fortran_chseqr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_clacn2 "F_FUNC(clacn2,CLACN2)"(int *n, c *v, c *x, s *est, int *kase, int *isave) nogil
cdef clacn2_t *clacn2_f = &_fortran_clacn2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_clacon "F_FUNC(clacon,CLACON)"(int *n, c *v, c *x, s *est, int *kase) nogil
cdef clacon_t *clacon_f = &_fortran_clacon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_clarf "F_FUNC(clarf,CLARF)"(char *side, int *m, int *n, c *v, int *incv, c *tau, c *c, int *ldc, c *work) nogil
cdef clarf_t *clarf_f = &_fortran_clarf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_clarz "F_FUNC(clarz,CLARZ)"(char *side, int *m, int *n, int *l, c *v, int *incv, c *tau, c *c, int *ldc, c *work) nogil
cdef clarz_t *clarz_f = &_fortran_clarz

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_claswp "F_FUNC(claswp,CLASWP)"(int *n, c *a, int *lda, int *k1, int *k2, int *ipiv, int *incx) nogil
cdef claswp_t *claswp_f = &_fortran_claswp

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_clauum "F_FUNC(clauum,CLAUUM)"(char *uplo, int *n, c *a, int *lda, int *info) nogil
cdef clauum_t *clauum_f = &_fortran_clauum

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cpbcon "F_FUNC(cpbcon,CPBCON)"(char *uplo, int *n, int *kd, c *ab, int *ldab, s *anorm, s *rcond, c *work, s *rwork, int *info) nogil
cdef cpbcon_t *cpbcon_f = &_fortran_cpbcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cpbequ "F_FUNC(cpbequ,CPBEQU)"(char *uplo, int *n, int *kd, c *ab, int *ldab, s *s, s *scond, s *amax, int *info) nogil
cdef cpbequ_t *cpbequ_f = &_fortran_cpbequ

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cpbrfs "F_FUNC(cpbrfs,CPBRFS)"(char *uplo, int *n, int *kd, int *nrhs, c *ab, int *ldab, c *afb, int *ldafb, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cpbrfs_t *cpbrfs_f = &_fortran_cpbrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cpbstf "F_FUNC(cpbstf,CPBSTF)"(char *uplo, int *n, int *kd, c *ab, int *ldab, int *info) nogil
cdef cpbstf_t *cpbstf_f = &_fortran_cpbstf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cpbsv "F_FUNC(cpbsv,CPBSV)"(char *uplo, int *n, int *kd, int *nrhs, c *ab, int *ldab, c *b, int *ldb, int *info) nogil
cdef cpbsv_t *cpbsv_f = &_fortran_cpbsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cpbsvx "F_FUNC(cpbsvx,CPBSVX)"(char *fact, char *uplo, int *n, int *kd, int *nrhs, c *ab, int *ldab, c *afb, int *ldafb, char *equed, s *s, c *b, int *ldb, c *x, int *ldx, s *rcond, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cpbsvx_t *cpbsvx_f = &_fortran_cpbsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cpbtf2 "F_FUNC(cpbtf2,CPBTF2)"(char *uplo, int *n, int *kd, c *ab, int *ldab, int *info) nogil
cdef cpbtf2_t *cpbtf2_f = &_fortran_cpbtf2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cpbtrf "F_FUNC(cpbtrf,CPBTRF)"(char *uplo, int *n, int *kd, c *ab, int *ldab, int *info) nogil
cdef cpbtrf_t *cpbtrf_f = &_fortran_cpbtrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cpbtrs "F_FUNC(cpbtrs,CPBTRS)"(char *uplo, int *n, int *kd, int *nrhs, c *ab, int *ldab, c *b, int *ldb, int *info) nogil
cdef cpbtrs_t *cpbtrs_f = &_fortran_cpbtrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cpocon "F_FUNC(cpocon,CPOCON)"(char *uplo, int *n, c *a, int *lda, s *anorm, s *rcond, c *work, s *rwork, int *info) nogil
cdef cpocon_t *cpocon_f = &_fortran_cpocon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cpoequ "F_FUNC(cpoequ,CPOEQU)"(int *n, c *a, int *lda, s *s, s *scond, s *amax, int *info) nogil
cdef cpoequ_t *cpoequ_f = &_fortran_cpoequ

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cporfs "F_FUNC(cporfs,CPORFS)"(char *uplo, int *n, int *nrhs, c *a, int *lda, c *af, int *ldaf, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cporfs_t *cporfs_f = &_fortran_cporfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cposv "F_FUNC(cposv,CPOSV)"(char *uplo, int *n, int *nrhs, c *a, int *lda, c *b, int *ldb, int *info) nogil
cdef cposv_t *cposv_f = &_fortran_cposv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cposvx "F_FUNC(cposvx,CPOSVX)"(char *fact, char *uplo, int *n, int *nrhs, c *a, int *lda, c *af, int *ldaf, char *equed, s *s, c *b, int *ldb, c *x, int *ldx, s *rcond, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cposvx_t *cposvx_f = &_fortran_cposvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cpotf2 "F_FUNC(cpotf2,CPOTF2)"(char *uplo, int *n, c *a, int *lda, int *info) nogil
cdef cpotf2_t *cpotf2_f = &_fortran_cpotf2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cpotrf "F_FUNC(cpotrf,CPOTRF)"(char *uplo, int *n, c *a, int *lda, int *info) nogil
cdef cpotrf_t *cpotrf_f = &_fortran_cpotrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cpotri "F_FUNC(cpotri,CPOTRI)"(char *uplo, int *n, c *a, int *lda, int *info) nogil
cdef cpotri_t *cpotri_f = &_fortran_cpotri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cpotrs "F_FUNC(cpotrs,CPOTRS)"(char *uplo, int *n, int *nrhs, c *a, int *lda, c *b, int *ldb, int *info) nogil
cdef cpotrs_t *cpotrs_f = &_fortran_cpotrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cppcon "F_FUNC(cppcon,CPPCON)"(char *uplo, int *n, c *ap, s *anorm, s *rcond, c *work, s *rwork, int *info) nogil
cdef cppcon_t *cppcon_f = &_fortran_cppcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cppequ "F_FUNC(cppequ,CPPEQU)"(char *uplo, int *n, c *ap, s *s, s *scond, s *amax, int *info) nogil
cdef cppequ_t *cppequ_f = &_fortran_cppequ

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cpprfs "F_FUNC(cpprfs,CPPRFS)"(char *uplo, int *n, int *nrhs, c *ap, c *afp, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cpprfs_t *cpprfs_f = &_fortran_cpprfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cppsv "F_FUNC(cppsv,CPPSV)"(char *uplo, int *n, int *nrhs, c *ap, c *b, int *ldb, int *info) nogil
cdef cppsv_t *cppsv_f = &_fortran_cppsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cppsvx "F_FUNC(cppsvx,CPPSVX)"(char *fact, char *uplo, int *n, int *nrhs, c *ap, c *afp, char *equed, s *s, c *b, int *ldb, c *x, int *ldx, s *rcond, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cppsvx_t *cppsvx_f = &_fortran_cppsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cpptrf "F_FUNC(cpptrf,CPPTRF)"(char *uplo, int *n, c *ap, int *info) nogil
cdef cpptrf_t *cpptrf_f = &_fortran_cpptrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cpptri "F_FUNC(cpptri,CPPTRI)"(char *uplo, int *n, c *ap, int *info) nogil
cdef cpptri_t *cpptri_f = &_fortran_cpptri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cpptrs "F_FUNC(cpptrs,CPPTRS)"(char *uplo, int *n, int *nrhs, c *ap, c *b, int *ldb, int *info) nogil
cdef cpptrs_t *cpptrs_f = &_fortran_cpptrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cptcon "F_FUNC(cptcon,CPTCON)"(int *n, s *d, c *e, s *anorm, s *rcond, s *rwork, int *info) nogil
cdef cptcon_t *cptcon_f = &_fortran_cptcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cpteqr "F_FUNC(cpteqr,CPTEQR)"(char *compz, int *n, s *d, s *e, c *z, int *ldz, s *work, int *info) nogil
cdef cpteqr_t *cpteqr_f = &_fortran_cpteqr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cptrfs "F_FUNC(cptrfs,CPTRFS)"(char *uplo, int *n, int *nrhs, s *d, c *e, s *df, c *ef, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cptrfs_t *cptrfs_f = &_fortran_cptrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cptsv "F_FUNC(cptsv,CPTSV)"(int *n, int *nrhs, s *d, c *e, c *b, int *ldb, int *info) nogil
cdef cptsv_t *cptsv_f = &_fortran_cptsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cptsvx "F_FUNC(cptsvx,CPTSVX)"(char *fact, int *n, int *nrhs, s *d, c *e, s *df, c *ef, c *b, int *ldb, c *x, int *ldx, s *rcond, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cptsvx_t *cptsvx_f = &_fortran_cptsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cpttrf "F_FUNC(cpttrf,CPTTRF)"(int *n, s *d, c *e, int *info) nogil
cdef cpttrf_t *cpttrf_f = &_fortran_cpttrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cpttrs "F_FUNC(cpttrs,CPTTRS)"(char *uplo, int *n, int *nrhs, s *d, c *e, c *b, int *ldb, int *info) nogil
cdef cpttrs_t *cpttrs_f = &_fortran_cpttrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cptts2 "F_FUNC(cptts2,CPTTS2)"(int *iuplo, int *n, int *nrhs, s *d, c *e, c *b, int *ldb) nogil
cdef cptts2_t *cptts2_f = &_fortran_cptts2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_crot "F_FUNC(crot,CROT)"(int *n, c *cx, int *incx, c *cy, int *incy, s *c, c *s) nogil
cdef crot_t *crot_f = &_fortran_crot

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cspcon "F_FUNC(cspcon,CSPCON)"(char *uplo, int *n, c *ap, int *ipiv, s *anorm, s *rcond, c *work, int *info) nogil
cdef cspcon_t *cspcon_f = &_fortran_cspcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cspmv "F_FUNC(cspmv,CSPMV)"(char *uplo, int *n, c *alpha, c *ap, c *x, int *incx, c *beta, c *y, int *incy) nogil
cdef cspmv_t *cspmv_f = &_fortran_cspmv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cspr "F_FUNC(cspr,CSPR)"(char *uplo, int *n, c *alpha, c *x, int *incx, c *ap) nogil
cdef cspr_t *cspr_f = &_fortran_cspr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_csprfs "F_FUNC(csprfs,CSPRFS)"(char *uplo, int *n, int *nrhs, c *ap, c *afp, int *ipiv, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef csprfs_t *csprfs_f = &_fortran_csprfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cspsv "F_FUNC(cspsv,CSPSV)"(char *uplo, int *n, int *nrhs, c *ap, int *ipiv, c *b, int *ldb, int *info) nogil
cdef cspsv_t *cspsv_f = &_fortran_cspsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cspsvx "F_FUNC(cspsvx,CSPSVX)"(char *fact, char *uplo, int *n, int *nrhs, c *ap, c *afp, int *ipiv, c *b, int *ldb, c *x, int *ldx, s *rcond, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef cspsvx_t *cspsvx_f = &_fortran_cspsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_csptrf "F_FUNC(csptrf,CSPTRF)"(char *uplo, int *n, c *ap, int *ipiv, int *info) nogil
cdef csptrf_t *csptrf_f = &_fortran_csptrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_csptri "F_FUNC(csptri,CSPTRI)"(char *uplo, int *n, c *ap, int *ipiv, c *work, int *info) nogil
cdef csptri_t *csptri_f = &_fortran_csptri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_csptrs "F_FUNC(csptrs,CSPTRS)"(char *uplo, int *n, int *nrhs, c *ap, int *ipiv, c *b, int *ldb, int *info) nogil
cdef csptrs_t *csptrs_f = &_fortran_csptrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_csrscl "F_FUNC(csrscl,CSRSCL)"(int *n, s *sa, c *sx, int *incx) nogil
cdef csrscl_t *csrscl_f = &_fortran_csrscl

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cstedc "F_FUNC(cstedc,CSTEDC)"(char *compz, int *n, s *d, s *e, c *z, int *ldz, c *work, int *lwork, s *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef cstedc_t *cstedc_f = &_fortran_cstedc

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cstegr "F_FUNC(cstegr,CSTEGR)"(char *jobz, char *range, int *n, s *d, s *e, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, c *z, int *ldz, int *isuppz, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef cstegr_t *cstegr_f = &_fortran_cstegr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cstein "F_FUNC(cstein,CSTEIN)"(int *n, s *d, s *e, int *m, s *w, int *iblock, int *isplit, c *z, int *ldz, s *work, int *iwork, int *ifail, int *info) nogil
cdef cstein_t *cstein_f = &_fortran_cstein

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cstemr "F_FUNC(cstemr,CSTEMR)"(char *jobz, char *range, int *n, s *d, s *e, s *vl, s *vu, int *il, int *iu, int *m, s *w, c *z, int *ldz, int *nzc, int *isuppz, bint *tryrac, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef cstemr_t *cstemr_f = &_fortran_cstemr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_csteqr "F_FUNC(csteqr,CSTEQR)"(char *compz, int *n, s *d, s *e, c *z, int *ldz, s *work, int *info) nogil
cdef csteqr_t *csteqr_f = &_fortran_csteqr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_csycon "F_FUNC(csycon,CSYCON)"(char *uplo, int *n, c *a, int *lda, int *ipiv, s *anorm, s *rcond, c *work, int *info) nogil
cdef csycon_t *csycon_f = &_fortran_csycon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_csymv "F_FUNC(csymv,CSYMV)"(char *uplo, int *n, c *alpha, c *a, int *lda, c *x, int *incx, c *beta, c *y, int *incy) nogil
cdef csymv_t *csymv_f = &_fortran_csymv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_csyr "F_FUNC(csyr,CSYR)"(char *uplo, int *n, c *alpha, c *x, int *incx, c *a, int *lda) nogil
cdef csyr_t *csyr_f = &_fortran_csyr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_csyrfs "F_FUNC(csyrfs,CSYRFS)"(char *uplo, int *n, int *nrhs, c *a, int *lda, c *af, int *ldaf, int *ipiv, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef csyrfs_t *csyrfs_f = &_fortran_csyrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_csysv "F_FUNC(csysv,CSYSV)"(char *uplo, int *n, int *nrhs, c *a, int *lda, int *ipiv, c *b, int *ldb, c *work, int *lwork, int *info) nogil
cdef csysv_t *csysv_f = &_fortran_csysv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_csysvx "F_FUNC(csysvx,CSYSVX)"(char *fact, char *uplo, int *n, int *nrhs, c *a, int *lda, c *af, int *ldaf, int *ipiv, c *b, int *ldb, c *x, int *ldx, s *rcond, s *ferr, s *berr, c *work, int *lwork, s *rwork, int *info) nogil
cdef csysvx_t *csysvx_f = &_fortran_csysvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_csytf2 "F_FUNC(csytf2,CSYTF2)"(char *uplo, int *n, c *a, int *lda, int *ipiv, int *info) nogil
cdef csytf2_t *csytf2_f = &_fortran_csytf2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_csytrf "F_FUNC(csytrf,CSYTRF)"(char *uplo, int *n, c *a, int *lda, int *ipiv, c *work, int *lwork, int *info) nogil
cdef csytrf_t *csytrf_f = &_fortran_csytrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_csytri "F_FUNC(csytri,CSYTRI)"(char *uplo, int *n, c *a, int *lda, int *ipiv, c *work, int *info) nogil
cdef csytri_t *csytri_f = &_fortran_csytri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_csytrs "F_FUNC(csytrs,CSYTRS)"(char *uplo, int *n, int *nrhs, c *a, int *lda, int *ipiv, c *b, int *ldb, int *info) nogil
cdef csytrs_t *csytrs_f = &_fortran_csytrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctbcon "F_FUNC(ctbcon,CTBCON)"(char *norm, char *uplo, char *diag, int *n, int *kd, c *ab, int *ldab, s *rcond, c *work, s *rwork, int *info) nogil
cdef ctbcon_t *ctbcon_f = &_fortran_ctbcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctbrfs "F_FUNC(ctbrfs,CTBRFS)"(char *uplo, char *trans, char *diag, int *n, int *kd, int *nrhs, c *ab, int *ldab, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef ctbrfs_t *ctbrfs_f = &_fortran_ctbrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctbtrs "F_FUNC(ctbtrs,CTBTRS)"(char *uplo, char *trans, char *diag, int *n, int *kd, int *nrhs, c *ab, int *ldab, c *b, int *ldb, int *info) nogil
cdef ctbtrs_t *ctbtrs_f = &_fortran_ctbtrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctgevc "F_FUNC(ctgevc,CTGEVC)"(char *side, char *howmny, bint *select, int *n, c *s, int *lds, c *p, int *ldp, c *vl, int *ldvl, c *vr, int *ldvr, int *mm, int *m, c *work, s *rwork, int *info) nogil
cdef ctgevc_t *ctgevc_f = &_fortran_ctgevc

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctgex2 "F_FUNC(ctgex2,CTGEX2)"(bint *wantq, bint *wantz, int *n, c *a, int *lda, c *b, int *ldb, c *q, int *ldq, c *z, int *ldz, int *j1, int *info) nogil
cdef ctgex2_t *ctgex2_f = &_fortran_ctgex2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctgexc "F_FUNC(ctgexc,CTGEXC)"(bint *wantq, bint *wantz, int *n, c *a, int *lda, c *b, int *ldb, c *q, int *ldq, c *z, int *ldz, int *ifst, int *ilst, int *info) nogil
cdef ctgexc_t *ctgexc_f = &_fortran_ctgexc

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctgsen "F_FUNC(ctgsen,CTGSEN)"(int *ijob, bint *wantq, bint *wantz, bint *select, int *n, c *a, int *lda, c *b, int *ldb, c *alpha, c *beta, c *q, int *ldq, c *z, int *ldz, int *m, s *pl, s *pr, s *dif, c *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef ctgsen_t *ctgsen_f = &_fortran_ctgsen

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctgsja "F_FUNC(ctgsja,CTGSJA)"(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, int *k, int *l, c *a, int *lda, c *b, int *ldb, s *tola, s *tolb, s *alpha, s *beta, c *u, int *ldu, c *v, int *ldv, c *q, int *ldq, c *work, int *ncycle, int *info) nogil
cdef ctgsja_t *ctgsja_f = &_fortran_ctgsja

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctgsna "F_FUNC(ctgsna,CTGSNA)"(char *job, char *howmny, bint *select, int *n, c *a, int *lda, c *b, int *ldb, c *vl, int *ldvl, c *vr, int *ldvr, s *s, s *dif, int *mm, int *m, c *work, int *lwork, int *iwork, int *info) nogil
cdef ctgsna_t *ctgsna_f = &_fortran_ctgsna

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctgsy2 "F_FUNC(ctgsy2,CTGSY2)"(char *trans, int *ijob, int *m, int *n, c *a, int *lda, c *b, int *ldb, c *c, int *ldc, c *d, int *ldd, c *e, int *lde, c *f, int *ldf, s *scale, s *rdsum, s *rdscal, int *info) nogil
cdef ctgsy2_t *ctgsy2_f = &_fortran_ctgsy2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctgsyl "F_FUNC(ctgsyl,CTGSYL)"(char *trans, int *ijob, int *m, int *n, c *a, int *lda, c *b, int *ldb, c *c, int *ldc, c *d, int *ldd, c *e, int *lde, c *f, int *ldf, s *scale, s *dif, c *work, int *lwork, int *iwork, int *info) nogil
cdef ctgsyl_t *ctgsyl_f = &_fortran_ctgsyl

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctpcon "F_FUNC(ctpcon,CTPCON)"(char *norm, char *uplo, char *diag, int *n, c *ap, s *rcond, c *work, s *rwork, int *info) nogil
cdef ctpcon_t *ctpcon_f = &_fortran_ctpcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctprfs "F_FUNC(ctprfs,CTPRFS)"(char *uplo, char *trans, char *diag, int *n, int *nrhs, c *ap, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef ctprfs_t *ctprfs_f = &_fortran_ctprfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctptri "F_FUNC(ctptri,CTPTRI)"(char *uplo, char *diag, int *n, c *ap, int *info) nogil
cdef ctptri_t *ctptri_f = &_fortran_ctptri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctptrs "F_FUNC(ctptrs,CTPTRS)"(char *uplo, char *trans, char *diag, int *n, int *nrhs, c *ap, c *b, int *ldb, int *info) nogil
cdef ctptrs_t *ctptrs_f = &_fortran_ctptrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctrcon "F_FUNC(ctrcon,CTRCON)"(char *norm, char *uplo, char *diag, int *n, c *a, int *lda, s *rcond, c *work, s *rwork, int *info) nogil
cdef ctrcon_t *ctrcon_f = &_fortran_ctrcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctrevc "F_FUNC(ctrevc,CTREVC)"(char *side, char *howmny, bint *select, int *n, c *t, int *ldt, c *vl, int *ldvl, c *vr, int *ldvr, int *mm, int *m, c *work, s *rwork, int *info) nogil
cdef ctrevc_t *ctrevc_f = &_fortran_ctrevc

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctrexc "F_FUNC(ctrexc,CTREXC)"(char *compq, int *n, c *t, int *ldt, c *q, int *ldq, int *ifst, int *ilst, int *info) nogil
cdef ctrexc_t *ctrexc_f = &_fortran_ctrexc

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctrrfs "F_FUNC(ctrrfs,CTRRFS)"(char *uplo, char *trans, char *diag, int *n, int *nrhs, c *a, int *lda, c *b, int *ldb, c *x, int *ldx, s *ferr, s *berr, c *work, s *rwork, int *info) nogil
cdef ctrrfs_t *ctrrfs_f = &_fortran_ctrrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctrsen "F_FUNC(ctrsen,CTRSEN)"(char *job, char *compq, bint *select, int *n, c *t, int *ldt, c *q, int *ldq, c *w, int *m, s *s, s *sep, c *work, int *lwork, int *info) nogil
cdef ctrsen_t *ctrsen_f = &_fortran_ctrsen

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctrsna "F_FUNC(ctrsna,CTRSNA)"(char *job, char *howmny, bint *select, int *n, c *t, int *ldt, c *vl, int *ldvl, c *vr, int *ldvr, s *s, s *sep, int *mm, int *m, c *work, int *ldwork, s *rwork, int *info) nogil
cdef ctrsna_t *ctrsna_f = &_fortran_ctrsna

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctrsyl "F_FUNC(ctrsyl,CTRSYL)"(char *trana, char *tranb, int *isgn, int *m, int *n, c *a, int *lda, c *b, int *ldb, c *c, int *ldc, s *scale, int *info) nogil
cdef ctrsyl_t *ctrsyl_f = &_fortran_ctrsyl

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctrti2 "F_FUNC(ctrti2,CTRTI2)"(char *uplo, char *diag, int *n, c *a, int *lda, int *info) nogil
cdef ctrti2_t *ctrti2_f = &_fortran_ctrti2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctrtri "F_FUNC(ctrtri,CTRTRI)"(char *uplo, char *diag, int *n, c *a, int *lda, int *info) nogil
cdef ctrtri_t *ctrtri_f = &_fortran_ctrtri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctrtrs "F_FUNC(ctrtrs,CTRTRS)"(char *uplo, char *trans, char *diag, int *n, int *nrhs, c *a, int *lda, c *b, int *ldb, int *info) nogil
cdef ctrtrs_t *ctrtrs_f = &_fortran_ctrtrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctzrqf "F_FUNC(ctzrqf,CTZRQF)"(int *m, int *n, c *a, int *lda, c *tau, int *info) nogil
cdef ctzrqf_t *ctzrqf_f = &_fortran_ctzrqf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ctzrzf "F_FUNC(ctzrzf,CTZRZF)"(int *m, int *n, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef ctzrzf_t *ctzrzf_f = &_fortran_ctzrzf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cung2l "F_FUNC(cung2l,CUNG2L)"(int *m, int *n, int *k, c *a, int *lda, c *tau, c *work, int *info) nogil
cdef cung2l_t *cung2l_f = &_fortran_cung2l

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cung2r "F_FUNC(cung2r,CUNG2R)"(int *m, int *n, int *k, c *a, int *lda, c *tau, c *work, int *info) nogil
cdef cung2r_t *cung2r_f = &_fortran_cung2r

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cungbr "F_FUNC(cungbr,CUNGBR)"(char *vect, int *m, int *n, int *k, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef cungbr_t *cungbr_f = &_fortran_cungbr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cunghr "F_FUNC(cunghr,CUNGHR)"(int *n, int *ilo, int *ihi, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef cunghr_t *cunghr_f = &_fortran_cunghr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cungl2 "F_FUNC(cungl2,CUNGL2)"(int *m, int *n, int *k, c *a, int *lda, c *tau, c *work, int *info) nogil
cdef cungl2_t *cungl2_f = &_fortran_cungl2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cunglq "F_FUNC(cunglq,CUNGLQ)"(int *m, int *n, int *k, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef cunglq_t *cunglq_f = &_fortran_cunglq

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cungql "F_FUNC(cungql,CUNGQL)"(int *m, int *n, int *k, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef cungql_t *cungql_f = &_fortran_cungql

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cungqr "F_FUNC(cungqr,CUNGQR)"(int *m, int *n, int *k, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef cungqr_t *cungqr_f = &_fortran_cungqr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cungr2 "F_FUNC(cungr2,CUNGR2)"(int *m, int *n, int *k, c *a, int *lda, c *tau, c *work, int *info) nogil
cdef cungr2_t *cungr2_f = &_fortran_cungr2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cungrq "F_FUNC(cungrq,CUNGRQ)"(int *m, int *n, int *k, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef cungrq_t *cungrq_f = &_fortran_cungrq

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cungtr "F_FUNC(cungtr,CUNGTR)"(char *uplo, int *n, c *a, int *lda, c *tau, c *work, int *lwork, int *info) nogil
cdef cungtr_t *cungtr_f = &_fortran_cungtr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cunm2l "F_FUNC(cunm2l,CUNM2L)"(char *side, char *trans, int *m, int *n, int *k, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *info) nogil
cdef cunm2l_t *cunm2l_f = &_fortran_cunm2l

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cunm2r "F_FUNC(cunm2r,CUNM2R)"(char *side, char *trans, int *m, int *n, int *k, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *info) nogil
cdef cunm2r_t *cunm2r_f = &_fortran_cunm2r

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cunmbr "F_FUNC(cunmbr,CUNMBR)"(char *vect, char *side, char *trans, int *m, int *n, int *k, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *lwork, int *info) nogil
cdef cunmbr_t *cunmbr_f = &_fortran_cunmbr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cunmhr "F_FUNC(cunmhr,CUNMHR)"(char *side, char *trans, int *m, int *n, int *ilo, int *ihi, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *lwork, int *info) nogil
cdef cunmhr_t *cunmhr_f = &_fortran_cunmhr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cunml2 "F_FUNC(cunml2,CUNML2)"(char *side, char *trans, int *m, int *n, int *k, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *info) nogil
cdef cunml2_t *cunml2_f = &_fortran_cunml2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cunmlq "F_FUNC(cunmlq,CUNMLQ)"(char *side, char *trans, int *m, int *n, int *k, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *lwork, int *info) nogil
cdef cunmlq_t *cunmlq_f = &_fortran_cunmlq

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cunmql "F_FUNC(cunmql,CUNMQL)"(char *side, char *trans, int *m, int *n, int *k, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *lwork, int *info) nogil
cdef cunmql_t *cunmql_f = &_fortran_cunmql

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cunmqr "F_FUNC(cunmqr,CUNMQR)"(char *side, char *trans, int *m, int *n, int *k, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *lwork, int *info) nogil
cdef cunmqr_t *cunmqr_f = &_fortran_cunmqr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cunmr2 "F_FUNC(cunmr2,CUNMR2)"(char *side, char *trans, int *m, int *n, int *k, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *info) nogil
cdef cunmr2_t *cunmr2_f = &_fortran_cunmr2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cunmr3 "F_FUNC(cunmr3,CUNMR3)"(char *side, char *trans, int *m, int *n, int *k, int *l, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *info) nogil
cdef cunmr3_t *cunmr3_f = &_fortran_cunmr3

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cunmrq "F_FUNC(cunmrq,CUNMRQ)"(char *side, char *trans, int *m, int *n, int *k, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *lwork, int *info) nogil
cdef cunmrq_t *cunmrq_f = &_fortran_cunmrq

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cunmrz "F_FUNC(cunmrz,CUNMRZ)"(char *side, char *trans, int *m, int *n, int *k, int *l, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *lwork, int *info) nogil
cdef cunmrz_t *cunmrz_f = &_fortran_cunmrz

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cunmtr "F_FUNC(cunmtr,CUNMTR)"(char *side, char *uplo, char *trans, int *m, int *n, c *a, int *lda, c *tau, c *c, int *ldc, c *work, int *lwork, int *info) nogil
cdef cunmtr_t *cunmtr_f = &_fortran_cunmtr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cupgtr "F_FUNC(cupgtr,CUPGTR)"(char *uplo, int *n, c *ap, c *tau, c *q, int *ldq, c *work, int *info) nogil
cdef cupgtr_t *cupgtr_f = &_fortran_cupgtr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_cupmtr "F_FUNC(cupmtr,CUPMTR)"(char *side, char *uplo, char *trans, int *m, int *n, c *ap, c *tau, c *c, int *ldc, c *work, int *info) nogil
cdef cupmtr_t *cupmtr_f = &_fortran_cupmtr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dbdsdc "F_FUNC(dbdsdc,DBDSDC)"(char *uplo, char *compq, int *n, d *d, d *e, d *u, int *ldu, d *vt, int *ldvt, d *q, int *iq, d *work, int *iwork, int *info) nogil
cdef dbdsdc_t *dbdsdc_f = &_fortran_dbdsdc

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dbdsqr "F_FUNC(dbdsqr,DBDSQR)"(char *uplo, int *n, int *ncvt, int *nru, int *ncc, d *d, d *e, d *vt, int *ldvt, d *u, int *ldu, d *c, int *ldc, d *work, int *info) nogil
cdef dbdsqr_t *dbdsqr_f = &_fortran_dbdsqr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ddisna "F_FUNC(ddisna,DDISNA)"(char *job, int *m, int *n, d *d, d *sep, int *info) nogil
cdef ddisna_t *ddisna_f = &_fortran_ddisna

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgbbrd "F_FUNC(dgbbrd,DGBBRD)"(char *vect, int *m, int *n, int *ncc, int *kl, int *ku, d *ab, int *ldab, d *d, d *e, d *q, int *ldq, d *pt, int *ldpt, d *c, int *ldc, d *work, int *info) nogil
cdef dgbbrd_t *dgbbrd_f = &_fortran_dgbbrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgbcon "F_FUNC(dgbcon,DGBCON)"(char *norm, int *n, int *kl, int *ku, d *ab, int *ldab, int *ipiv, d *anorm, d *rcond, d *work, int *iwork, int *info) nogil
cdef dgbcon_t *dgbcon_f = &_fortran_dgbcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgbequ "F_FUNC(dgbequ,DGBEQU)"(int *m, int *n, int *kl, int *ku, d *ab, int *ldab, d *r, d *c, d *rowcnd, d *colcnd, d *amax, int *info) nogil
cdef dgbequ_t *dgbequ_f = &_fortran_dgbequ

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgbrfs "F_FUNC(dgbrfs,DGBRFS)"(char *trans, int *n, int *kl, int *ku, int *nrhs, d *ab, int *ldab, d *afb, int *ldafb, int *ipiv, d *b, int *ldb, d *x, int *ldx, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dgbrfs_t *dgbrfs_f = &_fortran_dgbrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgbsv "F_FUNC(dgbsv,DGBSV)"(int *n, int *kl, int *ku, int *nrhs, d *ab, int *ldab, int *ipiv, d *b, int *ldb, int *info) nogil
cdef dgbsv_t *dgbsv_f = &_fortran_dgbsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgbsvx "F_FUNC(dgbsvx,DGBSVX)"(char *fact, char *trans, int *n, int *kl, int *ku, int *nrhs, d *ab, int *ldab, d *afb, int *ldafb, int *ipiv, char *equed, d *r, d *c, d *b, int *ldb, d *x, int *ldx, d *rcond, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dgbsvx_t *dgbsvx_f = &_fortran_dgbsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgbtf2 "F_FUNC(dgbtf2,DGBTF2)"(int *m, int *n, int *kl, int *ku, d *ab, int *ldab, int *ipiv, int *info) nogil
cdef dgbtf2_t *dgbtf2_f = &_fortran_dgbtf2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgbtrf "F_FUNC(dgbtrf,DGBTRF)"(int *m, int *n, int *kl, int *ku, d *ab, int *ldab, int *ipiv, int *info) nogil
cdef dgbtrf_t *dgbtrf_f = &_fortran_dgbtrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgbtrs "F_FUNC(dgbtrs,DGBTRS)"(char *trans, int *n, int *kl, int *ku, int *nrhs, d *ab, int *ldab, int *ipiv, d *b, int *ldb, int *info) nogil
cdef dgbtrs_t *dgbtrs_f = &_fortran_dgbtrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgebak "F_FUNC(dgebak,DGEBAK)"(char *job, char *side, int *n, int *ilo, int *ihi, d *scale, int *m, d *v, int *ldv, int *info) nogil
cdef dgebak_t *dgebak_f = &_fortran_dgebak

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgebal "F_FUNC(dgebal,DGEBAL)"(char *job, int *n, d *a, int *lda, int *ilo, int *ihi, d *scale, int *info) nogil
cdef dgebal_t *dgebal_f = &_fortran_dgebal

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgebd2 "F_FUNC(dgebd2,DGEBD2)"(int *m, int *n, d *a, int *lda, d *d, d *e, d *tauq, d *taup, d *work, int *info) nogil
cdef dgebd2_t *dgebd2_f = &_fortran_dgebd2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgebrd "F_FUNC(dgebrd,DGEBRD)"(int *m, int *n, d *a, int *lda, d *d, d *e, d *tauq, d *taup, d *work, int *lwork, int *info) nogil
cdef dgebrd_t *dgebrd_f = &_fortran_dgebrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgecon "F_FUNC(dgecon,DGECON)"(char *norm, int *n, d *a, int *lda, d *anorm, d *rcond, d *work, int *iwork, int *info) nogil
cdef dgecon_t *dgecon_f = &_fortran_dgecon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgeequ "F_FUNC(dgeequ,DGEEQU)"(int *m, int *n, d *a, int *lda, d *r, d *c, d *rowcnd, d *colcnd, d *amax, int *info) nogil
cdef dgeequ_t *dgeequ_f = &_fortran_dgeequ

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgees "F_FUNC(dgees,DGEES)"(char *jobvs, char *sort, dselect2 *select, int *n, d *a, int *lda, int *sdim, d *wr, d *wi, d *vs, int *ldvs, d *work, int *lwork, bint *bwork, int *info) nogil
cdef dgees_t *dgees_f = &_fortran_dgees

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgeesx "F_FUNC(dgeesx,DGEESX)"(char *jobvs, char *sort, dselect2 *select, char *sense, int *n, d *a, int *lda, int *sdim, d *wr, d *wi, d *vs, int *ldvs, d *rconde, d *rcondv, d *work, int *lwork, int *iwork, int *liwork, bint *bwork, int *info) nogil
cdef dgeesx_t *dgeesx_f = &_fortran_dgeesx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgeev "F_FUNC(dgeev,DGEEV)"(char *jobvl, char *jobvr, int *n, d *a, int *lda, d *wr, d *wi, d *vl, int *ldvl, d *vr, int *ldvr, d *work, int *lwork, int *info) nogil
cdef dgeev_t *dgeev_f = &_fortran_dgeev

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgeevx "F_FUNC(dgeevx,DGEEVX)"(char *balanc, char *jobvl, char *jobvr, char *sense, int *n, d *a, int *lda, d *wr, d *wi, d *vl, int *ldvl, d *vr, int *ldvr, int *ilo, int *ihi, d *scale, d *abnrm, d *rconde, d *rcondv, d *work, int *lwork, int *iwork, int *info) nogil
cdef dgeevx_t *dgeevx_f = &_fortran_dgeevx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgegs "F_FUNC(dgegs,DGEGS)"(char *jobvsl, char *jobvsr, int *n, d *a, int *lda, d *b, int *ldb, d *alphar, d *alphai, d *beta, d *vsl, int *ldvsl, d *vsr, int *ldvsr, d *work, int *lwork, int *info) nogil
cdef dgegs_t *dgegs_f = &_fortran_dgegs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgegv "F_FUNC(dgegv,DGEGV)"(char *jobvl, char *jobvr, int *n, d *a, int *lda, d *b, int *ldb, d *alphar, d *alphai, d *beta, d *vl, int *ldvl, d *vr, int *ldvr, d *work, int *lwork, int *info) nogil
cdef dgegv_t *dgegv_f = &_fortran_dgegv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgehd2 "F_FUNC(dgehd2,DGEHD2)"(int *n, int *ilo, int *ihi, d *a, int *lda, d *tau, d *work, int *info) nogil
cdef dgehd2_t *dgehd2_f = &_fortran_dgehd2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgehrd "F_FUNC(dgehrd,DGEHRD)"(int *n, int *ilo, int *ihi, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dgehrd_t *dgehrd_f = &_fortran_dgehrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgelq2 "F_FUNC(dgelq2,DGELQ2)"(int *m, int *n, d *a, int *lda, d *tau, d *work, int *info) nogil
cdef dgelq2_t *dgelq2_f = &_fortran_dgelq2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgelqf "F_FUNC(dgelqf,DGELQF)"(int *m, int *n, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dgelqf_t *dgelqf_f = &_fortran_dgelqf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgels "F_FUNC(dgels,DGELS)"(char *trans, int *m, int *n, int *nrhs, d *a, int *lda, d *b, int *ldb, d *work, int *lwork, int *info) nogil
cdef dgels_t *dgels_f = &_fortran_dgels

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgelsd "F_FUNC(dgelsd,DGELSD)"(int *m, int *n, int *nrhs, d *a, int *lda, d *b, int *ldb, d *s, d *rcond, int *rank, d *work, int *lwork, int *iwork, int *info) nogil
cdef dgelsd_t *dgelsd_f = &_fortran_dgelsd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgelss "F_FUNC(dgelss,DGELSS)"(int *m, int *n, int *nrhs, d *a, int *lda, d *b, int *ldb, d *s, d *rcond, int *rank, d *work, int *lwork, int *info) nogil
cdef dgelss_t *dgelss_f = &_fortran_dgelss

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgelsx "F_FUNC(dgelsx,DGELSX)"(int *m, int *n, int *nrhs, d *a, int *lda, d *b, int *ldb, int *jpvt, d *rcond, int *rank, d *work, int *info) nogil
cdef dgelsx_t *dgelsx_f = &_fortran_dgelsx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgelsy "F_FUNC(dgelsy,DGELSY)"(int *m, int *n, int *nrhs, d *a, int *lda, d *b, int *ldb, int *jpvt, d *rcond, int *rank, d *work, int *lwork, int *info) nogil
cdef dgelsy_t *dgelsy_f = &_fortran_dgelsy

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgeql2 "F_FUNC(dgeql2,DGEQL2)"(int *m, int *n, d *a, int *lda, d *tau, d *work, int *info) nogil
cdef dgeql2_t *dgeql2_f = &_fortran_dgeql2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgeqlf "F_FUNC(dgeqlf,DGEQLF)"(int *m, int *n, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dgeqlf_t *dgeqlf_f = &_fortran_dgeqlf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgeqp3 "F_FUNC(dgeqp3,DGEQP3)"(int *m, int *n, d *a, int *lda, int *jpvt, d *tau, d *work, int *lwork, int *info) nogil
cdef dgeqp3_t *dgeqp3_f = &_fortran_dgeqp3

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgeqpf "F_FUNC(dgeqpf,DGEQPF)"(int *m, int *n, d *a, int *lda, int *jpvt, d *tau, d *work, int *info) nogil
cdef dgeqpf_t *dgeqpf_f = &_fortran_dgeqpf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgeqr2 "F_FUNC(dgeqr2,DGEQR2)"(int *m, int *n, d *a, int *lda, d *tau, d *work, int *info) nogil
cdef dgeqr2_t *dgeqr2_f = &_fortran_dgeqr2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgeqrf "F_FUNC(dgeqrf,DGEQRF)"(int *m, int *n, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dgeqrf_t *dgeqrf_f = &_fortran_dgeqrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgerfs "F_FUNC(dgerfs,DGERFS)"(char *trans, int *n, int *nrhs, d *a, int *lda, d *af, int *ldaf, int *ipiv, d *b, int *ldb, d *x, int *ldx, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dgerfs_t *dgerfs_f = &_fortran_dgerfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgerq2 "F_FUNC(dgerq2,DGERQ2)"(int *m, int *n, d *a, int *lda, d *tau, d *work, int *info) nogil
cdef dgerq2_t *dgerq2_f = &_fortran_dgerq2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgerqf "F_FUNC(dgerqf,DGERQF)"(int *m, int *n, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dgerqf_t *dgerqf_f = &_fortran_dgerqf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgesc2 "F_FUNC(dgesc2,DGESC2)"(int *n, d *a, int *lda, d *rhs, int *ipiv, int *jpiv, d *scale) nogil
cdef dgesc2_t *dgesc2_f = &_fortran_dgesc2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgesdd "F_FUNC(dgesdd,DGESDD)"(char *jobz, int *m, int *n, d *a, int *lda, d *s, d *u, int *ldu, d *vt, int *ldvt, d *work, int *lwork, int *iwork, int *info) nogil
cdef dgesdd_t *dgesdd_f = &_fortran_dgesdd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgesv "F_FUNC(dgesv,DGESV)"(int *n, int *nrhs, d *a, int *lda, int *ipiv, d *b, int *ldb, int *info) nogil
cdef dgesv_t *dgesv_f = &_fortran_dgesv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgesvd "F_FUNC(dgesvd,DGESVD)"(char *jobu, char *jobvt, int *m, int *n, d *a, int *lda, d *s, d *u, int *ldu, d *vt, int *ldvt, d *work, int *lwork, int *info) nogil
cdef dgesvd_t *dgesvd_f = &_fortran_dgesvd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgesvx "F_FUNC(dgesvx,DGESVX)"(char *fact, char *trans, int *n, int *nrhs, d *a, int *lda, d *af, int *ldaf, int *ipiv, char *equed, d *r, d *c, d *b, int *ldb, d *x, int *ldx, d *rcond, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dgesvx_t *dgesvx_f = &_fortran_dgesvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgetc2 "F_FUNC(dgetc2,DGETC2)"(int *n, d *a, int *lda, int *ipiv, int *jpiv, int *info) nogil
cdef dgetc2_t *dgetc2_f = &_fortran_dgetc2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgetf2 "F_FUNC(dgetf2,DGETF2)"(int *m, int *n, d *a, int *lda, int *ipiv, int *info) nogil
cdef dgetf2_t *dgetf2_f = &_fortran_dgetf2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgetrf "F_FUNC(dgetrf,DGETRF)"(int *m, int *n, d *a, int *lda, int *ipiv, int *info) nogil
cdef dgetrf_t *dgetrf_f = &_fortran_dgetrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgetri "F_FUNC(dgetri,DGETRI)"(int *n, d *a, int *lda, int *ipiv, d *work, int *lwork, int *info) nogil
cdef dgetri_t *dgetri_f = &_fortran_dgetri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgetrs "F_FUNC(dgetrs,DGETRS)"(char *trans, int *n, int *nrhs, d *a, int *lda, int *ipiv, d *b, int *ldb, int *info) nogil
cdef dgetrs_t *dgetrs_f = &_fortran_dgetrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dggbak "F_FUNC(dggbak,DGGBAK)"(char *job, char *side, int *n, int *ilo, int *ihi, d *lscale, d *rscale, int *m, d *v, int *ldv, int *info) nogil
cdef dggbak_t *dggbak_f = &_fortran_dggbak

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dggbal "F_FUNC(dggbal,DGGBAL)"(char *job, int *n, d *a, int *lda, d *b, int *ldb, int *ilo, int *ihi, d *lscale, d *rscale, d *work, int *info) nogil
cdef dggbal_t *dggbal_f = &_fortran_dggbal

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgges "F_FUNC(dgges,DGGES)"(char *jobvsl, char *jobvsr, char *sort, dselect3 *selctg, int *n, d *a, int *lda, d *b, int *ldb, int *sdim, d *alphar, d *alphai, d *beta, d *vsl, int *ldvsl, d *vsr, int *ldvsr, d *work, int *lwork, bint *bwork, int *info) nogil
cdef dgges_t *dgges_f = &_fortran_dgges

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dggesx "F_FUNC(dggesx,DGGESX)"(char *jobvsl, char *jobvsr, char *sort, dselect3 *selctg, char *sense, int *n, d *a, int *lda, d *b, int *ldb, int *sdim, d *alphar, d *alphai, d *beta, d *vsl, int *ldvsl, d *vsr, int *ldvsr, d *rconde, d *rcondv, d *work, int *lwork, int *iwork, int *liwork, bint *bwork, int *info) nogil
cdef dggesx_t *dggesx_f = &_fortran_dggesx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dggev "F_FUNC(dggev,DGGEV)"(char *jobvl, char *jobvr, int *n, d *a, int *lda, d *b, int *ldb, d *alphar, d *alphai, d *beta, d *vl, int *ldvl, d *vr, int *ldvr, d *work, int *lwork, int *info) nogil
cdef dggev_t *dggev_f = &_fortran_dggev

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dggevx "F_FUNC(dggevx,DGGEVX)"(char *balanc, char *jobvl, char *jobvr, char *sense, int *n, d *a, int *lda, d *b, int *ldb, d *alphar, d *alphai, d *beta, d *vl, int *ldvl, d *vr, int *ldvr, int *ilo, int *ihi, d *lscale, d *rscale, d *abnrm, d *bbnrm, d *rconde, d *rcondv, d *work, int *lwork, int *iwork, bint *bwork, int *info) nogil
cdef dggevx_t *dggevx_f = &_fortran_dggevx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dggglm "F_FUNC(dggglm,DGGGLM)"(int *n, int *m, int *p, d *a, int *lda, d *b, int *ldb, d *d, d *x, d *y, d *work, int *lwork, int *info) nogil
cdef dggglm_t *dggglm_f = &_fortran_dggglm

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgghrd "F_FUNC(dgghrd,DGGHRD)"(char *compq, char *compz, int *n, int *ilo, int *ihi, d *a, int *lda, d *b, int *ldb, d *q, int *ldq, d *z, int *ldz, int *info) nogil
cdef dgghrd_t *dgghrd_f = &_fortran_dgghrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgglse "F_FUNC(dgglse,DGGLSE)"(int *m, int *n, int *p, d *a, int *lda, d *b, int *ldb, d *c, d *d, d *x, d *work, int *lwork, int *info) nogil
cdef dgglse_t *dgglse_f = &_fortran_dgglse

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dggqrf "F_FUNC(dggqrf,DGGQRF)"(int *n, int *m, int *p, d *a, int *lda, d *taua, d *b, int *ldb, d *taub, d *work, int *lwork, int *info) nogil
cdef dggqrf_t *dggqrf_f = &_fortran_dggqrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dggrqf "F_FUNC(dggrqf,DGGRQF)"(int *m, int *p, int *n, d *a, int *lda, d *taua, d *b, int *ldb, d *taub, d *work, int *lwork, int *info) nogil
cdef dggrqf_t *dggrqf_f = &_fortran_dggrqf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dggsvd "F_FUNC(dggsvd,DGGSVD)"(char *jobu, char *jobv, char *jobq, int *m, int *n, int *p, int *k, int *l, d *a, int *lda, d *b, int *ldb, d *alpha, d *beta, d *u, int *ldu, d *v, int *ldv, d *q, int *ldq, d *work, int *iwork, int *info) nogil
cdef dggsvd_t *dggsvd_f = &_fortran_dggsvd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dggsvp "F_FUNC(dggsvp,DGGSVP)"(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, d *a, int *lda, d *b, int *ldb, d *tola, d *tolb, int *k, int *l, d *u, int *ldu, d *v, int *ldv, d *q, int *ldq, int *iwork, d *tau, d *work, int *info) nogil
cdef dggsvp_t *dggsvp_f = &_fortran_dggsvp

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgtcon "F_FUNC(dgtcon,DGTCON)"(char *norm, int *n, d *dl, d *d, d *du, d *du2, int *ipiv, d *anorm, d *rcond, d *work, int *iwork, int *info) nogil
cdef dgtcon_t *dgtcon_f = &_fortran_dgtcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgtrfs "F_FUNC(dgtrfs,DGTRFS)"(char *trans, int *n, int *nrhs, d *dl, d *d, d *du, d *dlf, d *df, d *duf, d *du2, int *ipiv, d *b, int *ldb, d *x, int *ldx, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dgtrfs_t *dgtrfs_f = &_fortran_dgtrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgtsv "F_FUNC(dgtsv,DGTSV)"(int *n, int *nrhs, d *dl, d *d, d *du, d *b, int *ldb, int *info) nogil
cdef dgtsv_t *dgtsv_f = &_fortran_dgtsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgtsvx "F_FUNC(dgtsvx,DGTSVX)"(char *fact, char *trans, int *n, int *nrhs, d *dl, d *d, d *du, d *dlf, d *df, d *duf, d *du2, int *ipiv, d *b, int *ldb, d *x, int *ldx, d *rcond, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dgtsvx_t *dgtsvx_f = &_fortran_dgtsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgttrf "F_FUNC(dgttrf,DGTTRF)"(int *n, d *dl, d *d, d *du, d *du2, int *ipiv, int *info) nogil
cdef dgttrf_t *dgttrf_f = &_fortran_dgttrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgttrs "F_FUNC(dgttrs,DGTTRS)"(char *trans, int *n, int *nrhs, d *dl, d *d, d *du, d *du2, int *ipiv, d *b, int *ldb, int *info) nogil
cdef dgttrs_t *dgttrs_f = &_fortran_dgttrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dgtts2 "F_FUNC(dgtts2,DGTTS2)"(int *itrans, int *n, int *nrhs, d *dl, d *d, d *du, d *du2, int *ipiv, d *b, int *ldb) nogil
cdef dgtts2_t *dgtts2_f = &_fortran_dgtts2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dhgeqz "F_FUNC(dhgeqz,DHGEQZ)"(char *job, char *compq, char *compz, int *n, int *ilo, int *ihi, d *h, int *ldh, d *t, int *ldt, d *alphar, d *alphai, d *beta, d *q, int *ldq, d *z, int *ldz, d *work, int *lwork, int *info) nogil
cdef dhgeqz_t *dhgeqz_f = &_fortran_dhgeqz

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dhsein "F_FUNC(dhsein,DHSEIN)"(char *side, char *eigsrc, char *initv, bint *select, int *n, d *h, int *ldh, d *wr, d *wi, d *vl, int *ldvl, d *vr, int *ldvr, int *mm, int *m, d *work, int *ifaill, int *ifailr, int *info) nogil
cdef dhsein_t *dhsein_f = &_fortran_dhsein

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dhseqr "F_FUNC(dhseqr,DHSEQR)"(char *job, char *compz, int *n, int *ilo, int *ihi, d *h, int *ldh, d *wr, d *wi, d *z, int *ldz, d *work, int *lwork, int *info) nogil
cdef dhseqr_t *dhseqr_f = &_fortran_dhseqr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dlacn2 "F_FUNC(dlacn2,DLACN2)"(int *n, d *v, d *x, int *isgn, d *est, int *kase, int *isave) nogil
cdef dlacn2_t *dlacn2_f = &_fortran_dlacn2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dlacon "F_FUNC(dlacon,DLACON)"(int *n, d *v, d *x, int *isgn, d *est, int *kase) nogil
cdef dlacon_t *dlacon_f = &_fortran_dlacon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dlanv2 "F_FUNC(dlanv2,DLANV2)"(d *a, d *b, d *c, d *d, d *rt1r, d *rt1i, d *rt2r, d *rt2i, d *cs, d *sn) nogil
cdef dlanv2_t *dlanv2_f = &_fortran_dlanv2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dlarf "F_FUNC(dlarf,DLARF)"(char *side, int *m, int *n, d *v, int *incv, d *tau, d *c, int *ldc, d *work) nogil
cdef dlarf_t *dlarf_f = &_fortran_dlarf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dlarz "F_FUNC(dlarz,DLARZ)"(char *side, int *m, int *n, int *l, d *v, int *incv, d *tau, d *c, int *ldc, d *work) nogil
cdef dlarz_t *dlarz_f = &_fortran_dlarz

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dlaswp "F_FUNC(dlaswp,DLASWP)"(int *n, d *a, int *lda, int *k1, int *k2, int *ipiv, int *incx) nogil
cdef dlaswp_t *dlaswp_f = &_fortran_dlaswp

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dlauum "F_FUNC(dlauum,DLAUUM)"(char *uplo, int *n, d *a, int *lda, int *info) nogil
cdef dlauum_t *dlauum_f = &_fortran_dlauum

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dopgtr "F_FUNC(dopgtr,DOPGTR)"(char *uplo, int *n, d *ap, d *tau, d *q, int *ldq, d *work, int *info) nogil
cdef dopgtr_t *dopgtr_f = &_fortran_dopgtr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dopmtr "F_FUNC(dopmtr,DOPMTR)"(char *side, char *uplo, char *trans, int *m, int *n, d *ap, d *tau, d *c, int *ldc, d *work, int *info) nogil
cdef dopmtr_t *dopmtr_f = &_fortran_dopmtr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dorg2l "F_FUNC(dorg2l,DORG2L)"(int *m, int *n, int *k, d *a, int *lda, d *tau, d *work, int *info) nogil
cdef dorg2l_t *dorg2l_f = &_fortran_dorg2l

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dorg2r "F_FUNC(dorg2r,DORG2R)"(int *m, int *n, int *k, d *a, int *lda, d *tau, d *work, int *info) nogil
cdef dorg2r_t *dorg2r_f = &_fortran_dorg2r

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dorgbr "F_FUNC(dorgbr,DORGBR)"(char *vect, int *m, int *n, int *k, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dorgbr_t *dorgbr_f = &_fortran_dorgbr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dorghr "F_FUNC(dorghr,DORGHR)"(int *n, int *ilo, int *ihi, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dorghr_t *dorghr_f = &_fortran_dorghr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dorgl2 "F_FUNC(dorgl2,DORGL2)"(int *m, int *n, int *k, d *a, int *lda, d *tau, d *work, int *info) nogil
cdef dorgl2_t *dorgl2_f = &_fortran_dorgl2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dorglq "F_FUNC(dorglq,DORGLQ)"(int *m, int *n, int *k, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dorglq_t *dorglq_f = &_fortran_dorglq

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dorgql "F_FUNC(dorgql,DORGQL)"(int *m, int *n, int *k, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dorgql_t *dorgql_f = &_fortran_dorgql

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dorgqr "F_FUNC(dorgqr,DORGQR)"(int *m, int *n, int *k, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dorgqr_t *dorgqr_f = &_fortran_dorgqr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dorgr2 "F_FUNC(dorgr2,DORGR2)"(int *m, int *n, int *k, d *a, int *lda, d *tau, d *work, int *info) nogil
cdef dorgr2_t *dorgr2_f = &_fortran_dorgr2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dorgrq "F_FUNC(dorgrq,DORGRQ)"(int *m, int *n, int *k, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dorgrq_t *dorgrq_f = &_fortran_dorgrq

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dorgtr "F_FUNC(dorgtr,DORGTR)"(char *uplo, int *n, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dorgtr_t *dorgtr_f = &_fortran_dorgtr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dorm2l "F_FUNC(dorm2l,DORM2L)"(char *side, char *trans, int *m, int *n, int *k, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *info) nogil
cdef dorm2l_t *dorm2l_f = &_fortran_dorm2l

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dorm2r "F_FUNC(dorm2r,DORM2R)"(char *side, char *trans, int *m, int *n, int *k, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *info) nogil
cdef dorm2r_t *dorm2r_f = &_fortran_dorm2r

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dormbr "F_FUNC(dormbr,DORMBR)"(char *vect, char *side, char *trans, int *m, int *n, int *k, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *lwork, int *info) nogil
cdef dormbr_t *dormbr_f = &_fortran_dormbr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dormhr "F_FUNC(dormhr,DORMHR)"(char *side, char *trans, int *m, int *n, int *ilo, int *ihi, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *lwork, int *info) nogil
cdef dormhr_t *dormhr_f = &_fortran_dormhr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dorml2 "F_FUNC(dorml2,DORML2)"(char *side, char *trans, int *m, int *n, int *k, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *info) nogil
cdef dorml2_t *dorml2_f = &_fortran_dorml2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dormlq "F_FUNC(dormlq,DORMLQ)"(char *side, char *trans, int *m, int *n, int *k, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *lwork, int *info) nogil
cdef dormlq_t *dormlq_f = &_fortran_dormlq

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dormql "F_FUNC(dormql,DORMQL)"(char *side, char *trans, int *m, int *n, int *k, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *lwork, int *info) nogil
cdef dormql_t *dormql_f = &_fortran_dormql

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dormqr "F_FUNC(dormqr,DORMQR)"(char *side, char *trans, int *m, int *n, int *k, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *lwork, int *info) nogil
cdef dormqr_t *dormqr_f = &_fortran_dormqr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dormr2 "F_FUNC(dormr2,DORMR2)"(char *side, char *trans, int *m, int *n, int *k, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *info) nogil
cdef dormr2_t *dormr2_f = &_fortran_dormr2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dormr3 "F_FUNC(dormr3,DORMR3)"(char *side, char *trans, int *m, int *n, int *k, int *l, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *info) nogil
cdef dormr3_t *dormr3_f = &_fortran_dormr3

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dormrq "F_FUNC(dormrq,DORMRQ)"(char *side, char *trans, int *m, int *n, int *k, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *lwork, int *info) nogil
cdef dormrq_t *dormrq_f = &_fortran_dormrq

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dormrz "F_FUNC(dormrz,DORMRZ)"(char *side, char *trans, int *m, int *n, int *k, int *l, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *lwork, int *info) nogil
cdef dormrz_t *dormrz_f = &_fortran_dormrz

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dormtr "F_FUNC(dormtr,DORMTR)"(char *side, char *uplo, char *trans, int *m, int *n, d *a, int *lda, d *tau, d *c, int *ldc, d *work, int *lwork, int *info) nogil
cdef dormtr_t *dormtr_f = &_fortran_dormtr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dpbcon "F_FUNC(dpbcon,DPBCON)"(char *uplo, int *n, int *kd, d *ab, int *ldab, d *anorm, d *rcond, d *work, int *iwork, int *info) nogil
cdef dpbcon_t *dpbcon_f = &_fortran_dpbcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dpbequ "F_FUNC(dpbequ,DPBEQU)"(char *uplo, int *n, int *kd, d *ab, int *ldab, d *s, d *scond, d *amax, int *info) nogil
cdef dpbequ_t *dpbequ_f = &_fortran_dpbequ

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dpbrfs "F_FUNC(dpbrfs,DPBRFS)"(char *uplo, int *n, int *kd, int *nrhs, d *ab, int *ldab, d *afb, int *ldafb, d *b, int *ldb, d *x, int *ldx, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dpbrfs_t *dpbrfs_f = &_fortran_dpbrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dpbstf "F_FUNC(dpbstf,DPBSTF)"(char *uplo, int *n, int *kd, d *ab, int *ldab, int *info) nogil
cdef dpbstf_t *dpbstf_f = &_fortran_dpbstf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dpbsv "F_FUNC(dpbsv,DPBSV)"(char *uplo, int *n, int *kd, int *nrhs, d *ab, int *ldab, d *b, int *ldb, int *info) nogil
cdef dpbsv_t *dpbsv_f = &_fortran_dpbsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dpbsvx "F_FUNC(dpbsvx,DPBSVX)"(char *fact, char *uplo, int *n, int *kd, int *nrhs, d *ab, int *ldab, d *afb, int *ldafb, char *equed, d *s, d *b, int *ldb, d *x, int *ldx, d *rcond, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dpbsvx_t *dpbsvx_f = &_fortran_dpbsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dpbtf2 "F_FUNC(dpbtf2,DPBTF2)"(char *uplo, int *n, int *kd, d *ab, int *ldab, int *info) nogil
cdef dpbtf2_t *dpbtf2_f = &_fortran_dpbtf2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dpbtrf "F_FUNC(dpbtrf,DPBTRF)"(char *uplo, int *n, int *kd, d *ab, int *ldab, int *info) nogil
cdef dpbtrf_t *dpbtrf_f = &_fortran_dpbtrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dpbtrs "F_FUNC(dpbtrs,DPBTRS)"(char *uplo, int *n, int *kd, int *nrhs, d *ab, int *ldab, d *b, int *ldb, int *info) nogil
cdef dpbtrs_t *dpbtrs_f = &_fortran_dpbtrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dpocon "F_FUNC(dpocon,DPOCON)"(char *uplo, int *n, d *a, int *lda, d *anorm, d *rcond, d *work, int *iwork, int *info) nogil
cdef dpocon_t *dpocon_f = &_fortran_dpocon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dpoequ "F_FUNC(dpoequ,DPOEQU)"(int *n, d *a, int *lda, d *s, d *scond, d *amax, int *info) nogil
cdef dpoequ_t *dpoequ_f = &_fortran_dpoequ

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dporfs "F_FUNC(dporfs,DPORFS)"(char *uplo, int *n, int *nrhs, d *a, int *lda, d *af, int *ldaf, d *b, int *ldb, d *x, int *ldx, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dporfs_t *dporfs_f = &_fortran_dporfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dposv "F_FUNC(dposv,DPOSV)"(char *uplo, int *n, int *nrhs, d *a, int *lda, d *b, int *ldb, int *info) nogil
cdef dposv_t *dposv_f = &_fortran_dposv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dposvx "F_FUNC(dposvx,DPOSVX)"(char *fact, char *uplo, int *n, int *nrhs, d *a, int *lda, d *af, int *ldaf, char *equed, d *s, d *b, int *ldb, d *x, int *ldx, d *rcond, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dposvx_t *dposvx_f = &_fortran_dposvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dpotf2 "F_FUNC(dpotf2,DPOTF2)"(char *uplo, int *n, d *a, int *lda, int *info) nogil
cdef dpotf2_t *dpotf2_f = &_fortran_dpotf2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dpotrf "F_FUNC(dpotrf,DPOTRF)"(char *uplo, int *n, d *a, int *lda, int *info) nogil
cdef dpotrf_t *dpotrf_f = &_fortran_dpotrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dpotri "F_FUNC(dpotri,DPOTRI)"(char *uplo, int *n, d *a, int *lda, int *info) nogil
cdef dpotri_t *dpotri_f = &_fortran_dpotri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dpotrs "F_FUNC(dpotrs,DPOTRS)"(char *uplo, int *n, int *nrhs, d *a, int *lda, d *b, int *ldb, int *info) nogil
cdef dpotrs_t *dpotrs_f = &_fortran_dpotrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dppcon "F_FUNC(dppcon,DPPCON)"(char *uplo, int *n, d *ap, d *anorm, d *rcond, d *work, int *iwork, int *info) nogil
cdef dppcon_t *dppcon_f = &_fortran_dppcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dppequ "F_FUNC(dppequ,DPPEQU)"(char *uplo, int *n, d *ap, d *s, d *scond, d *amax, int *info) nogil
cdef dppequ_t *dppequ_f = &_fortran_dppequ

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dpprfs "F_FUNC(dpprfs,DPPRFS)"(char *uplo, int *n, int *nrhs, d *ap, d *afp, d *b, int *ldb, d *x, int *ldx, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dpprfs_t *dpprfs_f = &_fortran_dpprfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dppsv "F_FUNC(dppsv,DPPSV)"(char *uplo, int *n, int *nrhs, d *ap, d *b, int *ldb, int *info) nogil
cdef dppsv_t *dppsv_f = &_fortran_dppsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dppsvx "F_FUNC(dppsvx,DPPSVX)"(char *fact, char *uplo, int *n, int *nrhs, d *ap, d *afp, char *equed, d *s, d *b, int *ldb, d *x, int *ldx, d *rcond, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dppsvx_t *dppsvx_f = &_fortran_dppsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dpptrf "F_FUNC(dpptrf,DPPTRF)"(char *uplo, int *n, d *ap, int *info) nogil
cdef dpptrf_t *dpptrf_f = &_fortran_dpptrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dpptri "F_FUNC(dpptri,DPPTRI)"(char *uplo, int *n, d *ap, int *info) nogil
cdef dpptri_t *dpptri_f = &_fortran_dpptri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dpptrs "F_FUNC(dpptrs,DPPTRS)"(char *uplo, int *n, int *nrhs, d *ap, d *b, int *ldb, int *info) nogil
cdef dpptrs_t *dpptrs_f = &_fortran_dpptrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dptcon "F_FUNC(dptcon,DPTCON)"(int *n, d *d, d *e, d *anorm, d *rcond, d *work, int *info) nogil
cdef dptcon_t *dptcon_f = &_fortran_dptcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dpteqr "F_FUNC(dpteqr,DPTEQR)"(char *compz, int *n, d *d, d *e, d *z, int *ldz, d *work, int *info) nogil
cdef dpteqr_t *dpteqr_f = &_fortran_dpteqr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dptrfs "F_FUNC(dptrfs,DPTRFS)"(int *n, int *nrhs, d *d, d *e, d *df, d *ef, d *b, int *ldb, d *x, int *ldx, d *ferr, d *berr, d *work, int *info) nogil
cdef dptrfs_t *dptrfs_f = &_fortran_dptrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dptsv "F_FUNC(dptsv,DPTSV)"(int *n, int *nrhs, d *d, d *e, d *b, int *ldb, int *info) nogil
cdef dptsv_t *dptsv_f = &_fortran_dptsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dptsvx "F_FUNC(dptsvx,DPTSVX)"(char *fact, int *n, int *nrhs, d *d, d *e, d *df, d *ef, d *b, int *ldb, d *x, int *ldx, d *rcond, d *ferr, d *berr, d *work, int *info) nogil
cdef dptsvx_t *dptsvx_f = &_fortran_dptsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dpttrf "F_FUNC(dpttrf,DPTTRF)"(int *n, d *d, d *e, int *info) nogil
cdef dpttrf_t *dpttrf_f = &_fortran_dpttrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dpttrs "F_FUNC(dpttrs,DPTTRS)"(int *n, int *nrhs, d *d, d *e, d *b, int *ldb, int *info) nogil
cdef dpttrs_t *dpttrs_f = &_fortran_dpttrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dptts2 "F_FUNC(dptts2,DPTTS2)"(int *n, int *nrhs, d *d, d *e, d *b, int *ldb) nogil
cdef dptts2_t *dptts2_f = &_fortran_dptts2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_drscl "F_FUNC(drscl,DRSCL)"(int *n, d *sa, d *sx, int *incx) nogil
cdef drscl_t *drscl_f = &_fortran_drscl

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsbev "F_FUNC(dsbev,DSBEV)"(char *jobz, char *uplo, int *n, int *kd, d *ab, int *ldab, d *w, d *z, int *ldz, d *work, int *info) nogil
cdef dsbev_t *dsbev_f = &_fortran_dsbev

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsbevd "F_FUNC(dsbevd,DSBEVD)"(char *jobz, char *uplo, int *n, int *kd, d *ab, int *ldab, d *w, d *z, int *ldz, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dsbevd_t *dsbevd_f = &_fortran_dsbevd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsbevx "F_FUNC(dsbevx,DSBEVX)"(char *jobz, char *range, char *uplo, int *n, int *kd, d *ab, int *ldab, d *q, int *ldq, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, d *z, int *ldz, d *work, int *iwork, int *ifail, int *info) nogil
cdef dsbevx_t *dsbevx_f = &_fortran_dsbevx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsbgst "F_FUNC(dsbgst,DSBGST)"(char *vect, char *uplo, int *n, int *ka, int *kb, d *ab, int *ldab, d *bb, int *ldbb, d *x, int *ldx, d *work, int *info) nogil
cdef dsbgst_t *dsbgst_f = &_fortran_dsbgst

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsbgv "F_FUNC(dsbgv,DSBGV)"(char *jobz, char *uplo, int *n, int *ka, int *kb, d *ab, int *ldab, d *bb, int *ldbb, d *w, d *z, int *ldz, d *work, int *info) nogil
cdef dsbgv_t *dsbgv_f = &_fortran_dsbgv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsbgvd "F_FUNC(dsbgvd,DSBGVD)"(char *jobz, char *uplo, int *n, int *ka, int *kb, d *ab, int *ldab, d *bb, int *ldbb, d *w, d *z, int *ldz, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dsbgvd_t *dsbgvd_f = &_fortran_dsbgvd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsbgvx "F_FUNC(dsbgvx,DSBGVX)"(char *jobz, char *range, char *uplo, int *n, int *ka, int *kb, d *ab, int *ldab, d *bb, int *ldbb, d *q, int *ldq, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, d *z, int *ldz, d *work, int *iwork, int *ifail, int *info) nogil
cdef dsbgvx_t *dsbgvx_f = &_fortran_dsbgvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsbtrd "F_FUNC(dsbtrd,DSBTRD)"(char *vect, char *uplo, int *n, int *kd, d *ab, int *ldab, d *d, d *e, d *q, int *ldq, d *work, int *info) nogil
cdef dsbtrd_t *dsbtrd_f = &_fortran_dsbtrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsgesv "F_FUNC(dsgesv,DSGESV)"(int *n, int *nrhs, d *a, int *lda, int *ipiv, d *b, int *ldb, d *x, int *ldx, d *work, s *swork, int *iter, int *info) nogil
cdef dsgesv_t *dsgesv_f = &_fortran_dsgesv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dspcon "F_FUNC(dspcon,DSPCON)"(char *uplo, int *n, d *ap, int *ipiv, d *anorm, d *rcond, d *work, int *iwork, int *info) nogil
cdef dspcon_t *dspcon_f = &_fortran_dspcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dspev "F_FUNC(dspev,DSPEV)"(char *jobz, char *uplo, int *n, d *ap, d *w, d *z, int *ldz, d *work, int *info) nogil
cdef dspev_t *dspev_f = &_fortran_dspev

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dspevd "F_FUNC(dspevd,DSPEVD)"(char *jobz, char *uplo, int *n, d *ap, d *w, d *z, int *ldz, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dspevd_t *dspevd_f = &_fortran_dspevd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dspevx "F_FUNC(dspevx,DSPEVX)"(char *jobz, char *range, char *uplo, int *n, d *ap, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, d *z, int *ldz, d *work, int *iwork, int *ifail, int *info) nogil
cdef dspevx_t *dspevx_f = &_fortran_dspevx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dspgst "F_FUNC(dspgst,DSPGST)"(int *itype, char *uplo, int *n, d *ap, d *bp, int *info) nogil
cdef dspgst_t *dspgst_f = &_fortran_dspgst

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dspgv "F_FUNC(dspgv,DSPGV)"(int *itype, char *jobz, char *uplo, int *n, d *ap, d *bp, d *w, d *z, int *ldz, d *work, int *info) nogil
cdef dspgv_t *dspgv_f = &_fortran_dspgv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dspgvd "F_FUNC(dspgvd,DSPGVD)"(int *itype, char *jobz, char *uplo, int *n, d *ap, d *bp, d *w, d *z, int *ldz, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dspgvd_t *dspgvd_f = &_fortran_dspgvd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dspgvx "F_FUNC(dspgvx,DSPGVX)"(int *itype, char *jobz, char *range, char *uplo, int *n, d *ap, d *bp, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, d *z, int *ldz, d *work, int *iwork, int *ifail, int *info) nogil
cdef dspgvx_t *dspgvx_f = &_fortran_dspgvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsprfs "F_FUNC(dsprfs,DSPRFS)"(char *uplo, int *n, int *nrhs, d *ap, d *afp, int *ipiv, d *b, int *ldb, d *x, int *ldx, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dsprfs_t *dsprfs_f = &_fortran_dsprfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dspsv "F_FUNC(dspsv,DSPSV)"(char *uplo, int *n, int *nrhs, d *ap, int *ipiv, d *b, int *ldb, int *info) nogil
cdef dspsv_t *dspsv_f = &_fortran_dspsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dspsvx "F_FUNC(dspsvx,DSPSVX)"(char *fact, char *uplo, int *n, int *nrhs, d *ap, d *afp, int *ipiv, d *b, int *ldb, d *x, int *ldx, d *rcond, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dspsvx_t *dspsvx_f = &_fortran_dspsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsptrd "F_FUNC(dsptrd,DSPTRD)"(char *uplo, int *n, d *ap, d *d, d *e, d *tau, int *info) nogil
cdef dsptrd_t *dsptrd_f = &_fortran_dsptrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsptrf "F_FUNC(dsptrf,DSPTRF)"(char *uplo, int *n, d *ap, int *ipiv, int *info) nogil
cdef dsptrf_t *dsptrf_f = &_fortran_dsptrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsptri "F_FUNC(dsptri,DSPTRI)"(char *uplo, int *n, d *ap, int *ipiv, d *work, int *info) nogil
cdef dsptri_t *dsptri_f = &_fortran_dsptri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsptrs "F_FUNC(dsptrs,DSPTRS)"(char *uplo, int *n, int *nrhs, d *ap, int *ipiv, d *b, int *ldb, int *info) nogil
cdef dsptrs_t *dsptrs_f = &_fortran_dsptrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dstebz "F_FUNC(dstebz,DSTEBZ)"(char *range, char *order, int *n, d *vl, d *vu, int *il, int *iu, d *abstol, d *d, d *e, int *m, int *nsplit, d *w, int *iblock, int *isplit, d *work, int *iwork, int *info) nogil
cdef dstebz_t *dstebz_f = &_fortran_dstebz

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dstedc "F_FUNC(dstedc,DSTEDC)"(char *compz, int *n, d *d, d *e, d *z, int *ldz, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dstedc_t *dstedc_f = &_fortran_dstedc

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dstegr "F_FUNC(dstegr,DSTEGR)"(char *jobz, char *range, int *n, d *d, d *e, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, d *z, int *ldz, int *isuppz, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dstegr_t *dstegr_f = &_fortran_dstegr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dstein "F_FUNC(dstein,DSTEIN)"(int *n, d *d, d *e, int *m, d *w, int *iblock, int *isplit, d *z, int *ldz, d *work, int *iwork, int *ifail, int *info) nogil
cdef dstein_t *dstein_f = &_fortran_dstein

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dstemr "F_FUNC(dstemr,DSTEMR)"(char *jobz, char *range, int *n, d *d, d *e, d *vl, d *vu, int *il, int *iu, int *m, d *w, d *z, int *ldz, int *nzc, int *isuppz, bint *tryrac, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dstemr_t *dstemr_f = &_fortran_dstemr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsteqr "F_FUNC(dsteqr,DSTEQR)"(char *compz, int *n, d *d, d *e, d *z, int *ldz, d *work, int *info) nogil
cdef dsteqr_t *dsteqr_f = &_fortran_dsteqr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsterf "F_FUNC(dsterf,DSTERF)"(int *n, d *d, d *e, int *info) nogil
cdef dsterf_t *dsterf_f = &_fortran_dsterf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dstev "F_FUNC(dstev,DSTEV)"(char *jobz, int *n, d *d, d *e, d *z, int *ldz, d *work, int *info) nogil
cdef dstev_t *dstev_f = &_fortran_dstev

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dstevd "F_FUNC(dstevd,DSTEVD)"(char *jobz, int *n, d *d, d *e, d *z, int *ldz, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dstevd_t *dstevd_f = &_fortran_dstevd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dstevr "F_FUNC(dstevr,DSTEVR)"(char *jobz, char *range, int *n, d *d, d *e, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, d *z, int *ldz, int *isuppz, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dstevr_t *dstevr_f = &_fortran_dstevr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dstevx "F_FUNC(dstevx,DSTEVX)"(char *jobz, char *range, int *n, d *d, d *e, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, d *z, int *ldz, d *work, int *iwork, int *ifail, int *info) nogil
cdef dstevx_t *dstevx_f = &_fortran_dstevx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsycon "F_FUNC(dsycon,DSYCON)"(char *uplo, int *n, d *a, int *lda, int *ipiv, d *anorm, d *rcond, d *work, int *iwork, int *info) nogil
cdef dsycon_t *dsycon_f = &_fortran_dsycon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsyev "F_FUNC(dsyev,DSYEV)"(char *jobz, char *uplo, int *n, d *a, int *lda, d *w, d *work, int *lwork, int *info) nogil
cdef dsyev_t *dsyev_f = &_fortran_dsyev

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsyevd "F_FUNC(dsyevd,DSYEVD)"(char *jobz, char *uplo, int *n, d *a, int *lda, d *w, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dsyevd_t *dsyevd_f = &_fortran_dsyevd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsyevr "F_FUNC(dsyevr,DSYEVR)"(char *jobz, char *range, char *uplo, int *n, d *a, int *lda, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, d *z, int *ldz, int *isuppz, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dsyevr_t *dsyevr_f = &_fortran_dsyevr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsyevx "F_FUNC(dsyevx,DSYEVX)"(char *jobz, char *range, char *uplo, int *n, d *a, int *lda, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, d *z, int *ldz, d *work, int *lwork, int *iwork, int *ifail, int *info) nogil
cdef dsyevx_t *dsyevx_f = &_fortran_dsyevx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsygs2 "F_FUNC(dsygs2,DSYGS2)"(int *itype, char *uplo, int *n, d *a, int *lda, d *b, int *ldb, int *info) nogil
cdef dsygs2_t *dsygs2_f = &_fortran_dsygs2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsygst "F_FUNC(dsygst,DSYGST)"(int *itype, char *uplo, int *n, d *a, int *lda, d *b, int *ldb, int *info) nogil
cdef dsygst_t *dsygst_f = &_fortran_dsygst

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsygv "F_FUNC(dsygv,DSYGV)"(int *itype, char *jobz, char *uplo, int *n, d *a, int *lda, d *b, int *ldb, d *w, d *work, int *lwork, int *info) nogil
cdef dsygv_t *dsygv_f = &_fortran_dsygv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsygvd "F_FUNC(dsygvd,DSYGVD)"(int *itype, char *jobz, char *uplo, int *n, d *a, int *lda, d *b, int *ldb, d *w, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dsygvd_t *dsygvd_f = &_fortran_dsygvd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsygvx "F_FUNC(dsygvx,DSYGVX)"(int *itype, char *jobz, char *range, char *uplo, int *n, d *a, int *lda, d *b, int *ldb, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, d *z, int *ldz, d *work, int *lwork, int *iwork, int *ifail, int *info) nogil
cdef dsygvx_t *dsygvx_f = &_fortran_dsygvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsyrfs "F_FUNC(dsyrfs,DSYRFS)"(char *uplo, int *n, int *nrhs, d *a, int *lda, d *af, int *ldaf, int *ipiv, d *b, int *ldb, d *x, int *ldx, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dsyrfs_t *dsyrfs_f = &_fortran_dsyrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsysv "F_FUNC(dsysv,DSYSV)"(char *uplo, int *n, int *nrhs, d *a, int *lda, int *ipiv, d *b, int *ldb, d *work, int *lwork, int *info) nogil
cdef dsysv_t *dsysv_f = &_fortran_dsysv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsysvx "F_FUNC(dsysvx,DSYSVX)"(char *fact, char *uplo, int *n, int *nrhs, d *a, int *lda, d *af, int *ldaf, int *ipiv, d *b, int *ldb, d *x, int *ldx, d *rcond, d *ferr, d *berr, d *work, int *lwork, int *iwork, int *info) nogil
cdef dsysvx_t *dsysvx_f = &_fortran_dsysvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsytd2 "F_FUNC(dsytd2,DSYTD2)"(char *uplo, int *n, d *a, int *lda, d *d, d *e, d *tau, int *info) nogil
cdef dsytd2_t *dsytd2_f = &_fortran_dsytd2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsytf2 "F_FUNC(dsytf2,DSYTF2)"(char *uplo, int *n, d *a, int *lda, int *ipiv, int *info) nogil
cdef dsytf2_t *dsytf2_f = &_fortran_dsytf2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsytrd "F_FUNC(dsytrd,DSYTRD)"(char *uplo, int *n, d *a, int *lda, d *d, d *e, d *tau, d *work, int *lwork, int *info) nogil
cdef dsytrd_t *dsytrd_f = &_fortran_dsytrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsytrf "F_FUNC(dsytrf,DSYTRF)"(char *uplo, int *n, d *a, int *lda, int *ipiv, d *work, int *lwork, int *info) nogil
cdef dsytrf_t *dsytrf_f = &_fortran_dsytrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsytri "F_FUNC(dsytri,DSYTRI)"(char *uplo, int *n, d *a, int *lda, int *ipiv, d *work, int *info) nogil
cdef dsytri_t *dsytri_f = &_fortran_dsytri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dsytrs "F_FUNC(dsytrs,DSYTRS)"(char *uplo, int *n, int *nrhs, d *a, int *lda, int *ipiv, d *b, int *ldb, int *info) nogil
cdef dsytrs_t *dsytrs_f = &_fortran_dsytrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtbcon "F_FUNC(dtbcon,DTBCON)"(char *norm, char *uplo, char *diag, int *n, int *kd, d *ab, int *ldab, d *rcond, d *work, int *iwork, int *info) nogil
cdef dtbcon_t *dtbcon_f = &_fortran_dtbcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtbrfs "F_FUNC(dtbrfs,DTBRFS)"(char *uplo, char *trans, char *diag, int *n, int *kd, int *nrhs, d *ab, int *ldab, d *b, int *ldb, d *x, int *ldx, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dtbrfs_t *dtbrfs_f = &_fortran_dtbrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtbtrs "F_FUNC(dtbtrs,DTBTRS)"(char *uplo, char *trans, char *diag, int *n, int *kd, int *nrhs, d *ab, int *ldab, d *b, int *ldb, int *info) nogil
cdef dtbtrs_t *dtbtrs_f = &_fortran_dtbtrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtgevc "F_FUNC(dtgevc,DTGEVC)"(char *side, char *howmny, bint *select, int *n, d *s, int *lds, d *p, int *ldp, d *vl, int *ldvl, d *vr, int *ldvr, int *mm, int *m, d *work, int *info) nogil
cdef dtgevc_t *dtgevc_f = &_fortran_dtgevc

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtgex2 "F_FUNC(dtgex2,DTGEX2)"(bint *wantq, bint *wantz, int *n, d *a, int *lda, d *b, int *ldb, d *q, int *ldq, d *z, int *ldz, int *j1, int *n1, int *n2, d *work, int *lwork, int *info) nogil
cdef dtgex2_t *dtgex2_f = &_fortran_dtgex2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtgexc "F_FUNC(dtgexc,DTGEXC)"(bint *wantq, bint *wantz, int *n, d *a, int *lda, d *b, int *ldb, d *q, int *ldq, d *z, int *ldz, int *ifst, int *ilst, d *work, int *lwork, int *info) nogil
cdef dtgexc_t *dtgexc_f = &_fortran_dtgexc

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtgsen "F_FUNC(dtgsen,DTGSEN)"(int *ijob, bint *wantq, bint *wantz, bint *select, int *n, d *a, int *lda, d *b, int *ldb, d *alphar, d *alphai, d *beta, d *q, int *ldq, d *z, int *ldz, int *m, d *pl, d *pr, d *dif, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dtgsen_t *dtgsen_f = &_fortran_dtgsen

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtgsja "F_FUNC(dtgsja,DTGSJA)"(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, int *k, int *l, d *a, int *lda, d *b, int *ldb, d *tola, d *tolb, d *alpha, d *beta, d *u, int *ldu, d *v, int *ldv, d *q, int *ldq, d *work, int *ncycle, int *info) nogil
cdef dtgsja_t *dtgsja_f = &_fortran_dtgsja

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtgsna "F_FUNC(dtgsna,DTGSNA)"(char *job, char *howmny, bint *select, int *n, d *a, int *lda, d *b, int *ldb, d *vl, int *ldvl, d *vr, int *ldvr, d *s, d *dif, int *mm, int *m, d *work, int *lwork, int *iwork, int *info) nogil
cdef dtgsna_t *dtgsna_f = &_fortran_dtgsna

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtgsy2 "F_FUNC(dtgsy2,DTGSY2)"(char *trans, int *ijob, int *m, int *n, d *a, int *lda, d *b, int *ldb, d *c, int *ldc, d *d, int *ldd, d *e, int *lde, d *f, int *ldf, d *scale, d *rdsum, d *rdscal, int *iwork, int *pq, int *info) nogil
cdef dtgsy2_t *dtgsy2_f = &_fortran_dtgsy2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtgsyl "F_FUNC(dtgsyl,DTGSYL)"(char *trans, int *ijob, int *m, int *n, d *a, int *lda, d *b, int *ldb, d *c, int *ldc, d *d, int *ldd, d *e, int *lde, d *f, int *ldf, d *scale, d *dif, d *work, int *lwork, int *iwork, int *info) nogil
cdef dtgsyl_t *dtgsyl_f = &_fortran_dtgsyl

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtpcon "F_FUNC(dtpcon,DTPCON)"(char *norm, char *uplo, char *diag, int *n, d *ap, d *rcond, d *work, int *iwork, int *info) nogil
cdef dtpcon_t *dtpcon_f = &_fortran_dtpcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtprfs "F_FUNC(dtprfs,DTPRFS)"(char *uplo, char *trans, char *diag, int *n, int *nrhs, d *ap, d *b, int *ldb, d *x, int *ldx, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dtprfs_t *dtprfs_f = &_fortran_dtprfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtptri "F_FUNC(dtptri,DTPTRI)"(char *uplo, char *diag, int *n, d *ap, int *info) nogil
cdef dtptri_t *dtptri_f = &_fortran_dtptri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtptrs "F_FUNC(dtptrs,DTPTRS)"(char *uplo, char *trans, char *diag, int *n, int *nrhs, d *ap, d *b, int *ldb, int *info) nogil
cdef dtptrs_t *dtptrs_f = &_fortran_dtptrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtrcon "F_FUNC(dtrcon,DTRCON)"(char *norm, char *uplo, char *diag, int *n, d *a, int *lda, d *rcond, d *work, int *iwork, int *info) nogil
cdef dtrcon_t *dtrcon_f = &_fortran_dtrcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtrevc "F_FUNC(dtrevc,DTREVC)"(char *side, char *howmny, bint *select, int *n, d *t, int *ldt, d *vl, int *ldvl, d *vr, int *ldvr, int *mm, int *m, d *work, int *info) nogil
cdef dtrevc_t *dtrevc_f = &_fortran_dtrevc

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtrexc "F_FUNC(dtrexc,DTREXC)"(char *compq, int *n, d *t, int *ldt, d *q, int *ldq, int *ifst, int *ilst, d *work, int *info) nogil
cdef dtrexc_t *dtrexc_f = &_fortran_dtrexc

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtrrfs "F_FUNC(dtrrfs,DTRRFS)"(char *uplo, char *trans, char *diag, int *n, int *nrhs, d *a, int *lda, d *b, int *ldb, d *x, int *ldx, d *ferr, d *berr, d *work, int *iwork, int *info) nogil
cdef dtrrfs_t *dtrrfs_f = &_fortran_dtrrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtrsen "F_FUNC(dtrsen,DTRSEN)"(char *job, char *compq, bint *select, int *n, d *t, int *ldt, d *q, int *ldq, d *wr, d *wi, int *m, d *s, d *sep, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef dtrsen_t *dtrsen_f = &_fortran_dtrsen

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtrsna "F_FUNC(dtrsna,DTRSNA)"(char *job, char *howmny, bint *select, int *n, d *t, int *ldt, d *vl, int *ldvl, d *vr, int *ldvr, d *s, d *sep, int *mm, int *m, d *work, int *ldwork, int *iwork, int *info) nogil
cdef dtrsna_t *dtrsna_f = &_fortran_dtrsna

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtrsyl "F_FUNC(dtrsyl,DTRSYL)"(char *trana, char *tranb, int *isgn, int *m, int *n, d *a, int *lda, d *b, int *ldb, d *c, int *ldc, d *scale, int *info) nogil
cdef dtrsyl_t *dtrsyl_f = &_fortran_dtrsyl

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtrti2 "F_FUNC(dtrti2,DTRTI2)"(char *uplo, char *diag, int *n, d *a, int *lda, int *info) nogil
cdef dtrti2_t *dtrti2_f = &_fortran_dtrti2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtrtri "F_FUNC(dtrtri,DTRTRI)"(char *uplo, char *diag, int *n, d *a, int *lda, int *info) nogil
cdef dtrtri_t *dtrtri_f = &_fortran_dtrtri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtrtrs "F_FUNC(dtrtrs,DTRTRS)"(char *uplo, char *trans, char *diag, int *n, int *nrhs, d *a, int *lda, d *b, int *ldb, int *info) nogil
cdef dtrtrs_t *dtrtrs_f = &_fortran_dtrtrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtzrqf "F_FUNC(dtzrqf,DTZRQF)"(int *m, int *n, d *a, int *lda, d *tau, int *info) nogil
cdef dtzrqf_t *dtzrqf_f = &_fortran_dtzrqf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_dtzrzf "F_FUNC(dtzrzf,DTZRZF)"(int *m, int *n, d *a, int *lda, d *tau, d *work, int *lwork, int *info) nogil
cdef dtzrzf_t *dtzrzf_f = &_fortran_dtzrzf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sbdsdc "F_FUNC(sbdsdc,SBDSDC)"(char *uplo, char *compq, int *n, s *d, s *e, s *u, int *ldu, s *vt, int *ldvt, s *q, int *iq, s *work, int *iwork, int *info) nogil
cdef sbdsdc_t *sbdsdc_f = &_fortran_sbdsdc

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sbdsqr "F_FUNC(sbdsqr,SBDSQR)"(char *uplo, int *n, int *ncvt, int *nru, int *ncc, s *d, s *e, s *vt, int *ldvt, s *u, int *ldu, s *c, int *ldc, s *work, int *info) nogil
cdef sbdsqr_t *sbdsqr_f = &_fortran_sbdsqr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sdisna "F_FUNC(sdisna,SDISNA)"(char *job, int *m, int *n, s *d, s *sep, int *info) nogil
cdef sdisna_t *sdisna_f = &_fortran_sdisna

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgbbrd "F_FUNC(sgbbrd,SGBBRD)"(char *vect, int *m, int *n, int *ncc, int *kl, int *ku, s *ab, int *ldab, s *d, s *e, s *q, int *ldq, s *pt, int *ldpt, s *c, int *ldc, s *work, int *info) nogil
cdef sgbbrd_t *sgbbrd_f = &_fortran_sgbbrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgbcon "F_FUNC(sgbcon,SGBCON)"(char *norm, int *n, int *kl, int *ku, s *ab, int *ldab, int *ipiv, s *anorm, s *rcond, s *work, int *iwork, int *info) nogil
cdef sgbcon_t *sgbcon_f = &_fortran_sgbcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgbequ "F_FUNC(sgbequ,SGBEQU)"(int *m, int *n, int *kl, int *ku, s *ab, int *ldab, s *r, s *c, s *rowcnd, s *colcnd, s *amax, int *info) nogil
cdef sgbequ_t *sgbequ_f = &_fortran_sgbequ

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgbrfs "F_FUNC(sgbrfs,SGBRFS)"(char *trans, int *n, int *kl, int *ku, int *nrhs, s *ab, int *ldab, s *afb, int *ldafb, int *ipiv, s *b, int *ldb, s *x, int *ldx, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef sgbrfs_t *sgbrfs_f = &_fortran_sgbrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgbsv "F_FUNC(sgbsv,SGBSV)"(int *n, int *kl, int *ku, int *nrhs, s *ab, int *ldab, int *ipiv, s *b, int *ldb, int *info) nogil
cdef sgbsv_t *sgbsv_f = &_fortran_sgbsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgbsvx "F_FUNC(sgbsvx,SGBSVX)"(char *fact, char *trans, int *n, int *kl, int *ku, int *nrhs, s *ab, int *ldab, s *afb, int *ldafb, int *ipiv, char *equed, s *r, s *c, s *b, int *ldb, s *x, int *ldx, s *rcond, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef sgbsvx_t *sgbsvx_f = &_fortran_sgbsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgbtf2 "F_FUNC(sgbtf2,SGBTF2)"(int *m, int *n, int *kl, int *ku, s *ab, int *ldab, int *ipiv, int *info) nogil
cdef sgbtf2_t *sgbtf2_f = &_fortran_sgbtf2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgbtrf "F_FUNC(sgbtrf,SGBTRF)"(int *m, int *n, int *kl, int *ku, s *ab, int *ldab, int *ipiv, int *info) nogil
cdef sgbtrf_t *sgbtrf_f = &_fortran_sgbtrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgbtrs "F_FUNC(sgbtrs,SGBTRS)"(char *trans, int *n, int *kl, int *ku, int *nrhs, s *ab, int *ldab, int *ipiv, s *b, int *ldb, int *info) nogil
cdef sgbtrs_t *sgbtrs_f = &_fortran_sgbtrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgebak "F_FUNC(sgebak,SGEBAK)"(char *job, char *side, int *n, int *ilo, int *ihi, s *scale, int *m, s *v, int *ldv, int *info) nogil
cdef sgebak_t *sgebak_f = &_fortran_sgebak

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgebal "F_FUNC(sgebal,SGEBAL)"(char *job, int *n, s *a, int *lda, int *ilo, int *ihi, s *scale, int *info) nogil
cdef sgebal_t *sgebal_f = &_fortran_sgebal

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgebd2 "F_FUNC(sgebd2,SGEBD2)"(int *m, int *n, s *a, int *lda, s *d, s *e, s *tauq, s *taup, s *work, int *info) nogil
cdef sgebd2_t *sgebd2_f = &_fortran_sgebd2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgebrd "F_FUNC(sgebrd,SGEBRD)"(int *m, int *n, s *a, int *lda, s *d, s *e, s *tauq, s *taup, s *work, int *lwork, int *info) nogil
cdef sgebrd_t *sgebrd_f = &_fortran_sgebrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgecon "F_FUNC(sgecon,SGECON)"(char *norm, int *n, s *a, int *lda, s *anorm, s *rcond, s *work, int *iwork, int *info) nogil
cdef sgecon_t *sgecon_f = &_fortran_sgecon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgeequ "F_FUNC(sgeequ,SGEEQU)"(int *m, int *n, s *a, int *lda, s *r, s *c, s *rowcnd, s *colcnd, s *amax, int *info) nogil
cdef sgeequ_t *sgeequ_f = &_fortran_sgeequ

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgees "F_FUNC(sgees,SGEES)"(char *jobvs, char *sort, sselect2 *select, int *n, s *a, int *lda, int *sdim, s *wr, s *wi, s *vs, int *ldvs, s *work, int *lwork, bint *bwork, int *info) nogil
cdef sgees_t *sgees_f = &_fortran_sgees

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgeesx "F_FUNC(sgeesx,SGEESX)"(char *jobvs, char *sort, sselect2 *select, char *sense, int *n, s *a, int *lda, int *sdim, s *wr, s *wi, s *vs, int *ldvs, s *rconde, s *rcondv, s *work, int *lwork, int *iwork, int *liwork, bint *bwork, int *info) nogil
cdef sgeesx_t *sgeesx_f = &_fortran_sgeesx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgeev "F_FUNC(sgeev,SGEEV)"(char *jobvl, char *jobvr, int *n, s *a, int *lda, s *wr, s *wi, s *vl, int *ldvl, s *vr, int *ldvr, s *work, int *lwork, int *info) nogil
cdef sgeev_t *sgeev_f = &_fortran_sgeev

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgeevx "F_FUNC(sgeevx,SGEEVX)"(char *balanc, char *jobvl, char *jobvr, char *sense, int *n, s *a, int *lda, s *wr, s *wi, s *vl, int *ldvl, s *vr, int *ldvr, int *ilo, int *ihi, s *scale, s *abnrm, s *rconde, s *rcondv, s *work, int *lwork, int *iwork, int *info) nogil
cdef sgeevx_t *sgeevx_f = &_fortran_sgeevx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgegs "F_FUNC(sgegs,SGEGS)"(char *jobvsl, char *jobvsr, int *n, s *a, int *lda, s *b, int *ldb, s *alphar, s *alphai, s *beta, s *vsl, int *ldvsl, s *vsr, int *ldvsr, s *work, int *lwork, int *info) nogil
cdef sgegs_t *sgegs_f = &_fortran_sgegs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgegv "F_FUNC(sgegv,SGEGV)"(char *jobvl, char *jobvr, int *n, s *a, int *lda, s *b, int *ldb, s *alphar, s *alphai, s *beta, s *vl, int *ldvl, s *vr, int *ldvr, s *work, int *lwork, int *info) nogil
cdef sgegv_t *sgegv_f = &_fortran_sgegv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgehd2 "F_FUNC(sgehd2,SGEHD2)"(int *n, int *ilo, int *ihi, s *a, int *lda, s *tau, s *work, int *info) nogil
cdef sgehd2_t *sgehd2_f = &_fortran_sgehd2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgehrd "F_FUNC(sgehrd,SGEHRD)"(int *n, int *ilo, int *ihi, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef sgehrd_t *sgehrd_f = &_fortran_sgehrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgelq2 "F_FUNC(sgelq2,SGELQ2)"(int *m, int *n, s *a, int *lda, s *tau, s *work, int *info) nogil
cdef sgelq2_t *sgelq2_f = &_fortran_sgelq2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgelqf "F_FUNC(sgelqf,SGELQF)"(int *m, int *n, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef sgelqf_t *sgelqf_f = &_fortran_sgelqf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgels "F_FUNC(sgels,SGELS)"(char *trans, int *m, int *n, int *nrhs, s *a, int *lda, s *b, int *ldb, s *work, int *lwork, int *info) nogil
cdef sgels_t *sgels_f = &_fortran_sgels

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgelsd "F_FUNC(sgelsd,SGELSD)"(int *m, int *n, int *nrhs, s *a, int *lda, s *b, int *ldb, s *s, s *rcond, int *rank, s *work, int *lwork, int *iwork, int *info) nogil
cdef sgelsd_t *sgelsd_f = &_fortran_sgelsd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgelss "F_FUNC(sgelss,SGELSS)"(int *m, int *n, int *nrhs, s *a, int *lda, s *b, int *ldb, s *s, s *rcond, int *rank, s *work, int *lwork, int *info) nogil
cdef sgelss_t *sgelss_f = &_fortran_sgelss

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgelsx "F_FUNC(sgelsx,SGELSX)"(int *m, int *n, int *nrhs, s *a, int *lda, s *b, int *ldb, int *jpvt, s *rcond, int *rank, s *work, int *info) nogil
cdef sgelsx_t *sgelsx_f = &_fortran_sgelsx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgelsy "F_FUNC(sgelsy,SGELSY)"(int *m, int *n, int *nrhs, s *a, int *lda, s *b, int *ldb, int *jpvt, s *rcond, int *rank, s *work, int *lwork, int *info) nogil
cdef sgelsy_t *sgelsy_f = &_fortran_sgelsy

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgeql2 "F_FUNC(sgeql2,SGEQL2)"(int *m, int *n, s *a, int *lda, s *tau, s *work, int *info) nogil
cdef sgeql2_t *sgeql2_f = &_fortran_sgeql2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgeqlf "F_FUNC(sgeqlf,SGEQLF)"(int *m, int *n, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef sgeqlf_t *sgeqlf_f = &_fortran_sgeqlf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgeqp3 "F_FUNC(sgeqp3,SGEQP3)"(int *m, int *n, s *a, int *lda, int *jpvt, s *tau, s *work, int *lwork, int *info) nogil
cdef sgeqp3_t *sgeqp3_f = &_fortran_sgeqp3

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgeqpf "F_FUNC(sgeqpf,SGEQPF)"(int *m, int *n, s *a, int *lda, int *jpvt, s *tau, s *work, int *info) nogil
cdef sgeqpf_t *sgeqpf_f = &_fortran_sgeqpf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgeqr2 "F_FUNC(sgeqr2,SGEQR2)"(int *m, int *n, s *a, int *lda, s *tau, s *work, int *info) nogil
cdef sgeqr2_t *sgeqr2_f = &_fortran_sgeqr2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgeqrf "F_FUNC(sgeqrf,SGEQRF)"(int *m, int *n, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef sgeqrf_t *sgeqrf_f = &_fortran_sgeqrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgerfs "F_FUNC(sgerfs,SGERFS)"(char *trans, int *n, int *nrhs, s *a, int *lda, s *af, int *ldaf, int *ipiv, s *b, int *ldb, s *x, int *ldx, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef sgerfs_t *sgerfs_f = &_fortran_sgerfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgerq2 "F_FUNC(sgerq2,SGERQ2)"(int *m, int *n, s *a, int *lda, s *tau, s *work, int *info) nogil
cdef sgerq2_t *sgerq2_f = &_fortran_sgerq2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgerqf "F_FUNC(sgerqf,SGERQF)"(int *m, int *n, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef sgerqf_t *sgerqf_f = &_fortran_sgerqf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgesc2 "F_FUNC(sgesc2,SGESC2)"(int *n, s *a, int *lda, s *rhs, int *ipiv, int *jpiv, s *scale) nogil
cdef sgesc2_t *sgesc2_f = &_fortran_sgesc2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgesdd "F_FUNC(sgesdd,SGESDD)"(char *jobz, int *m, int *n, s *a, int *lda, s *s, s *u, int *ldu, s *vt, int *ldvt, s *work, int *lwork, int *iwork, int *info) nogil
cdef sgesdd_t *sgesdd_f = &_fortran_sgesdd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgesv "F_FUNC(sgesv,SGESV)"(int *n, int *nrhs, s *a, int *lda, int *ipiv, s *b, int *ldb, int *info) nogil
cdef sgesv_t *sgesv_f = &_fortran_sgesv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgesvd "F_FUNC(sgesvd,SGESVD)"(char *jobu, char *jobvt, int *m, int *n, s *a, int *lda, s *s, s *u, int *ldu, s *vt, int *ldvt, s *work, int *lwork, int *info) nogil
cdef sgesvd_t *sgesvd_f = &_fortran_sgesvd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgesvx "F_FUNC(sgesvx,SGESVX)"(char *fact, char *trans, int *n, int *nrhs, s *a, int *lda, s *af, int *ldaf, int *ipiv, char *equed, s *r, s *c, s *b, int *ldb, s *x, int *ldx, s *rcond, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef sgesvx_t *sgesvx_f = &_fortran_sgesvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgetc2 "F_FUNC(sgetc2,SGETC2)"(int *n, s *a, int *lda, int *ipiv, int *jpiv, int *info) nogil
cdef sgetc2_t *sgetc2_f = &_fortran_sgetc2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgetf2 "F_FUNC(sgetf2,SGETF2)"(int *m, int *n, s *a, int *lda, int *ipiv, int *info) nogil
cdef sgetf2_t *sgetf2_f = &_fortran_sgetf2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgetrf "F_FUNC(sgetrf,SGETRF)"(int *m, int *n, s *a, int *lda, int *ipiv, int *info) nogil
cdef sgetrf_t *sgetrf_f = &_fortran_sgetrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgetri "F_FUNC(sgetri,SGETRI)"(int *n, s *a, int *lda, int *ipiv, s *work, int *lwork, int *info) nogil
cdef sgetri_t *sgetri_f = &_fortran_sgetri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgetrs "F_FUNC(sgetrs,SGETRS)"(char *trans, int *n, int *nrhs, s *a, int *lda, int *ipiv, s *b, int *ldb, int *info) nogil
cdef sgetrs_t *sgetrs_f = &_fortran_sgetrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sggbak "F_FUNC(sggbak,SGGBAK)"(char *job, char *side, int *n, int *ilo, int *ihi, s *lscale, s *rscale, int *m, s *v, int *ldv, int *info) nogil
cdef sggbak_t *sggbak_f = &_fortran_sggbak

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sggbal "F_FUNC(sggbal,SGGBAL)"(char *job, int *n, s *a, int *lda, s *b, int *ldb, int *ilo, int *ihi, s *lscale, s *rscale, s *work, int *info) nogil
cdef sggbal_t *sggbal_f = &_fortran_sggbal

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgges "F_FUNC(sgges,SGGES)"(char *jobvsl, char *jobvsr, char *sort, sselect3 *selctg, int *n, s *a, int *lda, s *b, int *ldb, int *sdim, s *alphar, s *alphai, s *beta, s *vsl, int *ldvsl, s *vsr, int *ldvsr, s *work, int *lwork, bint *bwork, int *info) nogil
cdef sgges_t *sgges_f = &_fortran_sgges

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sggesx "F_FUNC(sggesx,SGGESX)"(char *jobvsl, char *jobvsr, char *sort, sselect3 *selctg, char *sense, int *n, s *a, int *lda, s *b, int *ldb, int *sdim, s *alphar, s *alphai, s *beta, s *vsl, int *ldvsl, s *vsr, int *ldvsr, s *rconde, s *rcondv, s *work, int *lwork, int *iwork, int *liwork, bint *bwork, int *info) nogil
cdef sggesx_t *sggesx_f = &_fortran_sggesx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sggev "F_FUNC(sggev,SGGEV)"(char *jobvl, char *jobvr, int *n, s *a, int *lda, s *b, int *ldb, s *alphar, s *alphai, s *beta, s *vl, int *ldvl, s *vr, int *ldvr, s *work, int *lwork, int *info) nogil
cdef sggev_t *sggev_f = &_fortran_sggev

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sggevx "F_FUNC(sggevx,SGGEVX)"(char *balanc, char *jobvl, char *jobvr, char *sense, int *n, s *a, int *lda, s *b, int *ldb, s *alphar, s *alphai, s *beta, s *vl, int *ldvl, s *vr, int *ldvr, int *ilo, int *ihi, s *lscale, s *rscale, s *abnrm, s *bbnrm, s *rconde, s *rcondv, s *work, int *lwork, int *iwork, bint *bwork, int *info) nogil
cdef sggevx_t *sggevx_f = &_fortran_sggevx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sggglm "F_FUNC(sggglm,SGGGLM)"(int *n, int *m, int *p, s *a, int *lda, s *b, int *ldb, s *d, s *x, s *y, s *work, int *lwork, int *info) nogil
cdef sggglm_t *sggglm_f = &_fortran_sggglm

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgghrd "F_FUNC(sgghrd,SGGHRD)"(char *compq, char *compz, int *n, int *ilo, int *ihi, s *a, int *lda, s *b, int *ldb, s *q, int *ldq, s *z, int *ldz, int *info) nogil
cdef sgghrd_t *sgghrd_f = &_fortran_sgghrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgglse "F_FUNC(sgglse,SGGLSE)"(int *m, int *n, int *p, s *a, int *lda, s *b, int *ldb, s *c, s *d, s *x, s *work, int *lwork, int *info) nogil
cdef sgglse_t *sgglse_f = &_fortran_sgglse

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sggqrf "F_FUNC(sggqrf,SGGQRF)"(int *n, int *m, int *p, s *a, int *lda, s *taua, s *b, int *ldb, s *taub, s *work, int *lwork, int *info) nogil
cdef sggqrf_t *sggqrf_f = &_fortran_sggqrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sggrqf "F_FUNC(sggrqf,SGGRQF)"(int *m, int *p, int *n, s *a, int *lda, s *taua, s *b, int *ldb, s *taub, s *work, int *lwork, int *info) nogil
cdef sggrqf_t *sggrqf_f = &_fortran_sggrqf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sggsvd "F_FUNC(sggsvd,SGGSVD)"(char *jobu, char *jobv, char *jobq, int *m, int *n, int *p, int *k, int *l, s *a, int *lda, s *b, int *ldb, s *alpha, s *beta, s *u, int *ldu, s *v, int *ldv, s *q, int *ldq, s *work, int *iwork, int *info) nogil
cdef sggsvd_t *sggsvd_f = &_fortran_sggsvd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sggsvp "F_FUNC(sggsvp,SGGSVP)"(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, s *a, int *lda, s *b, int *ldb, s *tola, s *tolb, int *k, int *l, s *u, int *ldu, s *v, int *ldv, s *q, int *ldq, int *iwork, s *tau, s *work, int *info) nogil
cdef sggsvp_t *sggsvp_f = &_fortran_sggsvp

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgtcon "F_FUNC(sgtcon,SGTCON)"(char *norm, int *n, s *dl, s *d, s *du, s *du2, int *ipiv, s *anorm, s *rcond, s *work, int *iwork, int *info) nogil
cdef sgtcon_t *sgtcon_f = &_fortran_sgtcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgtrfs "F_FUNC(sgtrfs,SGTRFS)"(char *trans, int *n, int *nrhs, s *dl, s *d, s *du, s *dlf, s *df, s *duf, s *du2, int *ipiv, s *b, int *ldb, s *x, int *ldx, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef sgtrfs_t *sgtrfs_f = &_fortran_sgtrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgtsv "F_FUNC(sgtsv,SGTSV)"(int *n, int *nrhs, s *dl, s *d, s *du, s *b, int *ldb, int *info) nogil
cdef sgtsv_t *sgtsv_f = &_fortran_sgtsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgtsvx "F_FUNC(sgtsvx,SGTSVX)"(char *fact, char *trans, int *n, int *nrhs, s *dl, s *d, s *du, s *dlf, s *df, s *duf, s *du2, int *ipiv, s *b, int *ldb, s *x, int *ldx, s *rcond, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef sgtsvx_t *sgtsvx_f = &_fortran_sgtsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgttrf "F_FUNC(sgttrf,SGTTRF)"(int *n, s *dl, s *d, s *du, s *du2, int *ipiv, int *info) nogil
cdef sgttrf_t *sgttrf_f = &_fortran_sgttrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgttrs "F_FUNC(sgttrs,SGTTRS)"(char *trans, int *n, int *nrhs, s *dl, s *d, s *du, s *du2, int *ipiv, s *b, int *ldb, int *info) nogil
cdef sgttrs_t *sgttrs_f = &_fortran_sgttrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sgtts2 "F_FUNC(sgtts2,SGTTS2)"(int *itrans, int *n, int *nrhs, s *dl, s *d, s *du, s *du2, int *ipiv, s *b, int *ldb) nogil
cdef sgtts2_t *sgtts2_f = &_fortran_sgtts2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_shgeqz "F_FUNC(shgeqz,SHGEQZ)"(char *job, char *compq, char *compz, int *n, int *ilo, int *ihi, s *h, int *ldh, s *t, int *ldt, s *alphar, s *alphai, s *beta, s *q, int *ldq, s *z, int *ldz, s *work, int *lwork, int *info) nogil
cdef shgeqz_t *shgeqz_f = &_fortran_shgeqz

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_shsein "F_FUNC(shsein,SHSEIN)"(char *side, char *eigsrc, char *initv, bint *select, int *n, s *h, int *ldh, s *wr, s *wi, s *vl, int *ldvl, s *vr, int *ldvr, int *mm, int *m, s *work, int *ifaill, int *ifailr, int *info) nogil
cdef shsein_t *shsein_f = &_fortran_shsein

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_shseqr "F_FUNC(shseqr,SHSEQR)"(char *job, char *compz, int *n, int *ilo, int *ihi, s *h, int *ldh, s *wr, s *wi, s *z, int *ldz, s *work, int *lwork, int *info) nogil
cdef shseqr_t *shseqr_f = &_fortran_shseqr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_slacn2 "F_FUNC(slacn2,SLACN2)"(int *n, s *v, s *x, int *isgn, s *est, int *kase, int *isave) nogil
cdef slacn2_t *slacn2_f = &_fortran_slacn2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_slacon "F_FUNC(slacon,SLACON)"(int *n, s *v, s *x, int *isgn, s *est, int *kase) nogil
cdef slacon_t *slacon_f = &_fortran_slacon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_slanv2 "F_FUNC(slanv2,SLANV2)"(s *a, s *b, s *c, s *d, s *rt1r, s *rt1i, s *rt2r, s *rt2i, s *cs, s *sn) nogil
cdef slanv2_t *slanv2_f = &_fortran_slanv2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_slarf "F_FUNC(slarf,SLARF)"(char *side, int *m, int *n, s *v, int *incv, s *tau, s *c, int *ldc, s *work) nogil
cdef slarf_t *slarf_f = &_fortran_slarf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_slarz "F_FUNC(slarz,SLARZ)"(char *side, int *m, int *n, int *l, s *v, int *incv, s *tau, s *c, int *ldc, s *work) nogil
cdef slarz_t *slarz_f = &_fortran_slarz

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_slaswp "F_FUNC(slaswp,SLASWP)"(int *n, s *a, int *lda, int *k1, int *k2, int *ipiv, int *incx) nogil
cdef slaswp_t *slaswp_f = &_fortran_slaswp

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_slauum "F_FUNC(slauum,SLAUUM)"(char *uplo, int *n, s *a, int *lda, int *info) nogil
cdef slauum_t *slauum_f = &_fortran_slauum

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sopgtr "F_FUNC(sopgtr,SOPGTR)"(char *uplo, int *n, s *ap, s *tau, s *q, int *ldq, s *work, int *info) nogil
cdef sopgtr_t *sopgtr_f = &_fortran_sopgtr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sopmtr "F_FUNC(sopmtr,SOPMTR)"(char *side, char *uplo, char *trans, int *m, int *n, s *ap, s *tau, s *c, int *ldc, s *work, int *info) nogil
cdef sopmtr_t *sopmtr_f = &_fortran_sopmtr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sorg2l "F_FUNC(sorg2l,SORG2L)"(int *m, int *n, int *k, s *a, int *lda, s *tau, s *work, int *info) nogil
cdef sorg2l_t *sorg2l_f = &_fortran_sorg2l

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sorg2r "F_FUNC(sorg2r,SORG2R)"(int *m, int *n, int *k, s *a, int *lda, s *tau, s *work, int *info) nogil
cdef sorg2r_t *sorg2r_f = &_fortran_sorg2r

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sorgbr "F_FUNC(sorgbr,SORGBR)"(char *vect, int *m, int *n, int *k, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef sorgbr_t *sorgbr_f = &_fortran_sorgbr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sorghr "F_FUNC(sorghr,SORGHR)"(int *n, int *ilo, int *ihi, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef sorghr_t *sorghr_f = &_fortran_sorghr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sorgl2 "F_FUNC(sorgl2,SORGL2)"(int *m, int *n, int *k, s *a, int *lda, s *tau, s *work, int *info) nogil
cdef sorgl2_t *sorgl2_f = &_fortran_sorgl2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sorglq "F_FUNC(sorglq,SORGLQ)"(int *m, int *n, int *k, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef sorglq_t *sorglq_f = &_fortran_sorglq

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sorgql "F_FUNC(sorgql,SORGQL)"(int *m, int *n, int *k, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef sorgql_t *sorgql_f = &_fortran_sorgql

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sorgqr "F_FUNC(sorgqr,SORGQR)"(int *m, int *n, int *k, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef sorgqr_t *sorgqr_f = &_fortran_sorgqr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sorgr2 "F_FUNC(sorgr2,SORGR2)"(int *m, int *n, int *k, s *a, int *lda, s *tau, s *work, int *info) nogil
cdef sorgr2_t *sorgr2_f = &_fortran_sorgr2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sorgrq "F_FUNC(sorgrq,SORGRQ)"(int *m, int *n, int *k, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef sorgrq_t *sorgrq_f = &_fortran_sorgrq

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sorgtr "F_FUNC(sorgtr,SORGTR)"(char *uplo, int *n, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef sorgtr_t *sorgtr_f = &_fortran_sorgtr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sorm2l "F_FUNC(sorm2l,SORM2L)"(char *side, char *trans, int *m, int *n, int *k, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *info) nogil
cdef sorm2l_t *sorm2l_f = &_fortran_sorm2l

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sorm2r "F_FUNC(sorm2r,SORM2R)"(char *side, char *trans, int *m, int *n, int *k, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *info) nogil
cdef sorm2r_t *sorm2r_f = &_fortran_sorm2r

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sormbr "F_FUNC(sormbr,SORMBR)"(char *vect, char *side, char *trans, int *m, int *n, int *k, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *lwork, int *info) nogil
cdef sormbr_t *sormbr_f = &_fortran_sormbr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sormhr "F_FUNC(sormhr,SORMHR)"(char *side, char *trans, int *m, int *n, int *ilo, int *ihi, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *lwork, int *info) nogil
cdef sormhr_t *sormhr_f = &_fortran_sormhr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sorml2 "F_FUNC(sorml2,SORML2)"(char *side, char *trans, int *m, int *n, int *k, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *info) nogil
cdef sorml2_t *sorml2_f = &_fortran_sorml2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sormlq "F_FUNC(sormlq,SORMLQ)"(char *side, char *trans, int *m, int *n, int *k, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *lwork, int *info) nogil
cdef sormlq_t *sormlq_f = &_fortran_sormlq

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sormql "F_FUNC(sormql,SORMQL)"(char *side, char *trans, int *m, int *n, int *k, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *lwork, int *info) nogil
cdef sormql_t *sormql_f = &_fortran_sormql

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sormqr "F_FUNC(sormqr,SORMQR)"(char *side, char *trans, int *m, int *n, int *k, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *lwork, int *info) nogil
cdef sormqr_t *sormqr_f = &_fortran_sormqr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sormr2 "F_FUNC(sormr2,SORMR2)"(char *side, char *trans, int *m, int *n, int *k, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *info) nogil
cdef sormr2_t *sormr2_f = &_fortran_sormr2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sormr3 "F_FUNC(sormr3,SORMR3)"(char *side, char *trans, int *m, int *n, int *k, int *l, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *info) nogil
cdef sormr3_t *sormr3_f = &_fortran_sormr3

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sormrq "F_FUNC(sormrq,SORMRQ)"(char *side, char *trans, int *m, int *n, int *k, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *lwork, int *info) nogil
cdef sormrq_t *sormrq_f = &_fortran_sormrq

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sormrz "F_FUNC(sormrz,SORMRZ)"(char *side, char *trans, int *m, int *n, int *k, int *l, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *lwork, int *info) nogil
cdef sormrz_t *sormrz_f = &_fortran_sormrz

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sormtr "F_FUNC(sormtr,SORMTR)"(char *side, char *uplo, char *trans, int *m, int *n, s *a, int *lda, s *tau, s *c, int *ldc, s *work, int *lwork, int *info) nogil
cdef sormtr_t *sormtr_f = &_fortran_sormtr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_spbcon "F_FUNC(spbcon,SPBCON)"(char *uplo, int *n, int *kd, s *ab, int *ldab, s *anorm, s *rcond, s *work, int *iwork, int *info) nogil
cdef spbcon_t *spbcon_f = &_fortran_spbcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_spbequ "F_FUNC(spbequ,SPBEQU)"(char *uplo, int *n, int *kd, s *ab, int *ldab, s *s, s *scond, s *amax, int *info) nogil
cdef spbequ_t *spbequ_f = &_fortran_spbequ

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_spbrfs "F_FUNC(spbrfs,SPBRFS)"(char *uplo, int *n, int *kd, int *nrhs, s *ab, int *ldab, s *afb, int *ldafb, s *b, int *ldb, s *x, int *ldx, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef spbrfs_t *spbrfs_f = &_fortran_spbrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_spbstf "F_FUNC(spbstf,SPBSTF)"(char *uplo, int *n, int *kd, s *ab, int *ldab, int *info) nogil
cdef spbstf_t *spbstf_f = &_fortran_spbstf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_spbsv "F_FUNC(spbsv,SPBSV)"(char *uplo, int *n, int *kd, int *nrhs, s *ab, int *ldab, s *b, int *ldb, int *info) nogil
cdef spbsv_t *spbsv_f = &_fortran_spbsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_spbsvx "F_FUNC(spbsvx,SPBSVX)"(char *fact, char *uplo, int *n, int *kd, int *nrhs, s *ab, int *ldab, s *afb, int *ldafb, char *equed, s *s, s *b, int *ldb, s *x, int *ldx, s *rcond, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef spbsvx_t *spbsvx_f = &_fortran_spbsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_spbtf2 "F_FUNC(spbtf2,SPBTF2)"(char *uplo, int *n, int *kd, s *ab, int *ldab, int *info) nogil
cdef spbtf2_t *spbtf2_f = &_fortran_spbtf2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_spbtrf "F_FUNC(spbtrf,SPBTRF)"(char *uplo, int *n, int *kd, s *ab, int *ldab, int *info) nogil
cdef spbtrf_t *spbtrf_f = &_fortran_spbtrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_spbtrs "F_FUNC(spbtrs,SPBTRS)"(char *uplo, int *n, int *kd, int *nrhs, s *ab, int *ldab, s *b, int *ldb, int *info) nogil
cdef spbtrs_t *spbtrs_f = &_fortran_spbtrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_spocon "F_FUNC(spocon,SPOCON)"(char *uplo, int *n, s *a, int *lda, s *anorm, s *rcond, s *work, int *iwork, int *info) nogil
cdef spocon_t *spocon_f = &_fortran_spocon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_spoequ "F_FUNC(spoequ,SPOEQU)"(int *n, s *a, int *lda, s *s, s *scond, s *amax, int *info) nogil
cdef spoequ_t *spoequ_f = &_fortran_spoequ

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sporfs "F_FUNC(sporfs,SPORFS)"(char *uplo, int *n, int *nrhs, s *a, int *lda, s *af, int *ldaf, s *b, int *ldb, s *x, int *ldx, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef sporfs_t *sporfs_f = &_fortran_sporfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sposv "F_FUNC(sposv,SPOSV)"(char *uplo, int *n, int *nrhs, s *a, int *lda, s *b, int *ldb, int *info) nogil
cdef sposv_t *sposv_f = &_fortran_sposv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sposvx "F_FUNC(sposvx,SPOSVX)"(char *fact, char *uplo, int *n, int *nrhs, s *a, int *lda, s *af, int *ldaf, char *equed, s *s, s *b, int *ldb, s *x, int *ldx, s *rcond, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef sposvx_t *sposvx_f = &_fortran_sposvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_spotf2 "F_FUNC(spotf2,SPOTF2)"(char *uplo, int *n, s *a, int *lda, int *info) nogil
cdef spotf2_t *spotf2_f = &_fortran_spotf2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_spotrf "F_FUNC(spotrf,SPOTRF)"(char *uplo, int *n, s *a, int *lda, int *info) nogil
cdef spotrf_t *spotrf_f = &_fortran_spotrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_spotri "F_FUNC(spotri,SPOTRI)"(char *uplo, int *n, s *a, int *lda, int *info) nogil
cdef spotri_t *spotri_f = &_fortran_spotri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_spotrs "F_FUNC(spotrs,SPOTRS)"(char *uplo, int *n, int *nrhs, s *a, int *lda, s *b, int *ldb, int *info) nogil
cdef spotrs_t *spotrs_f = &_fortran_spotrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sppcon "F_FUNC(sppcon,SPPCON)"(char *uplo, int *n, s *ap, s *anorm, s *rcond, s *work, int *iwork, int *info) nogil
cdef sppcon_t *sppcon_f = &_fortran_sppcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sppequ "F_FUNC(sppequ,SPPEQU)"(char *uplo, int *n, s *ap, s *s, s *scond, s *amax, int *info) nogil
cdef sppequ_t *sppequ_f = &_fortran_sppequ

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_spprfs "F_FUNC(spprfs,SPPRFS)"(char *uplo, int *n, int *nrhs, s *ap, s *afp, s *b, int *ldb, s *x, int *ldx, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef spprfs_t *spprfs_f = &_fortran_spprfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sppsv "F_FUNC(sppsv,SPPSV)"(char *uplo, int *n, int *nrhs, s *ap, s *b, int *ldb, int *info) nogil
cdef sppsv_t *sppsv_f = &_fortran_sppsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sppsvx "F_FUNC(sppsvx,SPPSVX)"(char *fact, char *uplo, int *n, int *nrhs, s *ap, s *afp, char *equed, s *s, s *b, int *ldb, s *x, int *ldx, s *rcond, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef sppsvx_t *sppsvx_f = &_fortran_sppsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_spptrf "F_FUNC(spptrf,SPPTRF)"(char *uplo, int *n, s *ap, int *info) nogil
cdef spptrf_t *spptrf_f = &_fortran_spptrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_spptri "F_FUNC(spptri,SPPTRI)"(char *uplo, int *n, s *ap, int *info) nogil
cdef spptri_t *spptri_f = &_fortran_spptri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_spptrs "F_FUNC(spptrs,SPPTRS)"(char *uplo, int *n, int *nrhs, s *ap, s *b, int *ldb, int *info) nogil
cdef spptrs_t *spptrs_f = &_fortran_spptrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sptcon "F_FUNC(sptcon,SPTCON)"(int *n, s *d, s *e, s *anorm, s *rcond, s *work, int *info) nogil
cdef sptcon_t *sptcon_f = &_fortran_sptcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_spteqr "F_FUNC(spteqr,SPTEQR)"(char *compz, int *n, s *d, s *e, s *z, int *ldz, s *work, int *info) nogil
cdef spteqr_t *spteqr_f = &_fortran_spteqr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sptrfs "F_FUNC(sptrfs,SPTRFS)"(int *n, int *nrhs, s *d, s *e, s *df, s *ef, s *b, int *ldb, s *x, int *ldx, s *ferr, s *berr, s *work, int *info) nogil
cdef sptrfs_t *sptrfs_f = &_fortran_sptrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sptsv "F_FUNC(sptsv,SPTSV)"(int *n, int *nrhs, s *d, s *e, s *b, int *ldb, int *info) nogil
cdef sptsv_t *sptsv_f = &_fortran_sptsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sptsvx "F_FUNC(sptsvx,SPTSVX)"(char *fact, int *n, int *nrhs, s *d, s *e, s *df, s *ef, s *b, int *ldb, s *x, int *ldx, s *rcond, s *ferr, s *berr, s *work, int *info) nogil
cdef sptsvx_t *sptsvx_f = &_fortran_sptsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_spttrf "F_FUNC(spttrf,SPTTRF)"(int *n, s *d, s *e, int *info) nogil
cdef spttrf_t *spttrf_f = &_fortran_spttrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_spttrs "F_FUNC(spttrs,SPTTRS)"(int *n, int *nrhs, s *d, s *e, s *b, int *ldb, int *info) nogil
cdef spttrs_t *spttrs_f = &_fortran_spttrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sptts2 "F_FUNC(sptts2,SPTTS2)"(int *n, int *nrhs, s *d, s *e, s *b, int *ldb) nogil
cdef sptts2_t *sptts2_f = &_fortran_sptts2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_srscl "F_FUNC(srscl,SRSCL)"(int *n, s *sa, s *sx, int *incx) nogil
cdef srscl_t *srscl_f = &_fortran_srscl

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssbev "F_FUNC(ssbev,SSBEV)"(char *jobz, char *uplo, int *n, int *kd, s *ab, int *ldab, s *w, s *z, int *ldz, s *work, int *info) nogil
cdef ssbev_t *ssbev_f = &_fortran_ssbev

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssbevd "F_FUNC(ssbevd,SSBEVD)"(char *jobz, char *uplo, int *n, int *kd, s *ab, int *ldab, s *w, s *z, int *ldz, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef ssbevd_t *ssbevd_f = &_fortran_ssbevd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssbevx "F_FUNC(ssbevx,SSBEVX)"(char *jobz, char *range, char *uplo, int *n, int *kd, s *ab, int *ldab, s *q, int *ldq, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, s *z, int *ldz, s *work, int *iwork, int *ifail, int *info) nogil
cdef ssbevx_t *ssbevx_f = &_fortran_ssbevx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssbgst "F_FUNC(ssbgst,SSBGST)"(char *vect, char *uplo, int *n, int *ka, int *kb, s *ab, int *ldab, s *bb, int *ldbb, s *x, int *ldx, s *work, int *info) nogil
cdef ssbgst_t *ssbgst_f = &_fortran_ssbgst

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssbgv "F_FUNC(ssbgv,SSBGV)"(char *jobz, char *uplo, int *n, int *ka, int *kb, s *ab, int *ldab, s *bb, int *ldbb, s *w, s *z, int *ldz, s *work, int *info) nogil
cdef ssbgv_t *ssbgv_f = &_fortran_ssbgv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssbgvd "F_FUNC(ssbgvd,SSBGVD)"(char *jobz, char *uplo, int *n, int *ka, int *kb, s *ab, int *ldab, s *bb, int *ldbb, s *w, s *z, int *ldz, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef ssbgvd_t *ssbgvd_f = &_fortran_ssbgvd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssbgvx "F_FUNC(ssbgvx,SSBGVX)"(char *jobz, char *range, char *uplo, int *n, int *ka, int *kb, s *ab, int *ldab, s *bb, int *ldbb, s *q, int *ldq, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, s *z, int *ldz, s *work, int *iwork, int *ifail, int *info) nogil
cdef ssbgvx_t *ssbgvx_f = &_fortran_ssbgvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssbtrd "F_FUNC(ssbtrd,SSBTRD)"(char *vect, char *uplo, int *n, int *kd, s *ab, int *ldab, s *d, s *e, s *q, int *ldq, s *work, int *info) nogil
cdef ssbtrd_t *ssbtrd_f = &_fortran_ssbtrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sspcon "F_FUNC(sspcon,SSPCON)"(char *uplo, int *n, s *ap, int *ipiv, s *anorm, s *rcond, s *work, int *iwork, int *info) nogil
cdef sspcon_t *sspcon_f = &_fortran_sspcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sspev "F_FUNC(sspev,SSPEV)"(char *jobz, char *uplo, int *n, s *ap, s *w, s *z, int *ldz, s *work, int *info) nogil
cdef sspev_t *sspev_f = &_fortran_sspev

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sspevd "F_FUNC(sspevd,SSPEVD)"(char *jobz, char *uplo, int *n, s *ap, s *w, s *z, int *ldz, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef sspevd_t *sspevd_f = &_fortran_sspevd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sspevx "F_FUNC(sspevx,SSPEVX)"(char *jobz, char *range, char *uplo, int *n, s *ap, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, s *z, int *ldz, s *work, int *iwork, int *ifail, int *info) nogil
cdef sspevx_t *sspevx_f = &_fortran_sspevx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sspgst "F_FUNC(sspgst,SSPGST)"(int *itype, char *uplo, int *n, s *ap, s *bp, int *info) nogil
cdef sspgst_t *sspgst_f = &_fortran_sspgst

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sspgv "F_FUNC(sspgv,SSPGV)"(int *itype, char *jobz, char *uplo, int *n, s *ap, s *bp, s *w, s *z, int *ldz, s *work, int *info) nogil
cdef sspgv_t *sspgv_f = &_fortran_sspgv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sspgvd "F_FUNC(sspgvd,SSPGVD)"(int *itype, char *jobz, char *uplo, int *n, s *ap, s *bp, s *w, s *z, int *ldz, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef sspgvd_t *sspgvd_f = &_fortran_sspgvd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sspgvx "F_FUNC(sspgvx,SSPGVX)"(int *itype, char *jobz, char *range, char *uplo, int *n, s *ap, s *bp, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, s *z, int *ldz, s *work, int *iwork, int *ifail, int *info) nogil
cdef sspgvx_t *sspgvx_f = &_fortran_sspgvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssprfs "F_FUNC(ssprfs,SSPRFS)"(char *uplo, int *n, int *nrhs, s *ap, s *afp, int *ipiv, s *b, int *ldb, s *x, int *ldx, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef ssprfs_t *ssprfs_f = &_fortran_ssprfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sspsv "F_FUNC(sspsv,SSPSV)"(char *uplo, int *n, int *nrhs, s *ap, int *ipiv, s *b, int *ldb, int *info) nogil
cdef sspsv_t *sspsv_f = &_fortran_sspsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sspsvx "F_FUNC(sspsvx,SSPSVX)"(char *fact, char *uplo, int *n, int *nrhs, s *ap, s *afp, int *ipiv, s *b, int *ldb, s *x, int *ldx, s *rcond, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef sspsvx_t *sspsvx_f = &_fortran_sspsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssptrd "F_FUNC(ssptrd,SSPTRD)"(char *uplo, int *n, s *ap, s *d, s *e, s *tau, int *info) nogil
cdef ssptrd_t *ssptrd_f = &_fortran_ssptrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssptrf "F_FUNC(ssptrf,SSPTRF)"(char *uplo, int *n, s *ap, int *ipiv, int *info) nogil
cdef ssptrf_t *ssptrf_f = &_fortran_ssptrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssptri "F_FUNC(ssptri,SSPTRI)"(char *uplo, int *n, s *ap, int *ipiv, s *work, int *info) nogil
cdef ssptri_t *ssptri_f = &_fortran_ssptri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssptrs "F_FUNC(ssptrs,SSPTRS)"(char *uplo, int *n, int *nrhs, s *ap, int *ipiv, s *b, int *ldb, int *info) nogil
cdef ssptrs_t *ssptrs_f = &_fortran_ssptrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sstebz "F_FUNC(sstebz,SSTEBZ)"(char *range, char *order, int *n, s *vl, s *vu, int *il, int *iu, s *abstol, s *d, s *e, int *m, int *nsplit, s *w, int *iblock, int *isplit, s *work, int *iwork, int *info) nogil
cdef sstebz_t *sstebz_f = &_fortran_sstebz

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sstedc "F_FUNC(sstedc,SSTEDC)"(char *compz, int *n, s *d, s *e, s *z, int *ldz, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef sstedc_t *sstedc_f = &_fortran_sstedc

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sstegr "F_FUNC(sstegr,SSTEGR)"(char *jobz, char *range, int *n, s *d, s *e, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, s *z, int *ldz, int *isuppz, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef sstegr_t *sstegr_f = &_fortran_sstegr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sstein "F_FUNC(sstein,SSTEIN)"(int *n, s *d, s *e, int *m, s *w, int *iblock, int *isplit, s *z, int *ldz, s *work, int *iwork, int *ifail, int *info) nogil
cdef sstein_t *sstein_f = &_fortran_sstein

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sstemr "F_FUNC(sstemr,SSTEMR)"(char *jobz, char *range, int *n, s *d, s *e, s *vl, s *vu, int *il, int *iu, int *m, s *w, s *z, int *ldz, int *nzc, int *isuppz, bint *tryrac, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef sstemr_t *sstemr_f = &_fortran_sstemr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssteqr "F_FUNC(ssteqr,SSTEQR)"(char *compz, int *n, s *d, s *e, s *z, int *ldz, s *work, int *info) nogil
cdef ssteqr_t *ssteqr_f = &_fortran_ssteqr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssterf "F_FUNC(ssterf,SSTERF)"(int *n, s *d, s *e, int *info) nogil
cdef ssterf_t *ssterf_f = &_fortran_ssterf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sstev "F_FUNC(sstev,SSTEV)"(char *jobz, int *n, s *d, s *e, s *z, int *ldz, s *work, int *info) nogil
cdef sstev_t *sstev_f = &_fortran_sstev

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sstevd "F_FUNC(sstevd,SSTEVD)"(char *jobz, int *n, s *d, s *e, s *z, int *ldz, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef sstevd_t *sstevd_f = &_fortran_sstevd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sstevr "F_FUNC(sstevr,SSTEVR)"(char *jobz, char *range, int *n, s *d, s *e, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, s *z, int *ldz, int *isuppz, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef sstevr_t *sstevr_f = &_fortran_sstevr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_sstevx "F_FUNC(sstevx,SSTEVX)"(char *jobz, char *range, int *n, s *d, s *e, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, s *z, int *ldz, s *work, int *iwork, int *ifail, int *info) nogil
cdef sstevx_t *sstevx_f = &_fortran_sstevx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssycon "F_FUNC(ssycon,SSYCON)"(char *uplo, int *n, s *a, int *lda, int *ipiv, s *anorm, s *rcond, s *work, int *iwork, int *info) nogil
cdef ssycon_t *ssycon_f = &_fortran_ssycon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssyev "F_FUNC(ssyev,SSYEV)"(char *jobz, char *uplo, int *n, s *a, int *lda, s *w, s *work, int *lwork, int *info) nogil
cdef ssyev_t *ssyev_f = &_fortran_ssyev

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssyevd "F_FUNC(ssyevd,SSYEVD)"(char *jobz, char *uplo, int *n, s *a, int *lda, s *w, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef ssyevd_t *ssyevd_f = &_fortran_ssyevd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssyevr "F_FUNC(ssyevr,SSYEVR)"(char *jobz, char *range, char *uplo, int *n, s *a, int *lda, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, s *z, int *ldz, int *isuppz, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef ssyevr_t *ssyevr_f = &_fortran_ssyevr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssyevx "F_FUNC(ssyevx,SSYEVX)"(char *jobz, char *range, char *uplo, int *n, s *a, int *lda, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, s *z, int *ldz, s *work, int *lwork, int *iwork, int *ifail, int *info) nogil
cdef ssyevx_t *ssyevx_f = &_fortran_ssyevx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssygs2 "F_FUNC(ssygs2,SSYGS2)"(int *itype, char *uplo, int *n, s *a, int *lda, s *b, int *ldb, int *info) nogil
cdef ssygs2_t *ssygs2_f = &_fortran_ssygs2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssygst "F_FUNC(ssygst,SSYGST)"(int *itype, char *uplo, int *n, s *a, int *lda, s *b, int *ldb, int *info) nogil
cdef ssygst_t *ssygst_f = &_fortran_ssygst

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssygv "F_FUNC(ssygv,SSYGV)"(int *itype, char *jobz, char *uplo, int *n, s *a, int *lda, s *b, int *ldb, s *w, s *work, int *lwork, int *info) nogil
cdef ssygv_t *ssygv_f = &_fortran_ssygv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssygvd "F_FUNC(ssygvd,SSYGVD)"(int *itype, char *jobz, char *uplo, int *n, s *a, int *lda, s *b, int *ldb, s *w, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef ssygvd_t *ssygvd_f = &_fortran_ssygvd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssygvx "F_FUNC(ssygvx,SSYGVX)"(int *itype, char *jobz, char *range, char *uplo, int *n, s *a, int *lda, s *b, int *ldb, s *vl, s *vu, int *il, int *iu, s *abstol, int *m, s *w, s *z, int *ldz, s *work, int *lwork, int *iwork, int *ifail, int *info) nogil
cdef ssygvx_t *ssygvx_f = &_fortran_ssygvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssyrfs "F_FUNC(ssyrfs,SSYRFS)"(char *uplo, int *n, int *nrhs, s *a, int *lda, s *af, int *ldaf, int *ipiv, s *b, int *ldb, s *x, int *ldx, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef ssyrfs_t *ssyrfs_f = &_fortran_ssyrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssysv "F_FUNC(ssysv,SSYSV)"(char *uplo, int *n, int *nrhs, s *a, int *lda, int *ipiv, s *b, int *ldb, s *work, int *lwork, int *info) nogil
cdef ssysv_t *ssysv_f = &_fortran_ssysv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssysvx "F_FUNC(ssysvx,SSYSVX)"(char *fact, char *uplo, int *n, int *nrhs, s *a, int *lda, s *af, int *ldaf, int *ipiv, s *b, int *ldb, s *x, int *ldx, s *rcond, s *ferr, s *berr, s *work, int *lwork, int *iwork, int *info) nogil
cdef ssysvx_t *ssysvx_f = &_fortran_ssysvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssytd2 "F_FUNC(ssytd2,SSYTD2)"(char *uplo, int *n, s *a, int *lda, s *d, s *e, s *tau, int *info) nogil
cdef ssytd2_t *ssytd2_f = &_fortran_ssytd2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssytf2 "F_FUNC(ssytf2,SSYTF2)"(char *uplo, int *n, s *a, int *lda, int *ipiv, int *info) nogil
cdef ssytf2_t *ssytf2_f = &_fortran_ssytf2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssytrd "F_FUNC(ssytrd,SSYTRD)"(char *uplo, int *n, s *a, int *lda, s *d, s *e, s *tau, s *work, int *lwork, int *info) nogil
cdef ssytrd_t *ssytrd_f = &_fortran_ssytrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssytrf "F_FUNC(ssytrf,SSYTRF)"(char *uplo, int *n, s *a, int *lda, int *ipiv, s *work, int *lwork, int *info) nogil
cdef ssytrf_t *ssytrf_f = &_fortran_ssytrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssytri "F_FUNC(ssytri,SSYTRI)"(char *uplo, int *n, s *a, int *lda, int *ipiv, s *work, int *info) nogil
cdef ssytri_t *ssytri_f = &_fortran_ssytri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ssytrs "F_FUNC(ssytrs,SSYTRS)"(char *uplo, int *n, int *nrhs, s *a, int *lda, int *ipiv, s *b, int *ldb, int *info) nogil
cdef ssytrs_t *ssytrs_f = &_fortran_ssytrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_stbcon "F_FUNC(stbcon,STBCON)"(char *norm, char *uplo, char *diag, int *n, int *kd, s *ab, int *ldab, s *rcond, s *work, int *iwork, int *info) nogil
cdef stbcon_t *stbcon_f = &_fortran_stbcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_stbrfs "F_FUNC(stbrfs,STBRFS)"(char *uplo, char *trans, char *diag, int *n, int *kd, int *nrhs, s *ab, int *ldab, s *b, int *ldb, s *x, int *ldx, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef stbrfs_t *stbrfs_f = &_fortran_stbrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_stbtrs "F_FUNC(stbtrs,STBTRS)"(char *uplo, char *trans, char *diag, int *n, int *kd, int *nrhs, s *ab, int *ldab, s *b, int *ldb, int *info) nogil
cdef stbtrs_t *stbtrs_f = &_fortran_stbtrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_stgevc "F_FUNC(stgevc,STGEVC)"(char *side, char *howmny, bint *select, int *n, s *s, int *lds, s *p, int *ldp, s *vl, int *ldvl, s *vr, int *ldvr, int *mm, int *m, s *work, int *info) nogil
cdef stgevc_t *stgevc_f = &_fortran_stgevc

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_stgex2 "F_FUNC(stgex2,STGEX2)"(bint *wantq, bint *wantz, int *n, s *a, int *lda, s *b, int *ldb, s *q, int *ldq, s *z, int *ldz, int *j1, int *n1, int *n2, s *work, int *lwork, int *info) nogil
cdef stgex2_t *stgex2_f = &_fortran_stgex2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_stgexc "F_FUNC(stgexc,STGEXC)"(bint *wantq, bint *wantz, int *n, s *a, int *lda, s *b, int *ldb, s *q, int *ldq, s *z, int *ldz, int *ifst, int *ilst, s *work, int *lwork, int *info) nogil
cdef stgexc_t *stgexc_f = &_fortran_stgexc

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_stgsen "F_FUNC(stgsen,STGSEN)"(int *ijob, bint *wantq, bint *wantz, bint *select, int *n, s *a, int *lda, s *b, int *ldb, s *alphar, s *alphai, s *beta, s *q, int *ldq, s *z, int *ldz, int *m, s *pl, s *pr, s *dif, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef stgsen_t *stgsen_f = &_fortran_stgsen

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_stgsja "F_FUNC(stgsja,STGSJA)"(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, int *k, int *l, s *a, int *lda, s *b, int *ldb, s *tola, s *tolb, s *alpha, s *beta, s *u, int *ldu, s *v, int *ldv, s *q, int *ldq, s *work, int *ncycle, int *info) nogil
cdef stgsja_t *stgsja_f = &_fortran_stgsja

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_stgsna "F_FUNC(stgsna,STGSNA)"(char *job, char *howmny, bint *select, int *n, s *a, int *lda, s *b, int *ldb, s *vl, int *ldvl, s *vr, int *ldvr, s *s, s *dif, int *mm, int *m, s *work, int *lwork, int *iwork, int *info) nogil
cdef stgsna_t *stgsna_f = &_fortran_stgsna

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_stgsy2 "F_FUNC(stgsy2,STGSY2)"(char *trans, int *ijob, int *m, int *n, s *a, int *lda, s *b, int *ldb, s *c, int *ldc, s *d, int *ldd, s *e, int *lde, s *f, int *ldf, s *scale, s *rdsum, s *rdscal, int *iwork, int *pq, int *info) nogil
cdef stgsy2_t *stgsy2_f = &_fortran_stgsy2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_stgsyl "F_FUNC(stgsyl,STGSYL)"(char *trans, int *ijob, int *m, int *n, s *a, int *lda, s *b, int *ldb, s *c, int *ldc, s *d, int *ldd, s *e, int *lde, s *f, int *ldf, s *scale, s *dif, s *work, int *lwork, int *iwork, int *info) nogil
cdef stgsyl_t *stgsyl_f = &_fortran_stgsyl

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_stpcon "F_FUNC(stpcon,STPCON)"(char *norm, char *uplo, char *diag, int *n, s *ap, s *rcond, s *work, int *iwork, int *info) nogil
cdef stpcon_t *stpcon_f = &_fortran_stpcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_stprfs "F_FUNC(stprfs,STPRFS)"(char *uplo, char *trans, char *diag, int *n, int *nrhs, s *ap, s *b, int *ldb, s *x, int *ldx, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef stprfs_t *stprfs_f = &_fortran_stprfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_stptri "F_FUNC(stptri,STPTRI)"(char *uplo, char *diag, int *n, s *ap, int *info) nogil
cdef stptri_t *stptri_f = &_fortran_stptri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_stptrs "F_FUNC(stptrs,STPTRS)"(char *uplo, char *trans, char *diag, int *n, int *nrhs, s *ap, s *b, int *ldb, int *info) nogil
cdef stptrs_t *stptrs_f = &_fortran_stptrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_strcon "F_FUNC(strcon,STRCON)"(char *norm, char *uplo, char *diag, int *n, s *a, int *lda, s *rcond, s *work, int *iwork, int *info) nogil
cdef strcon_t *strcon_f = &_fortran_strcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_strevc "F_FUNC(strevc,STREVC)"(char *side, char *howmny, bint *select, int *n, s *t, int *ldt, s *vl, int *ldvl, s *vr, int *ldvr, int *mm, int *m, s *work, int *info) nogil
cdef strevc_t *strevc_f = &_fortran_strevc

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_strexc "F_FUNC(strexc,STREXC)"(char *compq, int *n, s *t, int *ldt, s *q, int *ldq, int *ifst, int *ilst, s *work, int *info) nogil
cdef strexc_t *strexc_f = &_fortran_strexc

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_strrfs "F_FUNC(strrfs,STRRFS)"(char *uplo, char *trans, char *diag, int *n, int *nrhs, s *a, int *lda, s *b, int *ldb, s *x, int *ldx, s *ferr, s *berr, s *work, int *iwork, int *info) nogil
cdef strrfs_t *strrfs_f = &_fortran_strrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_strsen "F_FUNC(strsen,STRSEN)"(char *job, char *compq, bint *select, int *n, s *t, int *ldt, s *q, int *ldq, s *wr, s *wi, int *m, s *s, s *sep, s *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef strsen_t *strsen_f = &_fortran_strsen

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_strsna "F_FUNC(strsna,STRSNA)"(char *job, char *howmny, bint *select, int *n, s *t, int *ldt, s *vl, int *ldvl, s *vr, int *ldvr, s *s, s *sep, int *mm, int *m, s *work, int *ldwork, int *iwork, int *info) nogil
cdef strsna_t *strsna_f = &_fortran_strsna

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_strsyl "F_FUNC(strsyl,STRSYL)"(char *trana, char *tranb, int *isgn, int *m, int *n, s *a, int *lda, s *b, int *ldb, s *c, int *ldc, s *scale, int *info) nogil
cdef strsyl_t *strsyl_f = &_fortran_strsyl

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_strti2 "F_FUNC(strti2,STRTI2)"(char *uplo, char *diag, int *n, s *a, int *lda, int *info) nogil
cdef strti2_t *strti2_f = &_fortran_strti2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_strtri "F_FUNC(strtri,STRTRI)"(char *uplo, char *diag, int *n, s *a, int *lda, int *info) nogil
cdef strtri_t *strtri_f = &_fortran_strtri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_strtrs "F_FUNC(strtrs,STRTRS)"(char *uplo, char *trans, char *diag, int *n, int *nrhs, s *a, int *lda, s *b, int *ldb, int *info) nogil
cdef strtrs_t *strtrs_f = &_fortran_strtrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_stzrqf "F_FUNC(stzrqf,STZRQF)"(int *m, int *n, s *a, int *lda, s *tau, int *info) nogil
cdef stzrqf_t *stzrqf_f = &_fortran_stzrqf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_stzrzf "F_FUNC(stzrzf,STZRZF)"(int *m, int *n, s *a, int *lda, s *tau, s *work, int *lwork, int *info) nogil
cdef stzrzf_t *stzrzf_f = &_fortran_stzrzf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_xerbla "F_FUNC(xerbla,XERBLA)"(char *srname, int *info) nogil
cdef xerbla_t *xerbla_f = &_fortran_xerbla

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zbdsqr "F_FUNC(zbdsqr,ZBDSQR)"(char *uplo, int *n, int *ncvt, int *nru, int *ncc, d *d, d *e, z *vt, int *ldvt, z *u, int *ldu, z *c, int *ldc, d *rwork, int *info) nogil
cdef zbdsqr_t *zbdsqr_f = &_fortran_zbdsqr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zcgesv "F_FUNC(zcgesv,ZCGESV)"(int *n, int *nrhs, z *a, int *lda, int *ipiv, z *b, int *ldb, z *x, int *ldx, z *work, c *swork, d *rwork, int *iter, int *info) nogil
cdef zcgesv_t *zcgesv_f = &_fortran_zcgesv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zdrscl "F_FUNC(zdrscl,ZDRSCL)"(int *n, d *sa, z *sx, int *incx) nogil
cdef zdrscl_t *zdrscl_f = &_fortran_zdrscl

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgbbrd "F_FUNC(zgbbrd,ZGBBRD)"(char *vect, int *m, int *n, int *ncc, int *kl, int *ku, z *ab, int *ldab, d *d, d *e, z *q, int *ldq, z *pt, int *ldpt, z *c, int *ldc, z *work, d *rwork, int *info) nogil
cdef zgbbrd_t *zgbbrd_f = &_fortran_zgbbrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgbcon "F_FUNC(zgbcon,ZGBCON)"(char *norm, int *n, int *kl, int *ku, z *ab, int *ldab, int *ipiv, d *anorm, d *rcond, z *work, d *rwork, int *info) nogil
cdef zgbcon_t *zgbcon_f = &_fortran_zgbcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgbequ "F_FUNC(zgbequ,ZGBEQU)"(int *m, int *n, int *kl, int *ku, z *ab, int *ldab, d *r, d *c, d *rowcnd, d *colcnd, d *amax, int *info) nogil
cdef zgbequ_t *zgbequ_f = &_fortran_zgbequ

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgbrfs "F_FUNC(zgbrfs,ZGBRFS)"(char *trans, int *n, int *kl, int *ku, int *nrhs, z *ab, int *ldab, z *afb, int *ldafb, int *ipiv, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zgbrfs_t *zgbrfs_f = &_fortran_zgbrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgbsv "F_FUNC(zgbsv,ZGBSV)"(int *n, int *kl, int *ku, int *nrhs, z *ab, int *ldab, int *ipiv, z *b, int *ldb, int *info) nogil
cdef zgbsv_t *zgbsv_f = &_fortran_zgbsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgbsvx "F_FUNC(zgbsvx,ZGBSVX)"(char *fact, char *trans, int *n, int *kl, int *ku, int *nrhs, z *ab, int *ldab, z *afb, int *ldafb, int *ipiv, char *equed, d *r, d *c, z *b, int *ldb, z *x, int *ldx, d *rcond, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zgbsvx_t *zgbsvx_f = &_fortran_zgbsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgbtf2 "F_FUNC(zgbtf2,ZGBTF2)"(int *m, int *n, int *kl, int *ku, z *ab, int *ldab, int *ipiv, int *info) nogil
cdef zgbtf2_t *zgbtf2_f = &_fortran_zgbtf2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgbtrf "F_FUNC(zgbtrf,ZGBTRF)"(int *m, int *n, int *kl, int *ku, z *ab, int *ldab, int *ipiv, int *info) nogil
cdef zgbtrf_t *zgbtrf_f = &_fortran_zgbtrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgbtrs "F_FUNC(zgbtrs,ZGBTRS)"(char *trans, int *n, int *kl, int *ku, int *nrhs, z *ab, int *ldab, int *ipiv, z *b, int *ldb, int *info) nogil
cdef zgbtrs_t *zgbtrs_f = &_fortran_zgbtrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgebak "F_FUNC(zgebak,ZGEBAK)"(char *job, char *side, int *n, int *ilo, int *ihi, d *scale, int *m, z *v, int *ldv, int *info) nogil
cdef zgebak_t *zgebak_f = &_fortran_zgebak

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgebal "F_FUNC(zgebal,ZGEBAL)"(char *job, int *n, z *a, int *lda, int *ilo, int *ihi, d *scale, int *info) nogil
cdef zgebal_t *zgebal_f = &_fortran_zgebal

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgebd2 "F_FUNC(zgebd2,ZGEBD2)"(int *m, int *n, z *a, int *lda, d *d, d *e, z *tauq, z *taup, z *work, int *info) nogil
cdef zgebd2_t *zgebd2_f = &_fortran_zgebd2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgebrd "F_FUNC(zgebrd,ZGEBRD)"(int *m, int *n, z *a, int *lda, d *d, d *e, z *tauq, z *taup, z *work, int *lwork, int *info) nogil
cdef zgebrd_t *zgebrd_f = &_fortran_zgebrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgecon "F_FUNC(zgecon,ZGECON)"(char *norm, int *n, z *a, int *lda, d *anorm, d *rcond, z *work, d *rwork, int *info) nogil
cdef zgecon_t *zgecon_f = &_fortran_zgecon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgeequ "F_FUNC(zgeequ,ZGEEQU)"(int *m, int *n, z *a, int *lda, d *r, d *c, d *rowcnd, d *colcnd, d *amax, int *info) nogil
cdef zgeequ_t *zgeequ_f = &_fortran_zgeequ

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgees "F_FUNC(zgees,ZGEES)"(char *jobvs, char *sort, zselect1 *select, int *n, z *a, int *lda, int *sdim, z *w, z *vs, int *ldvs, z *work, int *lwork, d *rwork, bint *bwork, int *info) nogil
cdef zgees_t *zgees_f = &_fortran_zgees

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgeesx "F_FUNC(zgeesx,ZGEESX)"(char *jobvs, char *sort, zselect1 *select, char *sense, int *n, z *a, int *lda, int *sdim, z *w, z *vs, int *ldvs, d *rconde, d *rcondv, z *work, int *lwork, d *rwork, bint *bwork, int *info) nogil
cdef zgeesx_t *zgeesx_f = &_fortran_zgeesx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgeev "F_FUNC(zgeev,ZGEEV)"(char *jobvl, char *jobvr, int *n, z *a, int *lda, z *w, z *vl, int *ldvl, z *vr, int *ldvr, z *work, int *lwork, d *rwork, int *info) nogil
cdef zgeev_t *zgeev_f = &_fortran_zgeev

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgeevx "F_FUNC(zgeevx,ZGEEVX)"(char *balanc, char *jobvl, char *jobvr, char *sense, int *n, z *a, int *lda, z *w, z *vl, int *ldvl, z *vr, int *ldvr, int *ilo, int *ihi, d *scale, d *abnrm, d *rconde, d *rcondv, z *work, int *lwork, d *rwork, int *info) nogil
cdef zgeevx_t *zgeevx_f = &_fortran_zgeevx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgegs "F_FUNC(zgegs,ZGEGS)"(char *jobvsl, char *jobvsr, int *n, z *a, int *lda, z *b, int *ldb, z *alpha, z *beta, z *vsl, int *ldvsl, z *vsr, int *ldvsr, z *work, int *lwork, d *rwork, int *info) nogil
cdef zgegs_t *zgegs_f = &_fortran_zgegs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgegv "F_FUNC(zgegv,ZGEGV)"(char *jobvl, char *jobvr, int *n, z *a, int *lda, z *b, int *ldb, z *alpha, z *beta, z *vl, int *ldvl, z *vr, int *ldvr, z *work, int *lwork, d *rwork, int *info) nogil
cdef zgegv_t *zgegv_f = &_fortran_zgegv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgehd2 "F_FUNC(zgehd2,ZGEHD2)"(int *n, int *ilo, int *ihi, z *a, int *lda, z *tau, z *work, int *info) nogil
cdef zgehd2_t *zgehd2_f = &_fortran_zgehd2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgehrd "F_FUNC(zgehrd,ZGEHRD)"(int *n, int *ilo, int *ihi, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef zgehrd_t *zgehrd_f = &_fortran_zgehrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgelq2 "F_FUNC(zgelq2,ZGELQ2)"(int *m, int *n, z *a, int *lda, z *tau, z *work, int *info) nogil
cdef zgelq2_t *zgelq2_f = &_fortran_zgelq2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgelqf "F_FUNC(zgelqf,ZGELQF)"(int *m, int *n, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef zgelqf_t *zgelqf_f = &_fortran_zgelqf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgels "F_FUNC(zgels,ZGELS)"(char *trans, int *m, int *n, int *nrhs, z *a, int *lda, z *b, int *ldb, z *work, int *lwork, int *info) nogil
cdef zgels_t *zgels_f = &_fortran_zgels

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgelsd "F_FUNC(zgelsd,ZGELSD)"(int *m, int *n, int *nrhs, z *a, int *lda, z *b, int *ldb, d *s, d *rcond, int *rank, z *work, int *lwork, d *rwork, int *iwork, int *info) nogil
cdef zgelsd_t *zgelsd_f = &_fortran_zgelsd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgelss "F_FUNC(zgelss,ZGELSS)"(int *m, int *n, int *nrhs, z *a, int *lda, z *b, int *ldb, d *s, d *rcond, int *rank, z *work, int *lwork, d *rwork, int *info) nogil
cdef zgelss_t *zgelss_f = &_fortran_zgelss

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgelsx "F_FUNC(zgelsx,ZGELSX)"(int *m, int *n, int *nrhs, z *a, int *lda, z *b, int *ldb, int *jpvt, d *rcond, int *rank, z *work, d *rwork, int *info) nogil
cdef zgelsx_t *zgelsx_f = &_fortran_zgelsx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgelsy "F_FUNC(zgelsy,ZGELSY)"(int *m, int *n, int *nrhs, z *a, int *lda, z *b, int *ldb, int *jpvt, d *rcond, int *rank, z *work, int *lwork, d *rwork, int *info) nogil
cdef zgelsy_t *zgelsy_f = &_fortran_zgelsy

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgeql2 "F_FUNC(zgeql2,ZGEQL2)"(int *m, int *n, z *a, int *lda, z *tau, z *work, int *info) nogil
cdef zgeql2_t *zgeql2_f = &_fortran_zgeql2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgeqlf "F_FUNC(zgeqlf,ZGEQLF)"(int *m, int *n, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef zgeqlf_t *zgeqlf_f = &_fortran_zgeqlf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgeqp3 "F_FUNC(zgeqp3,ZGEQP3)"(int *m, int *n, z *a, int *lda, int *jpvt, z *tau, z *work, int *lwork, d *rwork, int *info) nogil
cdef zgeqp3_t *zgeqp3_f = &_fortran_zgeqp3

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgeqpf "F_FUNC(zgeqpf,ZGEQPF)"(int *m, int *n, z *a, int *lda, int *jpvt, z *tau, z *work, d *rwork, int *info) nogil
cdef zgeqpf_t *zgeqpf_f = &_fortran_zgeqpf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgeqr2 "F_FUNC(zgeqr2,ZGEQR2)"(int *m, int *n, z *a, int *lda, z *tau, z *work, int *info) nogil
cdef zgeqr2_t *zgeqr2_f = &_fortran_zgeqr2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgeqrf "F_FUNC(zgeqrf,ZGEQRF)"(int *m, int *n, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef zgeqrf_t *zgeqrf_f = &_fortran_zgeqrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgerfs "F_FUNC(zgerfs,ZGERFS)"(char *trans, int *n, int *nrhs, z *a, int *lda, z *af, int *ldaf, int *ipiv, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zgerfs_t *zgerfs_f = &_fortran_zgerfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgerq2 "F_FUNC(zgerq2,ZGERQ2)"(int *m, int *n, z *a, int *lda, z *tau, z *work, int *info) nogil
cdef zgerq2_t *zgerq2_f = &_fortran_zgerq2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgerqf "F_FUNC(zgerqf,ZGERQF)"(int *m, int *n, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef zgerqf_t *zgerqf_f = &_fortran_zgerqf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgesc2 "F_FUNC(zgesc2,ZGESC2)"(int *n, z *a, int *lda, z *rhs, int *ipiv, int *jpiv, d *scale) nogil
cdef zgesc2_t *zgesc2_f = &_fortran_zgesc2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgesdd "F_FUNC(zgesdd,ZGESDD)"(char *jobz, int *m, int *n, z *a, int *lda, d *s, z *u, int *ldu, z *vt, int *ldvt, z *work, int *lwork, d *rwork, int *iwork, int *info) nogil
cdef zgesdd_t *zgesdd_f = &_fortran_zgesdd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgesv "F_FUNC(zgesv,ZGESV)"(int *n, int *nrhs, z *a, int *lda, int *ipiv, z *b, int *ldb, int *info) nogil
cdef zgesv_t *zgesv_f = &_fortran_zgesv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgesvd "F_FUNC(zgesvd,ZGESVD)"(char *jobu, char *jobvt, int *m, int *n, z *a, int *lda, d *s, z *u, int *ldu, z *vt, int *ldvt, z *work, int *lwork, d *rwork, int *info) nogil
cdef zgesvd_t *zgesvd_f = &_fortran_zgesvd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgesvx "F_FUNC(zgesvx,ZGESVX)"(char *fact, char *trans, int *n, int *nrhs, z *a, int *lda, z *af, int *ldaf, int *ipiv, char *equed, d *r, d *c, z *b, int *ldb, z *x, int *ldx, d *rcond, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zgesvx_t *zgesvx_f = &_fortran_zgesvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgetc2 "F_FUNC(zgetc2,ZGETC2)"(int *n, z *a, int *lda, int *ipiv, int *jpiv, int *info) nogil
cdef zgetc2_t *zgetc2_f = &_fortran_zgetc2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgetf2 "F_FUNC(zgetf2,ZGETF2)"(int *m, int *n, z *a, int *lda, int *ipiv, int *info) nogil
cdef zgetf2_t *zgetf2_f = &_fortran_zgetf2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgetrf "F_FUNC(zgetrf,ZGETRF)"(int *m, int *n, z *a, int *lda, int *ipiv, int *info) nogil
cdef zgetrf_t *zgetrf_f = &_fortran_zgetrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgetri "F_FUNC(zgetri,ZGETRI)"(int *n, z *a, int *lda, int *ipiv, z *work, int *lwork, int *info) nogil
cdef zgetri_t *zgetri_f = &_fortran_zgetri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgetrs "F_FUNC(zgetrs,ZGETRS)"(char *trans, int *n, int *nrhs, z *a, int *lda, int *ipiv, z *b, int *ldb, int *info) nogil
cdef zgetrs_t *zgetrs_f = &_fortran_zgetrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zggbak "F_FUNC(zggbak,ZGGBAK)"(char *job, char *side, int *n, int *ilo, int *ihi, d *lscale, d *rscale, int *m, z *v, int *ldv, int *info) nogil
cdef zggbak_t *zggbak_f = &_fortran_zggbak

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zggbal "F_FUNC(zggbal,ZGGBAL)"(char *job, int *n, z *a, int *lda, z *b, int *ldb, int *ilo, int *ihi, d *lscale, d *rscale, d *work, int *info) nogil
cdef zggbal_t *zggbal_f = &_fortran_zggbal

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgges "F_FUNC(zgges,ZGGES)"(char *jobvsl, char *jobvsr, char *sort, zselect2 *selctg, int *n, z *a, int *lda, z *b, int *ldb, int *sdim, z *alpha, z *beta, z *vsl, int *ldvsl, z *vsr, int *ldvsr, z *work, int *lwork, d *rwork, bint *bwork, int *info) nogil
cdef zgges_t *zgges_f = &_fortran_zgges

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zggesx "F_FUNC(zggesx,ZGGESX)"(char *jobvsl, char *jobvsr, char *sort, zselect2 *selctg, char *sense, int *n, z *a, int *lda, z *b, int *ldb, int *sdim, z *alpha, z *beta, z *vsl, int *ldvsl, z *vsr, int *ldvsr, d *rconde, d *rcondv, z *work, int *lwork, d *rwork, int *iwork, int *liwork, bint *bwork, int *info) nogil
cdef zggesx_t *zggesx_f = &_fortran_zggesx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zggev "F_FUNC(zggev,ZGGEV)"(char *jobvl, char *jobvr, int *n, z *a, int *lda, z *b, int *ldb, z *alpha, z *beta, z *vl, int *ldvl, z *vr, int *ldvr, z *work, int *lwork, d *rwork, int *info) nogil
cdef zggev_t *zggev_f = &_fortran_zggev

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zggevx "F_FUNC(zggevx,ZGGEVX)"(char *balanc, char *jobvl, char *jobvr, char *sense, int *n, z *a, int *lda, z *b, int *ldb, z *alpha, z *beta, z *vl, int *ldvl, z *vr, int *ldvr, int *ilo, int *ihi, d *lscale, d *rscale, d *abnrm, d *bbnrm, d *rconde, d *rcondv, z *work, int *lwork, d *rwork, int *iwork, bint *bwork, int *info) nogil
cdef zggevx_t *zggevx_f = &_fortran_zggevx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zggglm "F_FUNC(zggglm,ZGGGLM)"(int *n, int *m, int *p, z *a, int *lda, z *b, int *ldb, z *d, z *x, z *y, z *work, int *lwork, int *info) nogil
cdef zggglm_t *zggglm_f = &_fortran_zggglm

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgghrd "F_FUNC(zgghrd,ZGGHRD)"(char *compq, char *compz, int *n, int *ilo, int *ihi, z *a, int *lda, z *b, int *ldb, z *q, int *ldq, z *z, int *ldz, int *info) nogil
cdef zgghrd_t *zgghrd_f = &_fortran_zgghrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgglse "F_FUNC(zgglse,ZGGLSE)"(int *m, int *n, int *p, z *a, int *lda, z *b, int *ldb, z *c, z *d, z *x, z *work, int *lwork, int *info) nogil
cdef zgglse_t *zgglse_f = &_fortran_zgglse

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zggqrf "F_FUNC(zggqrf,ZGGQRF)"(int *n, int *m, int *p, z *a, int *lda, z *taua, z *b, int *ldb, z *taub, z *work, int *lwork, int *info) nogil
cdef zggqrf_t *zggqrf_f = &_fortran_zggqrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zggrqf "F_FUNC(zggrqf,ZGGRQF)"(int *m, int *p, int *n, z *a, int *lda, z *taua, z *b, int *ldb, z *taub, z *work, int *lwork, int *info) nogil
cdef zggrqf_t *zggrqf_f = &_fortran_zggrqf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zggsvd "F_FUNC(zggsvd,ZGGSVD)"(char *jobu, char *jobv, char *jobq, int *m, int *n, int *p, int *k, int *l, z *a, int *lda, z *b, int *ldb, d *alpha, d *beta, z *u, int *ldu, z *v, int *ldv, z *q, int *ldq, z *work, d *rwork, int *iwork, int *info) nogil
cdef zggsvd_t *zggsvd_f = &_fortran_zggsvd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zggsvp "F_FUNC(zggsvp,ZGGSVP)"(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, z *a, int *lda, z *b, int *ldb, d *tola, d *tolb, int *k, int *l, z *u, int *ldu, z *v, int *ldv, z *q, int *ldq, int *iwork, d *rwork, z *tau, z *work, int *info) nogil
cdef zggsvp_t *zggsvp_f = &_fortran_zggsvp

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgtcon "F_FUNC(zgtcon,ZGTCON)"(char *norm, int *n, z *dl, z *d, z *du, z *du2, int *ipiv, d *anorm, d *rcond, z *work, int *info) nogil
cdef zgtcon_t *zgtcon_f = &_fortran_zgtcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgtrfs "F_FUNC(zgtrfs,ZGTRFS)"(char *trans, int *n, int *nrhs, z *dl, z *d, z *du, z *dlf, z *df, z *duf, z *du2, int *ipiv, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zgtrfs_t *zgtrfs_f = &_fortran_zgtrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgtsv "F_FUNC(zgtsv,ZGTSV)"(int *n, int *nrhs, z *dl, z *d, z *du, z *b, int *ldb, int *info) nogil
cdef zgtsv_t *zgtsv_f = &_fortran_zgtsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgtsvx "F_FUNC(zgtsvx,ZGTSVX)"(char *fact, char *trans, int *n, int *nrhs, z *dl, z *d, z *du, z *dlf, z *df, z *duf, z *du2, int *ipiv, z *b, int *ldb, z *x, int *ldx, d *rcond, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zgtsvx_t *zgtsvx_f = &_fortran_zgtsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgttrf "F_FUNC(zgttrf,ZGTTRF)"(int *n, z *dl, z *d, z *du, z *du2, int *ipiv, int *info) nogil
cdef zgttrf_t *zgttrf_f = &_fortran_zgttrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgttrs "F_FUNC(zgttrs,ZGTTRS)"(char *trans, int *n, int *nrhs, z *dl, z *d, z *du, z *du2, int *ipiv, z *b, int *ldb, int *info) nogil
cdef zgttrs_t *zgttrs_f = &_fortran_zgttrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zgtts2 "F_FUNC(zgtts2,ZGTTS2)"(int *itrans, int *n, int *nrhs, z *dl, z *d, z *du, z *du2, int *ipiv, z *b, int *ldb) nogil
cdef zgtts2_t *zgtts2_f = &_fortran_zgtts2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhbev "F_FUNC(zhbev,ZHBEV)"(char *jobz, char *uplo, int *n, int *kd, z *ab, int *ldab, d *w, z *z, int *ldz, z *work, d *rwork, int *info) nogil
cdef zhbev_t *zhbev_f = &_fortran_zhbev

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhbevd "F_FUNC(zhbevd,ZHBEVD)"(char *jobz, char *uplo, int *n, int *kd, z *ab, int *ldab, d *w, z *z, int *ldz, z *work, int *lwork, d *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef zhbevd_t *zhbevd_f = &_fortran_zhbevd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhbevx "F_FUNC(zhbevx,ZHBEVX)"(char *jobz, char *range, char *uplo, int *n, int *kd, z *ab, int *ldab, z *q, int *ldq, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, z *z, int *ldz, z *work, d *rwork, int *iwork, int *ifail, int *info) nogil
cdef zhbevx_t *zhbevx_f = &_fortran_zhbevx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhbgst "F_FUNC(zhbgst,ZHBGST)"(char *vect, char *uplo, int *n, int *ka, int *kb, z *ab, int *ldab, z *bb, int *ldbb, z *x, int *ldx, z *work, d *rwork, int *info) nogil
cdef zhbgst_t *zhbgst_f = &_fortran_zhbgst

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhbgv "F_FUNC(zhbgv,ZHBGV)"(char *jobz, char *uplo, int *n, int *ka, int *kb, z *ab, int *ldab, z *bb, int *ldbb, d *w, z *z, int *ldz, z *work, d *rwork, int *info) nogil
cdef zhbgv_t *zhbgv_f = &_fortran_zhbgv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhbgvd "F_FUNC(zhbgvd,ZHBGVD)"(char *jobz, char *uplo, int *n, int *ka, int *kb, z *ab, int *ldab, z *bb, int *ldbb, d *w, z *z, int *ldz, z *work, int *lwork, d *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef zhbgvd_t *zhbgvd_f = &_fortran_zhbgvd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhbgvx "F_FUNC(zhbgvx,ZHBGVX)"(char *jobz, char *range, char *uplo, int *n, int *ka, int *kb, z *ab, int *ldab, z *bb, int *ldbb, z *q, int *ldq, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, z *z, int *ldz, z *work, d *rwork, int *iwork, int *ifail, int *info) nogil
cdef zhbgvx_t *zhbgvx_f = &_fortran_zhbgvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhbtrd "F_FUNC(zhbtrd,ZHBTRD)"(char *vect, char *uplo, int *n, int *kd, z *ab, int *ldab, d *d, d *e, z *q, int *ldq, z *work, int *info) nogil
cdef zhbtrd_t *zhbtrd_f = &_fortran_zhbtrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhecon "F_FUNC(zhecon,ZHECON)"(char *uplo, int *n, z *a, int *lda, int *ipiv, d *anorm, d *rcond, z *work, int *info) nogil
cdef zhecon_t *zhecon_f = &_fortran_zhecon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zheev "F_FUNC(zheev,ZHEEV)"(char *jobz, char *uplo, int *n, z *a, int *lda, d *w, z *work, int *lwork, d *rwork, int *info) nogil
cdef zheev_t *zheev_f = &_fortran_zheev

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zheevd "F_FUNC(zheevd,ZHEEVD)"(char *jobz, char *uplo, int *n, z *a, int *lda, d *w, z *work, int *lwork, d *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef zheevd_t *zheevd_f = &_fortran_zheevd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zheevr "F_FUNC(zheevr,ZHEEVR)"(char *jobz, char *range, char *uplo, int *n, z *a, int *lda, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, z *z, int *ldz, int *isuppz, z *work, int *lwork, d *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef zheevr_t *zheevr_f = &_fortran_zheevr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zheevx "F_FUNC(zheevx,ZHEEVX)"(char *jobz, char *range, char *uplo, int *n, z *a, int *lda, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, z *z, int *ldz, z *work, int *lwork, d *rwork, int *iwork, int *ifail, int *info) nogil
cdef zheevx_t *zheevx_f = &_fortran_zheevx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhegs2 "F_FUNC(zhegs2,ZHEGS2)"(int *itype, char *uplo, int *n, z *a, int *lda, z *b, int *ldb, int *info) nogil
cdef zhegs2_t *zhegs2_f = &_fortran_zhegs2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhegst "F_FUNC(zhegst,ZHEGST)"(int *itype, char *uplo, int *n, z *a, int *lda, z *b, int *ldb, int *info) nogil
cdef zhegst_t *zhegst_f = &_fortran_zhegst

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhegv "F_FUNC(zhegv,ZHEGV)"(int *itype, char *jobz, char *uplo, int *n, z *a, int *lda, z *b, int *ldb, d *w, z *work, int *lwork, d *rwork, int *info) nogil
cdef zhegv_t *zhegv_f = &_fortran_zhegv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhegvd "F_FUNC(zhegvd,ZHEGVD)"(int *itype, char *jobz, char *uplo, int *n, z *a, int *lda, z *b, int *ldb, d *w, z *work, int *lwork, d *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef zhegvd_t *zhegvd_f = &_fortran_zhegvd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhegvx "F_FUNC(zhegvx,ZHEGVX)"(int *itype, char *jobz, char *range, char *uplo, int *n, z *a, int *lda, z *b, int *ldb, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, z *z, int *ldz, z *work, int *lwork, d *rwork, int *iwork, int *ifail, int *info) nogil
cdef zhegvx_t *zhegvx_f = &_fortran_zhegvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zherfs "F_FUNC(zherfs,ZHERFS)"(char *uplo, int *n, int *nrhs, z *a, int *lda, z *af, int *ldaf, int *ipiv, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zherfs_t *zherfs_f = &_fortran_zherfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhesv "F_FUNC(zhesv,ZHESV)"(char *uplo, int *n, int *nrhs, z *a, int *lda, int *ipiv, z *b, int *ldb, z *work, int *lwork, int *info) nogil
cdef zhesv_t *zhesv_f = &_fortran_zhesv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhesvx "F_FUNC(zhesvx,ZHESVX)"(char *fact, char *uplo, int *n, int *nrhs, z *a, int *lda, z *af, int *ldaf, int *ipiv, z *b, int *ldb, z *x, int *ldx, d *rcond, d *ferr, d *berr, z *work, int *lwork, d *rwork, int *info) nogil
cdef zhesvx_t *zhesvx_f = &_fortran_zhesvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhetd2 "F_FUNC(zhetd2,ZHETD2)"(char *uplo, int *n, z *a, int *lda, d *d, d *e, z *tau, int *info) nogil
cdef zhetd2_t *zhetd2_f = &_fortran_zhetd2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhetf2 "F_FUNC(zhetf2,ZHETF2)"(char *uplo, int *n, z *a, int *lda, int *ipiv, int *info) nogil
cdef zhetf2_t *zhetf2_f = &_fortran_zhetf2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhetrd "F_FUNC(zhetrd,ZHETRD)"(char *uplo, int *n, z *a, int *lda, d *d, d *e, z *tau, z *work, int *lwork, int *info) nogil
cdef zhetrd_t *zhetrd_f = &_fortran_zhetrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhetrf "F_FUNC(zhetrf,ZHETRF)"(char *uplo, int *n, z *a, int *lda, int *ipiv, z *work, int *lwork, int *info) nogil
cdef zhetrf_t *zhetrf_f = &_fortran_zhetrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhetri "F_FUNC(zhetri,ZHETRI)"(char *uplo, int *n, z *a, int *lda, int *ipiv, z *work, int *info) nogil
cdef zhetri_t *zhetri_f = &_fortran_zhetri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhetrs "F_FUNC(zhetrs,ZHETRS)"(char *uplo, int *n, int *nrhs, z *a, int *lda, int *ipiv, z *b, int *ldb, int *info) nogil
cdef zhetrs_t *zhetrs_f = &_fortran_zhetrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhgeqz "F_FUNC(zhgeqz,ZHGEQZ)"(char *job, char *compq, char *compz, int *n, int *ilo, int *ihi, z *h, int *ldh, z *t, int *ldt, z *alpha, z *beta, z *q, int *ldq, z *z, int *ldz, z *work, int *lwork, d *rwork, int *info) nogil
cdef zhgeqz_t *zhgeqz_f = &_fortran_zhgeqz

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhpcon "F_FUNC(zhpcon,ZHPCON)"(char *uplo, int *n, z *ap, int *ipiv, d *anorm, d *rcond, z *work, int *info) nogil
cdef zhpcon_t *zhpcon_f = &_fortran_zhpcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhpev "F_FUNC(zhpev,ZHPEV)"(char *jobz, char *uplo, int *n, z *ap, d *w, z *z, int *ldz, z *work, d *rwork, int *info) nogil
cdef zhpev_t *zhpev_f = &_fortran_zhpev

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhpevd "F_FUNC(zhpevd,ZHPEVD)"(char *jobz, char *uplo, int *n, z *ap, d *w, z *z, int *ldz, z *work, int *lwork, d *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef zhpevd_t *zhpevd_f = &_fortran_zhpevd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhpevx "F_FUNC(zhpevx,ZHPEVX)"(char *jobz, char *range, char *uplo, int *n, z *ap, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, z *z, int *ldz, z *work, d *rwork, int *iwork, int *ifail, int *info) nogil
cdef zhpevx_t *zhpevx_f = &_fortran_zhpevx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhpgst "F_FUNC(zhpgst,ZHPGST)"(int *itype, char *uplo, int *n, z *ap, z *bp, int *info) nogil
cdef zhpgst_t *zhpgst_f = &_fortran_zhpgst

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhpgv "F_FUNC(zhpgv,ZHPGV)"(int *itype, char *jobz, char *uplo, int *n, z *ap, z *bp, d *w, z *z, int *ldz, z *work, d *rwork, int *info) nogil
cdef zhpgv_t *zhpgv_f = &_fortran_zhpgv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhpgvd "F_FUNC(zhpgvd,ZHPGVD)"(int *itype, char *jobz, char *uplo, int *n, z *ap, z *bp, d *w, z *z, int *ldz, z *work, int *lwork, d *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef zhpgvd_t *zhpgvd_f = &_fortran_zhpgvd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhpgvx "F_FUNC(zhpgvx,ZHPGVX)"(int *itype, char *jobz, char *range, char *uplo, int *n, z *ap, z *bp, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, z *z, int *ldz, z *work, d *rwork, int *iwork, int *ifail, int *info) nogil
cdef zhpgvx_t *zhpgvx_f = &_fortran_zhpgvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhprfs "F_FUNC(zhprfs,ZHPRFS)"(char *uplo, int *n, int *nrhs, z *ap, z *afp, int *ipiv, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zhprfs_t *zhprfs_f = &_fortran_zhprfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhpsv "F_FUNC(zhpsv,ZHPSV)"(char *uplo, int *n, int *nrhs, z *ap, int *ipiv, z *b, int *ldb, int *info) nogil
cdef zhpsv_t *zhpsv_f = &_fortran_zhpsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhpsvx "F_FUNC(zhpsvx,ZHPSVX)"(char *fact, char *uplo, int *n, int *nrhs, z *ap, z *afp, int *ipiv, z *b, int *ldb, z *x, int *ldx, d *rcond, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zhpsvx_t *zhpsvx_f = &_fortran_zhpsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhptrd "F_FUNC(zhptrd,ZHPTRD)"(char *uplo, int *n, z *ap, d *d, d *e, z *tau, int *info) nogil
cdef zhptrd_t *zhptrd_f = &_fortran_zhptrd

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhptrf "F_FUNC(zhptrf,ZHPTRF)"(char *uplo, int *n, z *ap, int *ipiv, int *info) nogil
cdef zhptrf_t *zhptrf_f = &_fortran_zhptrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhptri "F_FUNC(zhptri,ZHPTRI)"(char *uplo, int *n, z *ap, int *ipiv, z *work, int *info) nogil
cdef zhptri_t *zhptri_f = &_fortran_zhptri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhptrs "F_FUNC(zhptrs,ZHPTRS)"(char *uplo, int *n, int *nrhs, z *ap, int *ipiv, z *b, int *ldb, int *info) nogil
cdef zhptrs_t *zhptrs_f = &_fortran_zhptrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhsein "F_FUNC(zhsein,ZHSEIN)"(char *side, char *eigsrc, char *initv, bint *select, int *n, z *h, int *ldh, z *w, z *vl, int *ldvl, z *vr, int *ldvr, int *mm, int *m, z *work, d *rwork, int *ifaill, int *ifailr, int *info) nogil
cdef zhsein_t *zhsein_f = &_fortran_zhsein

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zhseqr "F_FUNC(zhseqr,ZHSEQR)"(char *job, char *compz, int *n, int *ilo, int *ihi, z *h, int *ldh, z *w, z *z, int *ldz, z *work, int *lwork, int *info) nogil
cdef zhseqr_t *zhseqr_f = &_fortran_zhseqr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zlacn2 "F_FUNC(zlacn2,ZLACN2)"(int *n, z *v, z *x, d *est, int *kase, int *isave) nogil
cdef zlacn2_t *zlacn2_f = &_fortran_zlacn2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zlacon "F_FUNC(zlacon,ZLACON)"(int *n, z *v, z *x, d *est, int *kase) nogil
cdef zlacon_t *zlacon_f = &_fortran_zlacon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zlarf "F_FUNC(zlarf,ZLARF)"(char *side, int *m, int *n, z *v, int *incv, z *tau, z *c, int *ldc, z *work) nogil
cdef zlarf_t *zlarf_f = &_fortran_zlarf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zlarz "F_FUNC(zlarz,ZLARZ)"(char *side, int *m, int *n, int *l, z *v, int *incv, z *tau, z *c, int *ldc, z *work) nogil
cdef zlarz_t *zlarz_f = &_fortran_zlarz

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zlaswp "F_FUNC(zlaswp,ZLASWP)"(int *n, z *a, int *lda, int *k1, int *k2, int *ipiv, int *incx) nogil
cdef zlaswp_t *zlaswp_f = &_fortran_zlaswp

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zlauum "F_FUNC(zlauum,ZLAUUM)"(char *uplo, int *n, z *a, int *lda, int *info) nogil
cdef zlauum_t *zlauum_f = &_fortran_zlauum

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zpbcon "F_FUNC(zpbcon,ZPBCON)"(char *uplo, int *n, int *kd, z *ab, int *ldab, d *anorm, d *rcond, z *work, d *rwork, int *info) nogil
cdef zpbcon_t *zpbcon_f = &_fortran_zpbcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zpbequ "F_FUNC(zpbequ,ZPBEQU)"(char *uplo, int *n, int *kd, z *ab, int *ldab, d *s, d *scond, d *amax, int *info) nogil
cdef zpbequ_t *zpbequ_f = &_fortran_zpbequ

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zpbrfs "F_FUNC(zpbrfs,ZPBRFS)"(char *uplo, int *n, int *kd, int *nrhs, z *ab, int *ldab, z *afb, int *ldafb, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zpbrfs_t *zpbrfs_f = &_fortran_zpbrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zpbstf "F_FUNC(zpbstf,ZPBSTF)"(char *uplo, int *n, int *kd, z *ab, int *ldab, int *info) nogil
cdef zpbstf_t *zpbstf_f = &_fortran_zpbstf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zpbsv "F_FUNC(zpbsv,ZPBSV)"(char *uplo, int *n, int *kd, int *nrhs, z *ab, int *ldab, z *b, int *ldb, int *info) nogil
cdef zpbsv_t *zpbsv_f = &_fortran_zpbsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zpbsvx "F_FUNC(zpbsvx,ZPBSVX)"(char *fact, char *uplo, int *n, int *kd, int *nrhs, z *ab, int *ldab, z *afb, int *ldafb, char *equed, d *s, z *b, int *ldb, z *x, int *ldx, d *rcond, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zpbsvx_t *zpbsvx_f = &_fortran_zpbsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zpbtf2 "F_FUNC(zpbtf2,ZPBTF2)"(char *uplo, int *n, int *kd, z *ab, int *ldab, int *info) nogil
cdef zpbtf2_t *zpbtf2_f = &_fortran_zpbtf2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zpbtrf "F_FUNC(zpbtrf,ZPBTRF)"(char *uplo, int *n, int *kd, z *ab, int *ldab, int *info) nogil
cdef zpbtrf_t *zpbtrf_f = &_fortran_zpbtrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zpbtrs "F_FUNC(zpbtrs,ZPBTRS)"(char *uplo, int *n, int *kd, int *nrhs, z *ab, int *ldab, z *b, int *ldb, int *info) nogil
cdef zpbtrs_t *zpbtrs_f = &_fortran_zpbtrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zpocon "F_FUNC(zpocon,ZPOCON)"(char *uplo, int *n, z *a, int *lda, d *anorm, d *rcond, z *work, d *rwork, int *info) nogil
cdef zpocon_t *zpocon_f = &_fortran_zpocon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zpoequ "F_FUNC(zpoequ,ZPOEQU)"(int *n, z *a, int *lda, d *s, d *scond, d *amax, int *info) nogil
cdef zpoequ_t *zpoequ_f = &_fortran_zpoequ

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zporfs "F_FUNC(zporfs,ZPORFS)"(char *uplo, int *n, int *nrhs, z *a, int *lda, z *af, int *ldaf, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zporfs_t *zporfs_f = &_fortran_zporfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zposv "F_FUNC(zposv,ZPOSV)"(char *uplo, int *n, int *nrhs, z *a, int *lda, z *b, int *ldb, int *info) nogil
cdef zposv_t *zposv_f = &_fortran_zposv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zposvx "F_FUNC(zposvx,ZPOSVX)"(char *fact, char *uplo, int *n, int *nrhs, z *a, int *lda, z *af, int *ldaf, char *equed, d *s, z *b, int *ldb, z *x, int *ldx, d *rcond, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zposvx_t *zposvx_f = &_fortran_zposvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zpotf2 "F_FUNC(zpotf2,ZPOTF2)"(char *uplo, int *n, z *a, int *lda, int *info) nogil
cdef zpotf2_t *zpotf2_f = &_fortran_zpotf2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zpotrf "F_FUNC(zpotrf,ZPOTRF)"(char *uplo, int *n, z *a, int *lda, int *info) nogil
cdef zpotrf_t *zpotrf_f = &_fortran_zpotrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zpotri "F_FUNC(zpotri,ZPOTRI)"(char *uplo, int *n, z *a, int *lda, int *info) nogil
cdef zpotri_t *zpotri_f = &_fortran_zpotri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zpotrs "F_FUNC(zpotrs,ZPOTRS)"(char *uplo, int *n, int *nrhs, z *a, int *lda, z *b, int *ldb, int *info) nogil
cdef zpotrs_t *zpotrs_f = &_fortran_zpotrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zppcon "F_FUNC(zppcon,ZPPCON)"(char *uplo, int *n, z *ap, d *anorm, d *rcond, z *work, d *rwork, int *info) nogil
cdef zppcon_t *zppcon_f = &_fortran_zppcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zppequ "F_FUNC(zppequ,ZPPEQU)"(char *uplo, int *n, z *ap, d *s, d *scond, d *amax, int *info) nogil
cdef zppequ_t *zppequ_f = &_fortran_zppequ

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zpprfs "F_FUNC(zpprfs,ZPPRFS)"(char *uplo, int *n, int *nrhs, z *ap, z *afp, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zpprfs_t *zpprfs_f = &_fortran_zpprfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zppsv "F_FUNC(zppsv,ZPPSV)"(char *uplo, int *n, int *nrhs, z *ap, z *b, int *ldb, int *info) nogil
cdef zppsv_t *zppsv_f = &_fortran_zppsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zppsvx "F_FUNC(zppsvx,ZPPSVX)"(char *fact, char *uplo, int *n, int *nrhs, z *ap, z *afp, char *equed, d *s, z *b, int *ldb, z *x, int *ldx, d *rcond, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zppsvx_t *zppsvx_f = &_fortran_zppsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zpptrf "F_FUNC(zpptrf,ZPPTRF)"(char *uplo, int *n, z *ap, int *info) nogil
cdef zpptrf_t *zpptrf_f = &_fortran_zpptrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zpptri "F_FUNC(zpptri,ZPPTRI)"(char *uplo, int *n, z *ap, int *info) nogil
cdef zpptri_t *zpptri_f = &_fortran_zpptri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zpptrs "F_FUNC(zpptrs,ZPPTRS)"(char *uplo, int *n, int *nrhs, z *ap, z *b, int *ldb, int *info) nogil
cdef zpptrs_t *zpptrs_f = &_fortran_zpptrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zptcon "F_FUNC(zptcon,ZPTCON)"(int *n, d *d, z *e, d *anorm, d *rcond, d *rwork, int *info) nogil
cdef zptcon_t *zptcon_f = &_fortran_zptcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zpteqr "F_FUNC(zpteqr,ZPTEQR)"(char *compz, int *n, d *d, d *e, z *z, int *ldz, d *work, int *info) nogil
cdef zpteqr_t *zpteqr_f = &_fortran_zpteqr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zptrfs "F_FUNC(zptrfs,ZPTRFS)"(char *uplo, int *n, int *nrhs, d *d, z *e, d *df, z *ef, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zptrfs_t *zptrfs_f = &_fortran_zptrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zptsv "F_FUNC(zptsv,ZPTSV)"(int *n, int *nrhs, d *d, z *e, z *b, int *ldb, int *info) nogil
cdef zptsv_t *zptsv_f = &_fortran_zptsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zptsvx "F_FUNC(zptsvx,ZPTSVX)"(char *fact, int *n, int *nrhs, d *d, z *e, d *df, z *ef, z *b, int *ldb, z *x, int *ldx, d *rcond, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zptsvx_t *zptsvx_f = &_fortran_zptsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zpttrf "F_FUNC(zpttrf,ZPTTRF)"(int *n, d *d, z *e, int *info) nogil
cdef zpttrf_t *zpttrf_f = &_fortran_zpttrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zpttrs "F_FUNC(zpttrs,ZPTTRS)"(char *uplo, int *n, int *nrhs, d *d, z *e, z *b, int *ldb, int *info) nogil
cdef zpttrs_t *zpttrs_f = &_fortran_zpttrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zptts2 "F_FUNC(zptts2,ZPTTS2)"(int *iuplo, int *n, int *nrhs, d *d, z *e, z *b, int *ldb) nogil
cdef zptts2_t *zptts2_f = &_fortran_zptts2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zrot "F_FUNC(zrot,ZROT)"(int *n, z *cx, int *incx, z *cy, int *incy, d *c, z *s) nogil
cdef zrot_t *zrot_f = &_fortran_zrot

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zspcon "F_FUNC(zspcon,ZSPCON)"(char *uplo, int *n, z *ap, int *ipiv, d *anorm, d *rcond, z *work, int *info) nogil
cdef zspcon_t *zspcon_f = &_fortran_zspcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zspmv "F_FUNC(zspmv,ZSPMV)"(char *uplo, int *n, z *alpha, z *ap, z *x, int *incx, z *beta, z *y, int *incy) nogil
cdef zspmv_t *zspmv_f = &_fortran_zspmv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zspr "F_FUNC(zspr,ZSPR)"(char *uplo, int *n, z *alpha, z *x, int *incx, z *ap) nogil
cdef zspr_t *zspr_f = &_fortran_zspr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zsprfs "F_FUNC(zsprfs,ZSPRFS)"(char *uplo, int *n, int *nrhs, z *ap, z *afp, int *ipiv, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zsprfs_t *zsprfs_f = &_fortran_zsprfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zspsv "F_FUNC(zspsv,ZSPSV)"(char *uplo, int *n, int *nrhs, z *ap, int *ipiv, z *b, int *ldb, int *info) nogil
cdef zspsv_t *zspsv_f = &_fortran_zspsv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zspsvx "F_FUNC(zspsvx,ZSPSVX)"(char *fact, char *uplo, int *n, int *nrhs, z *ap, z *afp, int *ipiv, z *b, int *ldb, z *x, int *ldx, d *rcond, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zspsvx_t *zspsvx_f = &_fortran_zspsvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zsptrf "F_FUNC(zsptrf,ZSPTRF)"(char *uplo, int *n, z *ap, int *ipiv, int *info) nogil
cdef zsptrf_t *zsptrf_f = &_fortran_zsptrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zsptri "F_FUNC(zsptri,ZSPTRI)"(char *uplo, int *n, z *ap, int *ipiv, z *work, int *info) nogil
cdef zsptri_t *zsptri_f = &_fortran_zsptri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zsptrs "F_FUNC(zsptrs,ZSPTRS)"(char *uplo, int *n, int *nrhs, z *ap, int *ipiv, z *b, int *ldb, int *info) nogil
cdef zsptrs_t *zsptrs_f = &_fortran_zsptrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zstedc "F_FUNC(zstedc,ZSTEDC)"(char *compz, int *n, d *d, d *e, z *z, int *ldz, z *work, int *lwork, d *rwork, int *lrwork, int *iwork, int *liwork, int *info) nogil
cdef zstedc_t *zstedc_f = &_fortran_zstedc

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zstegr "F_FUNC(zstegr,ZSTEGR)"(char *jobz, char *range, int *n, d *d, d *e, d *vl, d *vu, int *il, int *iu, d *abstol, int *m, d *w, z *z, int *ldz, int *isuppz, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef zstegr_t *zstegr_f = &_fortran_zstegr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zstein "F_FUNC(zstein,ZSTEIN)"(int *n, d *d, d *e, int *m, d *w, int *iblock, int *isplit, z *z, int *ldz, d *work, int *iwork, int *ifail, int *info) nogil
cdef zstein_t *zstein_f = &_fortran_zstein

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zstemr "F_FUNC(zstemr,ZSTEMR)"(char *jobz, char *range, int *n, d *d, d *e, d *vl, d *vu, int *il, int *iu, int *m, d *w, z *z, int *ldz, int *nzc, int *isuppz, bint *tryrac, d *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef zstemr_t *zstemr_f = &_fortran_zstemr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zsteqr "F_FUNC(zsteqr,ZSTEQR)"(char *compz, int *n, d *d, d *e, z *z, int *ldz, d *work, int *info) nogil
cdef zsteqr_t *zsteqr_f = &_fortran_zsteqr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zsycon "F_FUNC(zsycon,ZSYCON)"(char *uplo, int *n, z *a, int *lda, int *ipiv, d *anorm, d *rcond, z *work, int *info) nogil
cdef zsycon_t *zsycon_f = &_fortran_zsycon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zsymv "F_FUNC(zsymv,ZSYMV)"(char *uplo, int *n, z *alpha, z *a, int *lda, z *x, int *incx, z *beta, z *y, int *incy) nogil
cdef zsymv_t *zsymv_f = &_fortran_zsymv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zsyr "F_FUNC(zsyr,ZSYR)"(char *uplo, int *n, z *alpha, z *x, int *incx, z *a, int *lda) nogil
cdef zsyr_t *zsyr_f = &_fortran_zsyr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zsyrfs "F_FUNC(zsyrfs,ZSYRFS)"(char *uplo, int *n, int *nrhs, z *a, int *lda, z *af, int *ldaf, int *ipiv, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef zsyrfs_t *zsyrfs_f = &_fortran_zsyrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zsysv "F_FUNC(zsysv,ZSYSV)"(char *uplo, int *n, int *nrhs, z *a, int *lda, int *ipiv, z *b, int *ldb, z *work, int *lwork, int *info) nogil
cdef zsysv_t *zsysv_f = &_fortran_zsysv

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zsysvx "F_FUNC(zsysvx,ZSYSVX)"(char *fact, char *uplo, int *n, int *nrhs, z *a, int *lda, z *af, int *ldaf, int *ipiv, z *b, int *ldb, z *x, int *ldx, d *rcond, d *ferr, d *berr, z *work, int *lwork, d *rwork, int *info) nogil
cdef zsysvx_t *zsysvx_f = &_fortran_zsysvx

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zsytf2 "F_FUNC(zsytf2,ZSYTF2)"(char *uplo, int *n, z *a, int *lda, int *ipiv, int *info) nogil
cdef zsytf2_t *zsytf2_f = &_fortran_zsytf2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zsytrf "F_FUNC(zsytrf,ZSYTRF)"(char *uplo, int *n, z *a, int *lda, int *ipiv, z *work, int *lwork, int *info) nogil
cdef zsytrf_t *zsytrf_f = &_fortran_zsytrf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zsytri "F_FUNC(zsytri,ZSYTRI)"(char *uplo, int *n, z *a, int *lda, int *ipiv, z *work, int *info) nogil
cdef zsytri_t *zsytri_f = &_fortran_zsytri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zsytrs "F_FUNC(zsytrs,ZSYTRS)"(char *uplo, int *n, int *nrhs, z *a, int *lda, int *ipiv, z *b, int *ldb, int *info) nogil
cdef zsytrs_t *zsytrs_f = &_fortran_zsytrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztbcon "F_FUNC(ztbcon,ZTBCON)"(char *norm, char *uplo, char *diag, int *n, int *kd, z *ab, int *ldab, d *rcond, z *work, d *rwork, int *info) nogil
cdef ztbcon_t *ztbcon_f = &_fortran_ztbcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztbrfs "F_FUNC(ztbrfs,ZTBRFS)"(char *uplo, char *trans, char *diag, int *n, int *kd, int *nrhs, z *ab, int *ldab, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef ztbrfs_t *ztbrfs_f = &_fortran_ztbrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztbtrs "F_FUNC(ztbtrs,ZTBTRS)"(char *uplo, char *trans, char *diag, int *n, int *kd, int *nrhs, z *ab, int *ldab, z *b, int *ldb, int *info) nogil
cdef ztbtrs_t *ztbtrs_f = &_fortran_ztbtrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztgevc "F_FUNC(ztgevc,ZTGEVC)"(char *side, char *howmny, bint *select, int *n, z *s, int *lds, z *p, int *ldp, z *vl, int *ldvl, z *vr, int *ldvr, int *mm, int *m, z *work, d *rwork, int *info) nogil
cdef ztgevc_t *ztgevc_f = &_fortran_ztgevc

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztgex2 "F_FUNC(ztgex2,ZTGEX2)"(bint *wantq, bint *wantz, int *n, z *a, int *lda, z *b, int *ldb, z *q, int *ldq, z *z, int *ldz, int *j1, int *info) nogil
cdef ztgex2_t *ztgex2_f = &_fortran_ztgex2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztgexc "F_FUNC(ztgexc,ZTGEXC)"(bint *wantq, bint *wantz, int *n, z *a, int *lda, z *b, int *ldb, z *q, int *ldq, z *z, int *ldz, int *ifst, int *ilst, int *info) nogil
cdef ztgexc_t *ztgexc_f = &_fortran_ztgexc

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztgsen "F_FUNC(ztgsen,ZTGSEN)"(int *ijob, bint *wantq, bint *wantz, bint *select, int *n, z *a, int *lda, z *b, int *ldb, z *alpha, z *beta, z *q, int *ldq, z *z, int *ldz, int *m, d *pl, d *pr, d *dif, z *work, int *lwork, int *iwork, int *liwork, int *info) nogil
cdef ztgsen_t *ztgsen_f = &_fortran_ztgsen

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztgsja "F_FUNC(ztgsja,ZTGSJA)"(char *jobu, char *jobv, char *jobq, int *m, int *p, int *n, int *k, int *l, z *a, int *lda, z *b, int *ldb, d *tola, d *tolb, d *alpha, d *beta, z *u, int *ldu, z *v, int *ldv, z *q, int *ldq, z *work, int *ncycle, int *info) nogil
cdef ztgsja_t *ztgsja_f = &_fortran_ztgsja

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztgsna "F_FUNC(ztgsna,ZTGSNA)"(char *job, char *howmny, bint *select, int *n, z *a, int *lda, z *b, int *ldb, z *vl, int *ldvl, z *vr, int *ldvr, d *s, d *dif, int *mm, int *m, z *work, int *lwork, int *iwork, int *info) nogil
cdef ztgsna_t *ztgsna_f = &_fortran_ztgsna

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztgsy2 "F_FUNC(ztgsy2,ZTGSY2)"(char *trans, int *ijob, int *m, int *n, z *a, int *lda, z *b, int *ldb, z *c, int *ldc, z *d, int *ldd, z *e, int *lde, z *f, int *ldf, d *scale, d *rdsum, d *rdscal, int *info) nogil
cdef ztgsy2_t *ztgsy2_f = &_fortran_ztgsy2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztgsyl "F_FUNC(ztgsyl,ZTGSYL)"(char *trans, int *ijob, int *m, int *n, z *a, int *lda, z *b, int *ldb, z *c, int *ldc, z *d, int *ldd, z *e, int *lde, z *f, int *ldf, d *scale, d *dif, z *work, int *lwork, int *iwork, int *info) nogil
cdef ztgsyl_t *ztgsyl_f = &_fortran_ztgsyl

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztpcon "F_FUNC(ztpcon,ZTPCON)"(char *norm, char *uplo, char *diag, int *n, z *ap, d *rcond, z *work, d *rwork, int *info) nogil
cdef ztpcon_t *ztpcon_f = &_fortran_ztpcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztprfs "F_FUNC(ztprfs,ZTPRFS)"(char *uplo, char *trans, char *diag, int *n, int *nrhs, z *ap, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef ztprfs_t *ztprfs_f = &_fortran_ztprfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztptri "F_FUNC(ztptri,ZTPTRI)"(char *uplo, char *diag, int *n, z *ap, int *info) nogil
cdef ztptri_t *ztptri_f = &_fortran_ztptri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztptrs "F_FUNC(ztptrs,ZTPTRS)"(char *uplo, char *trans, char *diag, int *n, int *nrhs, z *ap, z *b, int *ldb, int *info) nogil
cdef ztptrs_t *ztptrs_f = &_fortran_ztptrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztrcon "F_FUNC(ztrcon,ZTRCON)"(char *norm, char *uplo, char *diag, int *n, z *a, int *lda, d *rcond, z *work, d *rwork, int *info) nogil
cdef ztrcon_t *ztrcon_f = &_fortran_ztrcon

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztrevc "F_FUNC(ztrevc,ZTREVC)"(char *side, char *howmny, bint *select, int *n, z *t, int *ldt, z *vl, int *ldvl, z *vr, int *ldvr, int *mm, int *m, z *work, d *rwork, int *info) nogil
cdef ztrevc_t *ztrevc_f = &_fortran_ztrevc

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztrexc "F_FUNC(ztrexc,ZTREXC)"(char *compq, int *n, z *t, int *ldt, z *q, int *ldq, int *ifst, int *ilst, int *info) nogil
cdef ztrexc_t *ztrexc_f = &_fortran_ztrexc

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztrrfs "F_FUNC(ztrrfs,ZTRRFS)"(char *uplo, char *trans, char *diag, int *n, int *nrhs, z *a, int *lda, z *b, int *ldb, z *x, int *ldx, d *ferr, d *berr, z *work, d *rwork, int *info) nogil
cdef ztrrfs_t *ztrrfs_f = &_fortran_ztrrfs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztrsen "F_FUNC(ztrsen,ZTRSEN)"(char *job, char *compq, bint *select, int *n, z *t, int *ldt, z *q, int *ldq, z *w, int *m, d *s, d *sep, z *work, int *lwork, int *info) nogil
cdef ztrsen_t *ztrsen_f = &_fortran_ztrsen

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztrsna "F_FUNC(ztrsna,ZTRSNA)"(char *job, char *howmny, bint *select, int *n, z *t, int *ldt, z *vl, int *ldvl, z *vr, int *ldvr, d *s, d *sep, int *mm, int *m, z *work, int *ldwork, d *rwork, int *info) nogil
cdef ztrsna_t *ztrsna_f = &_fortran_ztrsna

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztrsyl "F_FUNC(ztrsyl,ZTRSYL)"(char *trana, char *tranb, int *isgn, int *m, int *n, z *a, int *lda, z *b, int *ldb, z *c, int *ldc, d *scale, int *info) nogil
cdef ztrsyl_t *ztrsyl_f = &_fortran_ztrsyl

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztrti2 "F_FUNC(ztrti2,ZTRTI2)"(char *uplo, char *diag, int *n, z *a, int *lda, int *info) nogil
cdef ztrti2_t *ztrti2_f = &_fortran_ztrti2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztrtri "F_FUNC(ztrtri,ZTRTRI)"(char *uplo, char *diag, int *n, z *a, int *lda, int *info) nogil
cdef ztrtri_t *ztrtri_f = &_fortran_ztrtri

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztrtrs "F_FUNC(ztrtrs,ZTRTRS)"(char *uplo, char *trans, char *diag, int *n, int *nrhs, z *a, int *lda, z *b, int *ldb, int *info) nogil
cdef ztrtrs_t *ztrtrs_f = &_fortran_ztrtrs

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztzrqf "F_FUNC(ztzrqf,ZTZRQF)"(int *m, int *n, z *a, int *lda, z *tau, int *info) nogil
cdef ztzrqf_t *ztzrqf_f = &_fortran_ztzrqf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_ztzrzf "F_FUNC(ztzrzf,ZTZRZF)"(int *m, int *n, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef ztzrzf_t *ztzrzf_f = &_fortran_ztzrzf

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zung2l "F_FUNC(zung2l,ZUNG2L)"(int *m, int *n, int *k, z *a, int *lda, z *tau, z *work, int *info) nogil
cdef zung2l_t *zung2l_f = &_fortran_zung2l

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zung2r "F_FUNC(zung2r,ZUNG2R)"(int *m, int *n, int *k, z *a, int *lda, z *tau, z *work, int *info) nogil
cdef zung2r_t *zung2r_f = &_fortran_zung2r

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zungbr "F_FUNC(zungbr,ZUNGBR)"(char *vect, int *m, int *n, int *k, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef zungbr_t *zungbr_f = &_fortran_zungbr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zunghr "F_FUNC(zunghr,ZUNGHR)"(int *n, int *ilo, int *ihi, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef zunghr_t *zunghr_f = &_fortran_zunghr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zungl2 "F_FUNC(zungl2,ZUNGL2)"(int *m, int *n, int *k, z *a, int *lda, z *tau, z *work, int *info) nogil
cdef zungl2_t *zungl2_f = &_fortran_zungl2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zunglq "F_FUNC(zunglq,ZUNGLQ)"(int *m, int *n, int *k, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef zunglq_t *zunglq_f = &_fortran_zunglq

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zungql "F_FUNC(zungql,ZUNGQL)"(int *m, int *n, int *k, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef zungql_t *zungql_f = &_fortran_zungql

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zungqr "F_FUNC(zungqr,ZUNGQR)"(int *m, int *n, int *k, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef zungqr_t *zungqr_f = &_fortran_zungqr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zungr2 "F_FUNC(zungr2,ZUNGR2)"(int *m, int *n, int *k, z *a, int *lda, z *tau, z *work, int *info) nogil
cdef zungr2_t *zungr2_f = &_fortran_zungr2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zungrq "F_FUNC(zungrq,ZUNGRQ)"(int *m, int *n, int *k, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef zungrq_t *zungrq_f = &_fortran_zungrq

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zungtr "F_FUNC(zungtr,ZUNGTR)"(char *uplo, int *n, z *a, int *lda, z *tau, z *work, int *lwork, int *info) nogil
cdef zungtr_t *zungtr_f = &_fortran_zungtr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zunm2l "F_FUNC(zunm2l,ZUNM2L)"(char *side, char *trans, int *m, int *n, int *k, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *info) nogil
cdef zunm2l_t *zunm2l_f = &_fortran_zunm2l

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zunm2r "F_FUNC(zunm2r,ZUNM2R)"(char *side, char *trans, int *m, int *n, int *k, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *info) nogil
cdef zunm2r_t *zunm2r_f = &_fortran_zunm2r

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zunmbr "F_FUNC(zunmbr,ZUNMBR)"(char *vect, char *side, char *trans, int *m, int *n, int *k, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *lwork, int *info) nogil
cdef zunmbr_t *zunmbr_f = &_fortran_zunmbr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zunmhr "F_FUNC(zunmhr,ZUNMHR)"(char *side, char *trans, int *m, int *n, int *ilo, int *ihi, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *lwork, int *info) nogil
cdef zunmhr_t *zunmhr_f = &_fortran_zunmhr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zunml2 "F_FUNC(zunml2,ZUNML2)"(char *side, char *trans, int *m, int *n, int *k, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *info) nogil
cdef zunml2_t *zunml2_f = &_fortran_zunml2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zunmlq "F_FUNC(zunmlq,ZUNMLQ)"(char *side, char *trans, int *m, int *n, int *k, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *lwork, int *info) nogil
cdef zunmlq_t *zunmlq_f = &_fortran_zunmlq

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zunmql "F_FUNC(zunmql,ZUNMQL)"(char *side, char *trans, int *m, int *n, int *k, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *lwork, int *info) nogil
cdef zunmql_t *zunmql_f = &_fortran_zunmql

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zunmqr "F_FUNC(zunmqr,ZUNMQR)"(char *side, char *trans, int *m, int *n, int *k, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *lwork, int *info) nogil
cdef zunmqr_t *zunmqr_f = &_fortran_zunmqr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zunmr2 "F_FUNC(zunmr2,ZUNMR2)"(char *side, char *trans, int *m, int *n, int *k, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *info) nogil
cdef zunmr2_t *zunmr2_f = &_fortran_zunmr2

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zunmr3 "F_FUNC(zunmr3,ZUNMR3)"(char *side, char *trans, int *m, int *n, int *k, int *l, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *info) nogil
cdef zunmr3_t *zunmr3_f = &_fortran_zunmr3

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zunmrq "F_FUNC(zunmrq,ZUNMRQ)"(char *side, char *trans, int *m, int *n, int *k, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *lwork, int *info) nogil
cdef zunmrq_t *zunmrq_f = &_fortran_zunmrq

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zunmrz "F_FUNC(zunmrz,ZUNMRZ)"(char *side, char *trans, int *m, int *n, int *k, int *l, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *lwork, int *info) nogil
cdef zunmrz_t *zunmrz_f = &_fortran_zunmrz

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zunmtr "F_FUNC(zunmtr,ZUNMTR)"(char *side, char *uplo, char *trans, int *m, int *n, z *a, int *lda, z *tau, z *c, int *ldc, z *work, int *lwork, int *info) nogil
cdef zunmtr_t *zunmtr_f = &_fortran_zunmtr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zupgtr "F_FUNC(zupgtr,ZUPGTR)"(char *uplo, int *n, z *ap, z *tau, z *q, int *ldq, z *work, int *info) nogil
cdef zupgtr_t *zupgtr_f = &_fortran_zupgtr

cdef extern from "_lapack_subroutine_wrappers.h":
    void _fortran_zupmtr "F_FUNC(zupmtr,ZUPMTR)"(char *side, char *uplo, char *trans, int *m, int *n, z *ap, z *tau, z *c, int *ldc, z *work, int *info) nogil
cdef zupmtr_t *zupmtr_f = &_fortran_zupmtr


# Python accessible wrappers for testing:

def _test_dlamch(cmach):
    # This conversion is necessary to handle Python 3 strings.
    cmach_bytes = bytes(cmach)
    # Now that it is a bytes representation, a non-temporary variable
    # must be passed as a part of the function call.
    cdef char* cmach_char = cmach_bytes
    return dlamch_f(cmach_char)

def _test_slamch(cmach):
    # This conversion is necessary to handle Python 3 strings.
    cmach_bytes = bytes(cmach)
    # Now that it is a bytes representation, a non-temporary variable
    # must be passed as a part of the function call.
    cdef char* cmach_char = cmach_bytes
    return slamch_f(cmach_char)
