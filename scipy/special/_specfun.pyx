from libcpp.complex cimport complex as ccomplex

cimport numpy as cnp
cnp.import_array()

cdef extern from "xsf/airy.h" nogil:
    void specfun_airyzo 'xsf::airyzo'(int nt, int kf, double *xa, double *xb, double *xc, double *xd)

cdef extern from "xsf/fresnel.h" nogil:
    void specfun_fcszo 'xsf::fcszo'(int kf, int nt, ccomplex[double] *zo)

cdef extern from "xsf/kelvin.h" nogil:
    void specfun_klvnzo 'xsf::klvnzo'(int nt, int kd, double *zo)

cdef extern from "xsf/par_cyl.h" nogil:
    void specfun_pbdv 'xsf::detail::pbdv'(double x, double v, double *dv, double *dp, double *pdf, double *pdd)
    void specfun_pbvv 'xsf::detail::pbvv'(double x, double v, double *vv, double *vp, double *pvf, double *pvd)

cdef extern from "xsf/specfun/specfun.h" namespace "xsf::specfun":

    cpdef enum class Status:
        OK = 0
        NoMemory
        Other

cdef extern from "xsf/specfun/specfun.h" nogil:
    void specfun_bernob 'xsf::specfun::bernob'(int n, double *bn)
    void specfun_cerzo 'xsf::specfun::cerzo'(int nt, ccomplex[double] *zo)
    void specfun_cpbdn 'xsf::specfun::cpbdn'(int n, ccomplex[double] z, ccomplex[double] *cpb, ccomplex[double] *cpd)
    void specfun_cyzo 'xsf::specfun::cyzo'(int nt, int kf, int kc, ccomplex[double] *zo, ccomplex[double] *zv)
    void specfun_eulerb 'xsf::specfun::eulerb'(int n, double *en)
    void specfun_fcoef 'xsf::specfun::fcoef'(int kd, int m, double q, double a, double *fc)
    Status specfun_jdzo 'xsf::specfun::jdzo'(int nt, double *zo, int *n, int *m, int *p)
    void specfun_jyzo 'xsf::specfun::jyzo'(int n, int nt, double *rj0, double *rj1, double *ry0, double *ry1)
    void specfun_lamn 'xsf::specfun::lamn'(int n, double x, int *nm, double *bl, double *dl)
    void specfun_lamv 'xsf::specfun::lamv'(double v, double x, double *vm, double *vl, double *dl)
    void specfun_lqnb 'xsf::specfun::lqnb'(int n, double x, double* qn, double* qd)
    void specfun_pbdv 'xsf::specfun::pbdv'(double x, double v, double *dv, double *dp, double *pdf, double *pdd)
    void specfun_pbvv 'xsf::specfun::pbvv'(double x, double v, double *vv, double *vp, double *pvf, double *pvd)
    Status specfun_sdmn 'xsf::specfun::sdmn'(int m, int n, double c, double cv, double kd, double *df)
    Status specfun_segv 'xsf::specfun::segv'(int m, int n, double c, int kd, double *cv, double *eg)


def airyzo(int nt, int kf):
    """
    Compute the first NT zeros of Airy functions
    Ai(x) and Ai'(x), a and a', and the associated
    values of Ai(a') and Ai'(a); and the first NT
    zeros of Airy functions Bi(x) and Bi'(x), b and
    b', and the associated values of Bi(b') and
    Bi'(b).

    This is a wrapper for the function 'specfun_airyzo'.
    """
    cdef double *xxa
    cdef double *xxb
    cdef double *xxc
    cdef double *xxd
    cdef cnp.npy_intp dims[1]
    dims[0] = nt

    xa = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    xb = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    xc = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    xd = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)

    xxa = <cnp.float64_t *>cnp.PyArray_DATA(xa)
    xxb = <cnp.float64_t *>cnp.PyArray_DATA(xb)
    xxc = <cnp.float64_t *>cnp.PyArray_DATA(xc)
    xxd = <cnp.float64_t *>cnp.PyArray_DATA(xd)

    specfun_airyzo(nt, kf, xxa, xxb, xxc, xxd)

    return xa, xb, xc, xd


def bernob(int n):
    """
    Compute Bernoulli number Bn for n >= 2. This is a wrapper for the
    function 'specfun_bernob'.
    """
    cdef double *bbn
    cdef cnp.npy_intp dims[1]
    dims[0] = n + 1
    bn = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    bbn = <cnp.float64_t *>cnp.PyArray_DATA(bn)
    specfun_bernob(n, bbn)
    return bn


def cerzo(int nt):
    """
    Evaluate the complex zeros of error function erf(z) using
    the modified Newton's iteration method. This is a wrapper
    for the function 'specfun_cerzo'.
    """
    cdef ccomplex[double] *zzo
    cdef cnp.npy_intp dims[1]
    dims[0] = nt
    zo = cnp.PyArray_ZEROS(1, dims, cnp.NPY_COMPLEX128, 0)
    zzo = <ccomplex[double] *>cnp.PyArray_DATA(zo)
    specfun_cerzo(nt, zzo)
    return zo


def cpbdn(int n, ccomplex[double] z):
    """
    Compute the parabolic cylinder functions Dn(z) and Dn'(z)
    for a complex argument. This is a wrapper for the function
    'specfun_cpbdn'.
    """
    cdef ccomplex[double] *ccpb
    cdef ccomplex[double] *ccpd
    cdef cnp.npy_intp dims[1]
    dims[0] = abs(n) + 2

    cpb = cnp.PyArray_ZEROS(1, dims, cnp.NPY_COMPLEX128, 0)
    cpd = cnp.PyArray_ZEROS(1, dims, cnp.NPY_COMPLEX128, 0)
    ccpb = <ccomplex[double] *>cnp.PyArray_DATA(cpb)
    ccpd = <ccomplex[double] *>cnp.PyArray_DATA(cpd)
    specfun_cpbdn(n, <ccomplex[double]> z, ccpb, ccpd)
    return cpb, cpd


def cyzo(int nt, int kf, int kc):
    """
    Compute the complex zeros of Y0(z), Y1(z) and Y1'(z), and their
    associated values at the zeros using the modified Newton's
    iteration method. This is a wrapper for the function 'specfun_cyzo'.
    """
    cdef ccomplex[double] *zzo
    cdef ccomplex[double] *zzv
    cdef cnp.npy_intp dims[1]
    dims[0] = nt

    zo = cnp.PyArray_ZEROS(1, dims, cnp.NPY_COMPLEX128, 0)
    zv = cnp.PyArray_ZEROS(1, dims, cnp.NPY_COMPLEX128, 0)
    zzo = <ccomplex[double] *>cnp.PyArray_DATA(zo)
    zzv = <ccomplex[double] *>cnp.PyArray_DATA(zv)
    specfun_cyzo(nt, kf, kc, zzo, zzv)
    return zo, zv


def eulerb(int n):
    """
    Compute Bernoulli number Bn for n >= 2. This is a wrapper for the
    function 'specfun_bernob'.
    """
    cdef double *een
    cdef cnp.npy_intp dims[1]
    dims[0] = n + 1
    en = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    een = <cnp.float64_t *>cnp.PyArray_DATA(en)
    specfun_eulerb(n, een)
    return en


def fcoef(int kd, int m, double q, double a):
    """
    Compute expansion coefficients for Mathieu functions and modified
    Mathieu functions
    """
    cdef double *ffc

    cdef cnp.npy_intp dims[1]
    dims[0] = 251

    fc = cnp.PyArray_SimpleNew(1, dims, cnp.NPY_FLOAT64)
    ffc = <cnp.float64_t *>cnp.PyArray_DATA(fc)
    specfun_fcoef(kd, m, q, a, ffc)
    return fc


def fcszo(int kf, int nt):
    """
    Compute the complex zeros of Fresnel integral C(z) or S(z) using
    modified Newton's iteration method. This is a wrapper for the
    function 'specfun_fcszo'.
    """
    cdef ccomplex[double] *zzo
    cdef cnp.npy_intp dims[1]
    dims[0] = nt
    zo = cnp.PyArray_ZEROS(1, dims, cnp.NPY_COMPLEX128, 0)
    zzo = <ccomplex[double] *>cnp.PyArray_DATA(zo)
    specfun_fcszo(kf, nt, zzo)
    return zo


def jdzo(int nt):
    """
    Compute the zeros of Bessel functions Jn(x) and Jn'(x), and
    arrange them in the order of their magnitudes.

    This is a wrapper for the function 'specfun_jdzo'. The
    relationship between nt and the required array sizes is

         nt between     required array size
        -----------------------------------
          0  -  100   ->    (nt + 10)
        100  -  200   ->    (nt + 14)
        200  -  300   ->    (nt + 16)
        300  -  400   ->    (nt + 18)
        400  -  500   ->    (nt + 21)
        500  -  600   ->    (nt + 25)
        600  -  700   ->    (nt + 11)
        700  -  800   ->    (nt +  9)
        800  -  900   ->    (nt +  9)
        900  - 1000   ->    (nt + 10)
        1000 - 1100   ->    (nt + 10)
        1100 - 1200   ->    (nt + 11)

    It can be made a bit more granular but a generic +25 slack seems
    like an easy option instead of costly 1400-long arrays as Fortran
    code did originally, independent from the value of 'nt'.
    """

    cdef double *zzo
    cdef int *nn
    cdef int *mm
    cdef int *pp
    cdef cnp.npy_intp dims[1]

    dims[0] = nt + 25

    zo = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    n = cnp.PyArray_ZEROS(1, dims, cnp.NPY_INT32, 0)
    m = cnp.PyArray_ZEROS(1, dims, cnp.NPY_INT32, 0)
    p = cnp.PyArray_ZEROS(1, dims, cnp.NPY_INT32, 0)
    zzo = <cnp.float64_t *>cnp.PyArray_DATA(zo)
    nn = <int*>cnp.PyArray_DATA(n)
    mm = <int*>cnp.PyArray_DATA(m)
    pp = <int*>cnp.PyArray_DATA(p)

    if (specfun_jdzo(nt, zzo, nn, mm, pp) == Status.NoMemory):
        raise MemoryError('jnjnp_zeros: failed to allocate working memory in jdzo.')

    return n, m, p, zo


def jyzo(int n, int nt):
    """
    Compute the zeros of Bessel functions Jn(x), Yn(x), and their
    derivatives. This is a wrapper for the function 'specfun_jyzo'.
    """
    cdef double *rrj0
    cdef double *rrj1
    cdef double *rry0
    cdef double *rry1
    cdef cnp.npy_intp dims[1]
    dims[0] = nt

    rj0 = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    rj1 = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    ry0 = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    ry1 = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)

    rrj0 = <cnp.float64_t *>cnp.PyArray_DATA(rj0)
    rrj1 = <cnp.float64_t *>cnp.PyArray_DATA(rj1)
    rry0 = <cnp.float64_t *>cnp.PyArray_DATA(ry0)
    rry1 = <cnp.float64_t *>cnp.PyArray_DATA(ry1)

    specfun_jyzo(n, nt, rrj0, rrj1, rry0, rry1)

    return rj0, rj1, ry0, ry1


def klvnzo(int nt, int kd):
    """
    Compute the zeros of Kelvin functions. This is a wrapper for
    the function 'specfun_klvzo'.
    """
    cdef double *zzo
    cdef cnp.npy_intp dims[1]
    dims[0] = nt
    zo = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    zzo = <cnp.float64_t *>cnp.PyArray_DATA(zo)
    specfun_klvnzo(nt, kd, zzo)
    return zo


def lamn(int n, double x):
    """
    Compute lambda functions and their derivatives. This is a wrapper
    for the function 'specfun_lamn'.
    """
    cdef int nm
    cdef double *bbl
    cdef double *ddl
    cdef cnp.npy_intp dims[1]
    dims[0] = n + 1

    bl = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    dl = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    bbl = <cnp.float64_t *>cnp.PyArray_DATA(bl)
    ddl = <cnp.float64_t *>cnp.PyArray_DATA(dl)
    specfun_lamn(n, x, &nm, bbl, ddl)
    return nm, bl, dl


def lamv(double v, double x):
    """
    Compute lambda function with arbitrary order v, and their derivative.
    This is a wrapper for the function 'specfun_lamv'.
    """
    cdef double vm
    cdef double *vvl
    cdef double *ddl
    cdef cnp.npy_intp dims[1]
    dims[0] = int(v) + 1

    vl = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    dl = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    vvl = <cnp.float64_t *>cnp.PyArray_DATA(vl)
    ddl = <cnp.float64_t *>cnp.PyArray_DATA(dl)
    specfun_lamv(v, x, &vm, vvl, ddl)
    return vm, vl, dl


def pbdv(double v, double x):
    cdef double pdf
    cdef double pdd
    cdef double *ddv
    cdef double *ddp
    cdef cnp.npy_intp dims[1]
    dims[0] = abs(<int>v) + 2

    dv = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    dp = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    ddv = <cnp.float64_t *>cnp.PyArray_DATA(dv)
    ddp = <cnp.float64_t *>cnp.PyArray_DATA(dp)
    specfun_pbdv(x, v, ddv, ddp, &pdf, &pdd)
    return dv, dp, pdf, pdd


def pbvv(double v, double x):
    cdef double pvf
    cdef double pvd
    cdef double *dvv
    cdef double *dvp
    cdef cnp.npy_intp dims[1]
    dims[0] = abs(<int>v) + 2

    vv = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    vp = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    dvv = <cnp.float64_t *>cnp.PyArray_DATA(vv)
    dvp = <cnp.float64_t *>cnp.PyArray_DATA(vp)
    specfun_pbvv(x, v, dvv, dvp, &pvf, &pvd)
    return vv, vp, pvf, pvd


def sdmn(int m, int n, double c, double cv, int kd):
    """
    Compute the expansion coefficients of the prolate and oblate
    spheroidal functions, dk. This is a wrapper for the function
    'specfun_sdmn'.
    """
    cdef double *ddf
    cdef int nm = 25 + (int)(0.5 * (n - m) + c);
    cdef cnp.npy_intp dims[1]
    dims[0] = nm + 1

    df = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    ddf = <cnp.float64_t *>cnp.PyArray_DATA(df)
    if specfun_sdmn(m, n, c, cv, kd, ddf) == Status.NoMemory:
        raise MemoryError('sdmn: failed to allocate working memory.')
    return df


def segv(int m, int n, double c, int kd):
    """
    Compute the characteristic values of spheroidal wave functions.
    This is a wrapper for the function 'specfun_segv'.
    """
    cdef double cv
    cdef double *eeg
    cdef cnp.npy_intp dims[1]
    dims[0] = n - m + 1

    eg = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    eeg = <cnp.float64_t *>cnp.PyArray_DATA(eg)
    if specfun_segv(m, n, c, kd, &cv, eeg) == Status.NoMemory:
        # Note: segv is a private function that is called by either pro_cv_seq
        # or obl_cv_seq.  We make the error message useful by including the
        # approriate name in it.
        caller = 'pro_cv_seq' if kd == 1 else 'obl_cv_seq'
        msg = f'{caller}: failed to allocate working memory in segv.'
        raise MemoryError(msg)
    return cv, eg
