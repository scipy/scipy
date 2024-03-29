from libcpp.complex cimport complex as ccomplex

cimport numpy as cnp
cnp.import_array()

cdef extern from "special/airy.h" nogil:
    void specfun_airyzo 'special::airyzo'(int nt, int kf, double *xa, double *xb, double *xc, double *xd)

cdef extern from "special/fresnel.h" nogil:
    void specfun_fcszo 'special::fcszo'(int kf, int nt, ccomplex[double] *zo)

cdef extern from "special/kelvin.h" nogil:
    void specfun_klvnzo 'special::klvnzo'(int nt, int kd, double *zo)

cdef extern from "special/par_cyl.h" nogil:
    void specfun_pbdv 'special::detail::pbdv'(double x, double v, double *dv, double *dp, double *pdf, double *pdd)
    void specfun_pbvv 'special::detail::pbvv'(double x, double v, double *vv, double *vp, double *pvf, double *pvd)

cdef extern from "special/specfun/specfun.h" nogil:
    void specfun_bernob 'special::specfun::bernob'(int n, double *bn)
    void specfun_cerzo 'special::specfun::cerzo'(int nt, ccomplex[double] *zo)
    void specfun_clpmn 'special::specfun::clpmn'(ccomplex[double] z, int m, int n, int ntype, ccomplex[double] *cpm, ccomplex[double] *cpd)
    void specfun_clpn 'special::specfun::clpn'(int n, ccomplex[double] z, ccomplex[double] *cpn, ccomplex[double] *cpd)
    void specfun_clqmn 'special::specfun::clqmn'(ccomplex[double] z, int m, int n, ccomplex[double] *cqm, ccomplex[double] *cqd)
    void specfun_clqn 'special::specfun::clqn'(int n, ccomplex[double] z, ccomplex[double] *cqn, ccomplex[double] *cqd)
    void specfun_cpbdn 'special::specfun::cpbdn'(int n, ccomplex[double] z, ccomplex[double] *cpb, ccomplex[double] *cpd)
    void specfun_cyzo 'special::specfun::cyzo'(int nt, int kf, int kc, ccomplex[double] *zo, ccomplex[double] *zv)
    void specfun_eulerb 'special::specfun::eulerb'(int n, double *en)
    void specfun_fcoef 'special::specfun::fcoef'(int kd, int m, double q, double a, double *fc)
    void specfun_jdzo 'special::specfun::jdzo'(int nt, double *zo, int *n, int *m, int *p)
    void specfun_jyzo 'special::specfun::jyzo'(int n, int nt, double *rj0, double *rj1, double *ry0, double *ry1)
    void specfun_lamn 'special::specfun::lamn'(int n, double x, int *nm, double *bl, double *dl)
    void specfun_lamv 'special::specfun::lamv'(double v, double x, double *vm, double *vl, double *dl)
    void specfun_lpmn 'special::specfun::lpmn'(int m, int n, double x, double *pm, double *pd)
    void specfun_lpn 'special::specfun::lpn'(int n, double x, double *pn, double *pd)
    void specfun_lqmn 'special::specfun::lqmn'(double x, int m, int n, double *qm, double *qd)
    void specfun_lqnb 'special::specfun::lqnb'(int n, double x, double* qn, double* qd)
    void specfun_rctj 'special::specfun::rctj'(int n, double x, int *nm, double *rj, double *dj)
    void specfun_rcty 'special::specfun::rcty'(int n, double x, int *nm, double *ry, double *dy)
    void specfun_sdmn 'special::specfun::sdmn'(int m, int n, double c, double cv, double kd, double *df)
    void specfun_segv 'special::specfun::segv'(int m, int n, double c, int kd, double *cv, double *eg)


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


def clpmn(int m, int n, ccomplex[double] z, int ntype):
    """
    Compute the associated Legendre functions Pmn(z) and their derivatives
    Pmn'(z) for a complex argument. This is a wrapper for the function
    'specfun_clpmn'.
    """
    cdef ccomplex[double] *ccpm
    cdef ccomplex[double] *ccpd
    cdef cnp.npy_intp dims[2]
    dims[0] = m+1
    dims[1] = n+1

    # specfun_clpmn initializes the array internally
    cpm = cnp.PyArray_SimpleNew(2, dims, cnp.NPY_COMPLEX128)
    cpd = cnp.PyArray_SimpleNew(2, dims, cnp.NPY_COMPLEX128)
    ccpm = <ccomplex[double] *>cnp.PyArray_DATA(cpm)
    ccpd = <ccomplex[double] *>cnp.PyArray_DATA(cpd)
    specfun_clpmn(z, m, n, ntype, ccpm, ccpd)
    return cpm, cpd


def clpn(int n1, ccomplex[double] z):
    """
    Compute Legendre polynomials Pn(z) and their derivatives Pn'(z) for
    a complex argument. This is a wrapper for the function 'specfun_clpn'.
    """
    cdef ccomplex[double] *ccpn
    cdef ccomplex[double] *ccpd
    cdef cnp.npy_intp dims[1]
    dims[0] = n1 + 1

    # specfun_clpn initializes the array internally
    cpn = cnp.PyArray_SimpleNew(1, dims, cnp.NPY_COMPLEX128)
    cpd = cnp.PyArray_SimpleNew(1, dims, cnp.NPY_COMPLEX128)
    ccpn = <ccomplex[double] *>cnp.PyArray_DATA(cpn)
    ccpd = <ccomplex[double] *>cnp.PyArray_DATA(cpd)
    specfun_clpn(n1, z, ccpn, ccpd)
    return cpn, cpd


def clqmn(int m, int n, ccomplex[double] z):
    """
    Compute the associated Legendre functions of the second kind,
    Qmn(z) and Qmn'(z), for a complex argument. This is a wrapper
    for the function 'specfun_clqmn'.
    """
    cdef ccomplex[double] *ccqm
    cdef ccomplex[double] *ccqd
    cdef cnp.npy_intp dims[2]
    dims[0] = m+1
    dims[1] = n+1

    # specfun_clqmn initializes the array internally
    cqm = cnp.PyArray_SimpleNew(2, dims, cnp.NPY_COMPLEX128)
    cqd = cnp.PyArray_SimpleNew(2, dims, cnp.NPY_COMPLEX128)
    ccqm = <ccomplex[double] *>cnp.PyArray_DATA(cqm)
    ccqd = <ccomplex[double] *>cnp.PyArray_DATA(cqd)
    specfun_clqmn(z, m, n, ccqm, ccqd)
    return cqm, cqd


def clqn(int n, ccomplex[double] z):
    """
    Compute the Legendre functions Qn(z) and their derivatives
    Qn'(z) for a complex argument. This is a wrapper for the
    function 'specfun_clqn'.
    """
    cdef ccomplex[double] *ccqn
    cdef ccomplex[double] *ccqd
    cdef cnp.npy_intp dims[1]
    dims[0] = n + 1

    # specfun_clpn initializes the array internally
    cqn = cnp.PyArray_SimpleNew(1, dims, cnp.NPY_COMPLEX128)
    cqd = cnp.PyArray_SimpleNew(1, dims, cnp.NPY_COMPLEX128)
    ccqn = <ccomplex[double] *>cnp.PyArray_DATA(cqn)
    ccqd = <ccomplex[double] *>cnp.PyArray_DATA(cqd)
    specfun_clqn(n, <ccomplex[double]> z, ccqn, ccqd)
    return cqn, cqd


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

    specfun_jdzo(nt, zzo, nn, mm, pp)

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


def lpmn(int m, int n, double x):
    """
    Compute the associated Legendre functions Pmn(x) and their
    derivatives Pmn'(x) for real argument. This is a wrapper for
    the function 'specfun_lpmn'.
    """
    cdef double *cpm
    cdef double *cpd
    cdef cnp.npy_intp dims[2]
    dims[0] = m+1
    dims[1] = n+1

    # specfun_clpmn initializes the array internally
    pm = cnp.PyArray_SimpleNew(2, dims, cnp.NPY_FLOAT64)
    pd = cnp.PyArray_SimpleNew(2, dims, cnp.NPY_FLOAT64)
    cpm = <cnp.float64_t *>cnp.PyArray_DATA(pm)
    cpd = <cnp.float64_t *>cnp.PyArray_DATA(pd)
    specfun_lpmn(m, n, x, cpm, cpd)
    return pm, pd


def lpn(int n, double z):
    """
    Compute Legendre polynomials Pn(x) and their derivatives
    Pn'(x). This is a wrapper for the function 'specfun_lpn'.
    """
    cdef double *ppn
    cdef double *ppd
    cdef cnp.npy_intp dims[1]
    dims[0] = n + 1

    pn = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    pd = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    ppn = <cnp.float64_t *>cnp.PyArray_DATA(pn)
    ppd = <cnp.float64_t *>cnp.PyArray_DATA(pd)
    specfun_lpn(n, z, ppn, ppd)
    return pn, pd


def lqmn(int m, int n, double x):
    """
    Purpose: Compute the associated Legendre functions of the
    second kind, Qmn(x) and Qmn'(x). This is a wrapper for
    the function 'specfun_lqmn'.
    """
    cdef double *cqm
    cdef double *cqd
    cdef cnp.npy_intp dims[2]
    dims[0] = m+1
    dims[1] = n+1

    # specfun_clpmn initializes the array internally
    qm = cnp.PyArray_SimpleNew(2, dims, cnp.NPY_FLOAT64)
    qd = cnp.PyArray_SimpleNew(2, dims, cnp.NPY_FLOAT64)
    cqm = <cnp.float64_t *>cnp.PyArray_DATA(qm)
    cqd = <cnp.float64_t *>cnp.PyArray_DATA(qd)
    specfun_lqmn(x, m, n, cqm, cqd)
    return qm, qd


def lqnb(int n, double x):
    """
    Compute Legendre functions Qn(x) & Qn'(x). This is a wrapper for
    the function 'specfun_lqnb'.
    """
    cdef double *qqn
    cdef double *qqd
    cdef cnp.npy_intp dims[1]
    dims[0] = n + 1

    qn = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    qd = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    qqn = <cnp.float64_t *>cnp.PyArray_DATA(qn)
    qqd = <cnp.float64_t *>cnp.PyArray_DATA(qd)
    specfun_lqnb(n, x, qqn, qqd)
    return qn, qd


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

    
def rctj(int n, double x):
    """
    Compute Riccati-Bessel functions of the first kind and their
    derivatives. This is a wrapper for the function 'specfun_rctj'.
    """
    cdef int nm
    cdef double *rrj
    cdef double *ddj
    cdef cnp.npy_intp dims[1]
    dims[0] = n + 1

    rj = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    dj = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    rrj = <cnp.float64_t *>cnp.PyArray_DATA(rj)
    ddj = <cnp.float64_t *>cnp.PyArray_DATA(dj)
    specfun_rctj(n, x, &nm, rrj, ddj)
    return nm, rj, dj


def rcty(int n, double x):
    """
    Compute Riccati-Bessel functions of the second kind and their
    derivatives This is a wrapper for the function 'specfun_rcty'.
    """
    cdef int nm
    cdef double *rry
    cdef double *ddy
    cdef cnp.npy_intp dims[1]
    dims[0] = n + 1

    ry = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    dy = cnp.PyArray_ZEROS(1, dims, cnp.NPY_FLOAT64, 0)
    rry = <cnp.float64_t *>cnp.PyArray_DATA(ry)
    ddy = <cnp.float64_t *>cnp.PyArray_DATA(dy)
    specfun_rcty(n, x, &nm, rry, ddy)
    return nm, ry, dy


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
    specfun_sdmn(m, n, c, cv, kd, ddf)
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
    specfun_segv(m, n, c, kd, &cv, eeg)
    return cv, eg
