cdef extern from "gmp.h":
    # gmp integers
    ctypedef long mp_limb_t

    ctypedef struct __mpz_struct:
        int _mp_alloc
        int _mp_size
        mp_limb_t* _mp_d

    ctypedef __mpz_struct mpz_t[1]
    ctypedef __mpz_struct *mpz_ptr
    ctypedef const __mpz_struct *mpz_srcptr

    # gmp rationals
    ctypedef struct __mpq_struct:
        __mpz_struct _mp_num
        __mpz_struct _mp_den

    ctypedef __mpq_struct mpq_t[1]
    ctypedef __mpq_struct *mpq_ptr
    ctypedef const __mpq_struct *mpq_srcptr

    void mpz_set(mpz_t rop, mpz_t op)
    void mpq_set(mpq_ptr rop, mpq_srcptr op)
    void mpq_set_num(mpq_t rational, mpz_t numerator)
    void mpq_set_den(mpq_t rational, mpz_t denominator)


cdef extern from "mpfr.h":
    # mpfr reals
    ctypedef int mpfr_sign_t
    ctypedef long mpfr_prec_t
    ctypedef long mpfr_exp_t

    ctypedef struct __mpfr_struct:
        mpfr_prec_t _mpfr_prec
        mpfr_sign_t _mpfr_sign
        mpfr_exp_t _mpfr_exp
        mp_limb_t* _mpfr_d

    ctypedef __mpfr_struct mpfr_t[1]
    ctypedef __mpfr_struct *mpfr_ptr
    ctypedef const __mpfr_struct *mpfr_srcptr

    ctypedef enum mpfr_rnd_t:
        MPFR_RNDN
        MPFR_RNDZ
        MPFR_RNDU
        MPFR_RNDD
        MPFR_RNDA
        MPFR_RNDF
        MPFR_RNDNA

    mpfr_prec_t mpfr_get_prec(mpfr_t x)
    int mpfr_set(mpfr_t rop, mpfr_t op, mpfr_rnd_t rnd)


cdef extern from "mpc.h":
    # mpc complexes
    ctypedef struct __mpc_struct:
        mpfr_t re
        mpfr_t im

    ctypedef __mpc_struct mpc_t[1];
    ctypedef __mpc_struct *mpc_ptr;
    ctypedef const __mpc_struct *mpc_srcptr;
    ctypedef enum mpc_rnd_t:
        MPC_RNDNN
        MPC_RNDNZ
        MPC_RNDNU
        MPC_RNDND
        MPC_RNDZN
        MPC_RNDZZ
        MPC_RNDZU
        MPC_RNDZD
        MPC_RNDUN
        MPC_RNDUZ
        MPC_RNDUU
        MPC_RNDUD
        MPC_RNDDN
        MPC_RNDDZ
        MPC_RNDDU
        MPC_RNDDD

    mpfr_prec_t mpc_get_prec(mpc_srcptr x)
    void mpc_get_prec2(mpfr_prec_t *pr, mpfr_prec_t *pi, mpc_srcptr x)
    int mpc_set(mpc_ptr rop, mpc_srcptr op, mpc_rnd_t rnd)
    int mpc_set_fr_fr(mpc_ptr rop, mpfr_srcptr rp, mpfr_srcptr ip, mpc_rnd_t rnd)


cdef extern from "gmpy2.h":
    # Initialize the C-API
    # This must be called before any other functions, but not to access
    # the types.
    cdef int import_gmpy2() except -1

    # Object types
    ctypedef class gmpy2.mpz [object MPZ_Object]:
        cdef mpz_t z
    ctypedef class gmpy2.mpq [object MPQ_Object]:
        cdef mpq_t q
    ctypedef class gmpy2.mpfr [object MPFR_Object]:
        cdef mpfr_t f
        cdef int rc
    ctypedef class gmpy2.mpc [object MPC_Object]:
        cdef mpc_t c
        cdef int rc

    # Object creation
    cdef mpz GMPy_MPZ_New(void *)
    cdef mpq GMPy_MPQ_New(void *)
    cdef mpfr GMPy_MPFR_New(mpfr_prec_t prec, void *)
    cdef mpc GMPy_MPC_New(mpfr_prec_t rprec, mpfr_prec_t iprec, void *)

    # C field access
    cdef mpz_t MPZ(mpz)
    cdef mpq_t MPQ(mpq)
    cdef mpfr_t MPFR(mpfr)
    cdef mpc_t MPC(mpc)

    # Type check
    cdef bint MPZ_Check(object)
    cdef bint MPQ_Check(object)
    cdef bint MPFR_Check(object)
    cdef bint MPC_Check(object)


# Build a gmpy2 mpz from a gmp mpz
cdef inline mpz GMPy_MPZ_From_mpz(mpz_srcptr z):
    cdef mpz res = GMPy_MPZ_New(NULL)
    mpz_set(res.z, z)
    return res

# Build a gmpy2 mpq from a gmp mpq
cdef inline mpq GMPy_MPQ_From_mpq(mpq_srcptr q):
    cdef mpq res = GMPy_MPQ_New(NULL)
    mpq_set(res.q, q)
    return res

# Build a gmpy2 mpq from gmp mpz numerator and denominator
cdef inline mpq GMPy_MPQ_From_mpz(mpz_srcptr num, mpz_srcptr den):
    cdef mpq res = GMPy_MPQ_New(NULL)
    mpq_set_num(res.q, num)
    mpq_set_den(res.q, den)
    return res

# Build a gmpy2 mpfr from a mpfr
cdef inline mpfr GMPy_MPFR_From_mpfr(mpfr_srcptr x):
    cdef mpfr res = GMPy_MPFR_New(mpfr_get_prec(x), NULL)
    mpfr_set(res.f, x, MPFR_RNDN)
    return res

# Build a gmpy2 mpc from a mpc
cdef inline mpc GMPy_MPC_From_mpc(mpc_srcptr c):
    cdef mpfr_prec_t pr
    cdef mpfr_prec_t pi
    mpc_get_prec2(&pr, &pi, c)
    cdef mpc res = GMPy_MPC_New(pr, pi, NULL)
    mpc_set(res.c, c, MPC_RNDNN)
    return res

# Build a gmpy2 mpc from a real part mpfr and an imaginary part mpfr
cdef inline mpc GMPy_MPC_From_mpfr(mpfr_srcptr re, mpfr_srcptr im):
    cdef mpc res = GMPy_MPC_New(mpfr_get_prec(re), mpfr_get_prec(im), NULL)
    # We intentionally use MPFR funtions instead of MPC functions here
    # in order not to add an unneeded dependency on MPC. It's probably
    # faster too this way.
    mpfr_set(res.c.re, re, MPFR_RNDN)
    mpfr_set(res.c.im, im, MPFR_RNDN)
    return res
