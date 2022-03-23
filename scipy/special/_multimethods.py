import functools
import numpy as np
from numpy import dtype, ndarray
from scipy._lib.uarray import Dispatchable, all_of_type, create_multimethod
from scipy.special import _api


__all__ = [
    'ai_zeros',
    'assoc_laguerre',
    'bei_zeros',
    'beip_zeros',
    'ber_zeros',
    'bernoulli',
    'berp_zeros',
    'bi_zeros',
    'chebyc',
    'chebys',
    'chebyt',
    'chebyu',
    'clpmn',
    'comb',
    'diric',
    'ellip_harm',
    'ellip_harm_2',
    'ellip_normal',
    'erf_zeros',
    'euler',
    'factorial',
    'factorial2',
    'factorialk',
    'fresnel_zeros',
    'fresnelc_zeros',
    'fresnels_zeros',
    'gegenbauer',
    'genlaguerre',
    'h1vp',
    'h2vp',
    'hermite',
    'hermitenorm',
    'ivp',
    'jacobi',
    'jn_zeros',
    'jnjnp_zeros',
    'jnp_zeros',
    'jnyn_zeros',
    'jvp',
    'kei_zeros',
    'keip_zeros',
    'kelvin_zeros',
    'ker_zeros',
    'kerp_zeros',
    'kvp',
    'laguerre',
    'lambertw',
    'legendre',
    'lmbda',
    'log_softmax',
    'logsumexp',
    'lpmn',
    'lpn',
    'lqmn',
    'lqn',
    'mathieu_even_coef',
    'mathieu_odd_coef',
    'multigammaln',
    'obl_cv_seq',
    'p_roots',
    'pbdn_seq',
    'pbdv_seq',
    'pbvv_seq',
    'perm',
    'polygamma',
    'pro_cv_seq',
    'riccati_jn',
    'riccati_yn',
    'roots_chebyc',
    'roots_chebys',
    'roots_chebyt',
    'roots_chebyu',
    'roots_gegenbauer',
    'roots_genlaguerre',
    'roots_hermite',
    'roots_hermitenorm',
    'roots_jacobi',
    'roots_laguerre',
    'roots_legendre',
    'roots_sh_chebyt',
    'roots_sh_chebyu',
    'roots_sh_jacobi',
    'roots_sh_legendre',
    'sh_chebyt',
    'sh_chebyu',
    'sh_jacobi',
    'sh_legendre',
    'sinc',
    'softmax',
    'spherical_in',
    'spherical_jn',
    'spherical_kn',
    'spherical_yn',
    'y0_zeros',
    'y1_zeros',
    'y1p_zeros',
    'yn_zeros',
    'ynp_zeros',
    'yvp',
    'zeta'
]

_create_special = functools.partial(
    create_multimethod, domain="numpy.scipy.special"
)


def _get_docs(func):
    """
    Decorator to take the docstring from original
    function and assign to the multimethod.
    """
    func.__doc__ = getattr(_api, func.__name__).__doc__
    return func


def _identity_replacer(args, kwargs, arrays):
    return args, kwargs


def _x_n_kw_k_replacer(args, kwargs, dispatchables):
    def self_method(x, n, k=0.0, *args, **kwargs):
        kw_out = kwargs.copy()
        return (
            dispatchables[0],
            dispatchables[1],
            dispatchables[2],
        ) + args, kw_out

    return self_method(*args, **kwargs)


def _n_replacer(args, kwargs, dispatchables):
    def self_method(n, *args, **kwargs):
        kw_out = kwargs
        return (dispatchables[0],) + args, kw_out

    return self_method(*args, **kwargs)


def _m_n_z_replacer(args, kwargs, dispatchables):
    def self_method(m, n, z, *args, **kwargs):
        kw_out = kwargs
        return (m, dispatchables[0], dispatchables[1]) + args, kw_out

    return self_method(*args, **kwargs)


def _N_k_replacer(args, kwargs, dispatchables):
    def self_method(N, k, *args, **kwargs):
        kw_out = kwargs
        return (dispatchables[0], dispatchables[1]) + args, kw_out

    return self_method(*args, **kwargs)


def _x_n_replacer(args, kwargs, dispatchables):
    def self_method(x, n, *args, **kwargs):
        kw_out = kwargs
        return (dispatchables[0], n) + args, kw_out

    return self_method(*args, **kwargs)


def _n_k_replacer(args, kwargs, dispatchables):
    def self_method(n, k, *args, **kwargs):
        kw_out = kwargs
        return (dispatchables[0], dispatchables[1]) + args, kw_out

    return self_method(*args, **kwargs)


def _v_z_kw_n_replacer(args, kwargs, dispatchables):
    def self_method(v, z, n=1, *args, **kwargs):
        kw_out = kwargs.copy()
        return (
            dispatchables[0],
            dispatchables[1],
            n,
        ) + args, kw_out

    return self_method(*args, **kwargs)


def _z_kw_k_replacer(args, kwargs, dispatchables):
    def self_method(z, k=0, *args, **kwargs):
        kw_out = kwargs.copy()
        return (dispatchables[0], k) + args, kw_out

    return self_method(*args, **kwargs)


def _v_x_replacer(args, kwargs, dispatchables):
    def self_method(v, x, *args, **kwargs):
        kw_out = kwargs
        return (dispatchables[0], dispatchables[1]) + args, kw_out

    return self_method(*args, **kwargs)


def _x_replacer(args, kwargs, dispatchables):
    def self_method(x, *args, **kwargs):
        kw_out = kwargs
        return (dispatchables[0],) + args, kw_out

    return self_method(*args, **kwargs)


def _a_kw_axis_b_replacer(args, kwargs, dispatchables):
    def self_method(a, axis=None, b=None, *args, **kwargs):
        kw_out = kwargs.copy()
        return (dispatchables[0], axis, dispatchables[1]) + args, kw_out

    return self_method(*args, **kwargs)


def _n_z_replacer(args, kwargs, dispatchables):
    def self_method(n, z, *args, **kwargs):
        kw_out = kwargs
        return (dispatchables[0], dispatchables[1]) + args, kw_out

    return self_method(*args, **kwargs)


def _m_q_replacer(args, kwargs, dispatchables):
    def self_method(m, q, *args, **kwargs):
        kw_out = kwargs
        return (m, dispatchables[0]) + args, kw_out

    return self_method(*args, **kwargs)


def _a_replacer(args, kwargs, dispatchables):
    def self_method(a, *args, **kwargs):
        kw_out = kwargs
        return (dispatchables[0],) + args, kw_out

    return self_method(*args, **kwargs)


def _m_n_replacer(args, kwargs, dispatchables):
    def self_method(m, n, *args, **kwargs):
        kw_out = kwargs
        return (m, dispatchables[0]) + args, kw_out

    return self_method(*args, **kwargs)


def _n_x_replacer(args, kwargs, dispatchables):
    def self_method(n, x, *args, **kwargs):
        kw_out = kwargs
        return (dispatchables[0], dispatchables[1]) + args, kw_out

    return self_method(*args, **kwargs)


def _n_p_q_replacer(args, kwargs, dispatchables):
    def self_method(n, p, q, *args, **kwargs):
        kw_out = kwargs
        return (dispatchables[0], p, dispatchables[1]) + args, kw_out

    return self_method(*args, **kwargs)


def _x_kw_q_out_replacer(args, kwargs, dispatchables):
    def self_method(x, q=None, out=None, *args, **kwargs):
        kw_out = kwargs.copy()
        return (
            dispatchables[0],
            dispatchables[1],
            dispatchables[2],
        ) + args, kw_out

    return self_method(*args, **kwargs)


def _onearg_replacer(args, kwargs, dispatchables):
    def self_method(a, out=None, *args, **kwargs):
        kw_out = kwargs.copy()
        if out is not None:
            kw_out["out"] = dispatchables[1:]
        return (dispatchables[0],) + args, kw_out

    return self_method(*args, **kwargs)

def _twoargs_replacer(args, kwargs, dispatchables):
    def self_method(a, b, out=None, *args, **kwargs):
        kw_out = kwargs.copy()
        if out is not None:
            kw_out["out"] = dispatchables[2:]
        return (dispatchables[0], dispatchables[1]) + args, kw_out

    return self_method(*args, **kwargs)


def _threeargs_replacer(args, kwargs, dispatchables):
    def self_method(a, b, c, out=None, *args, **kwargs):
        kw_out = kwargs.copy()
        if out is not None:
            kw_out["out"] = dispatchables[3:]
        return (
                dispatchables[0],
                dispatchables[1],
                dispatchables[2]
            ) + args, kw_out

    return self_method(*args, **kwargs)


def _fourargs_replacer(args, kwargs, dispatchables):
    def self_method(a, b, c, d, out=None, *args, **kwargs):
        kw_out = kwargs.copy()
        if out is not None:
            kw_out["out"] = dispatchables[4:]
        return (
                dispatchables[0],
                dispatchables[1],
                dispatchables[2],
                dispatchables[3]
            ) + args, kw_out

    return self_method(*args, **kwargs)


def _fiveargs_replacer(args, kwargs, dispatchables):
    def self_method(a, b, c, d, e, out=None, *args, **kwargs):
        kw_out = kwargs.copy()
        if out is not None:
            kw_out["out"] = dispatchables[5:]
        return (
                dispatchables[0],
                dispatchables[1],
                dispatchables[2],
                dispatchables[3],
                dispatchables[4]
            ) + args, kw_out

    return self_method(*args, **kwargs)


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def ai_zeros(nt):
    return ()


@_create_special(_x_n_kw_k_replacer)
@all_of_type(ndarray)
@_get_docs
def assoc_laguerre(x, n, k=0.0):
    return (x, n, k)


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def bei_zeros(nt):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def beip_zeros(nt):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def ber_zeros(nt):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def bernoulli(n):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def berp_zeros(nt):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def bi_zeros(nt):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def roots_chebyc(n, mu=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def roots_gegenbauer(n, alpha, mu=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def chebyc(n, monic=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def chebys(n, monic=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def chebyt(n, monic=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def chebyu(n, monic=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def clpmn(m, n, z, type=3):
    return ()


@_create_special(_N_k_replacer)
@all_of_type(ndarray)
@_get_docs
def comb(N, k, exact=False, repetition=False, legacy=True):
    return (N, k)


@_create_special(_x_n_replacer)
@all_of_type(ndarray)
@_get_docs
def diric(x, n):
    return (x,)


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def ellip_harm(h2, k2, n, p, s, signm=1, signn=1):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def ellip_harm_2(h2, k2, n, p, s):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def ellip_normal(h2, k2, n, p):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def erf_zeros(nt):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def euler(n):
    return ()


@_create_special(_n_replacer)
@all_of_type(ndarray)
@_get_docs
def factorial(n, exact=False):
    return (n,)


@_create_special(_n_replacer)
@all_of_type(ndarray)
@_get_docs
def factorial2(n, exact=False):
    return (n,)


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def factorialk(n, k, exact=True):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def fresnel_zeros(nt):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def fresnelc_zeros(nt):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def fresnels_zeros(nt):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def gegenbauer(n, alpha, monic=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def genlaguerre(n, alpha, monic=False):
    return ()


@_create_special(_v_z_kw_n_replacer)
@all_of_type(ndarray)
@_get_docs
def h1vp(v, z, n=1):
    return (v, z)


@_create_special(_v_z_kw_n_replacer)
@all_of_type(ndarray)
@_get_docs
def h2vp(v, z, n=1):
    return (v, z, n)


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def roots_hermite(n, mu=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def roots_hermitenorm(n, mu=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def hermite(n, monic=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def hermitenorm(n, monic=False):
    return ()


@_create_special(_v_z_kw_n_replacer)
@all_of_type(ndarray)
@_get_docs
def ivp(v, z, n=1):
    return (v, z)


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def roots_jacobi(n, alpha, beta, mu=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def jacobi(n, alpha, beta, monic=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def jn_zeros(n, nt):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def jnjnp_zeros(nt):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def jnp_zeros(n, nt):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def jnyn_zeros(n, nt):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def roots_sh_jacobi(n, p1, q1, mu=False):
    return ()


@_create_special(_v_z_kw_n_replacer)
@all_of_type(ndarray)
@_get_docs
def jvp(v, z, n=1):
    return (v, z)


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def kei_zeros(nt):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def keip_zeros(nt):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def kelvin_zeros(nt):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def ker_zeros(nt):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def kerp_zeros(nt):
    return ()


@_create_special(_v_z_kw_n_replacer)
@all_of_type(ndarray)
@_get_docs
def kvp(v, z, n=1):
    return (v, z)


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def roots_laguerre(n, mu=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def roots_genlaguerre(n, alpha, mu=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def laguerre(n, monic=False):
    return ()


@_create_special(_z_kw_k_replacer)
@all_of_type(ndarray)
@_get_docs
def lambertw(z, k=0, tol=1e-08):
    return (z, )


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def legendre(n, monic=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def lmbda(v, x):
    return ()


@_create_special(_x_replacer)
@all_of_type(ndarray)
@_get_docs
def log_softmax(x, axis=None):
    return (x, )


@_create_special(_a_kw_axis_b_replacer)
@all_of_type(ndarray)
@_get_docs
def logsumexp(a, axis=None, b=None, keepdims=False, return_sign=False):
    return (a, b)


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def lpmn(m, n, z):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def lpn(n, z):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def lqmn(m, n, z):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def lqn(n, z):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def mathieu_even_coef(m, q):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def mathieu_odd_coef(m, q):
    return ()


@_create_special(_a_replacer)
@all_of_type(ndarray)
@_get_docs
def multigammaln(a, d):
    return (a, )


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def obl_cv_seq(m, n, c):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def roots_legendre(n, mu=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def p_roots(n, mu=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def pbdn_seq(n, z):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def pbdv_seq(v, x):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def pbvv_seq(v, x):
    return ()


@_create_special(_N_k_replacer)
@all_of_type(ndarray)
@_get_docs
def perm(N, k, exact=False):
    return (N, k)


@_create_special(_n_x_replacer)
@all_of_type(ndarray)
@_get_docs
def polygamma(n, x):
    return (n, x)


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def pro_cv_seq(m, n, c):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def roots_sh_legendre(n, mu=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def riccati_jn(n, x):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def riccati_yn(n, x):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def roots_chebyt(n, mu=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def roots_chebyu(n, mu=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def roots_sh_chebyt(n, mu=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def roots_sh_chebyu(n, mu=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def roots_chebys(n, mu=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def sh_chebyt(n, monic=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def sh_chebyu(n, monic=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def sh_jacobi(n, p, q, monic=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def sh_legendre(n, monic=False):
    return ()


@_create_special(_x_replacer)
@all_of_type(ndarray)
@_get_docs
def sinc(x):
    return (x,)


@_create_special(_x_replacer)
@all_of_type(ndarray)
@_get_docs
def softmax(x, axis=None):
    return (x,)


@_create_special(_n_z_replacer)
@all_of_type(ndarray)
@_get_docs
def spherical_in(n, z, derivative=False):
    return (n, z)


@_create_special(_n_z_replacer)
@all_of_type(ndarray)
@_get_docs
def spherical_jn(n, z, derivative=False):
    return (n, z)


@_create_special(_n_z_replacer)
@all_of_type(ndarray)
@_get_docs
def spherical_kn(n, z, derivative=False):
    return (n, z)


@_create_special(_n_z_replacer)
@all_of_type(ndarray)
@_get_docs
def spherical_yn(n, z, derivative=False):
    return (n, z)


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def y0_zeros(nt, complex=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def y1_zeros(nt, complex=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def y1p_zeros(nt, complex=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def yn_zeros(n, nt):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def ynp_zeros(n, nt):
    return ()


@_create_special(_v_z_kw_n_replacer)
@all_of_type(ndarray)
@_get_docs
def yvp(v, z, n=1):
    return (v, z)


@_create_special(_x_kw_q_out_replacer)
@all_of_type(ndarray)
@_get_docs
def zeta(x, q=None, out=None):
    return (x, q, out)
