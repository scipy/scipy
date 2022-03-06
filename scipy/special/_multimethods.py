import functools
import numpy as np
from numpy import dtype, ndarray
from scipy._lib.uarray import Dispatchable, all_of_type, create_multimethod
from scipy.special import _api
from scipy.special._backend import scalar_tuple_callable_array


__all__ = [
    'agm',
    'ai_zeros',
    'airy',
    'airye',
    'assoc_laguerre',
    'bdtr',
    'bdtrc',
    'bdtri',
    'bdtrik',
    'bdtrin',
    'bei',
    'bei_zeros',
    'beip',
    'beip_zeros',
    'ber',
    'ber_zeros',
    'bernoulli',
    'berp',
    'berp_zeros',
    'besselpoly',
    'beta',
    'betainc',
    'betaincinv',
    'betaln',
    'bi_zeros',
    'binom',
    'boxcox',
    'boxcox1p',
    'btdtr',
    'btdtri',
    'btdtria',
    'btdtrib',
    'c_roots',
    'cbrt',
    'cg_roots',
    'chdtr',
    'chdtrc',
    'chdtri',
    'chdtriv',
    'chebyc',
    'chebys',
    'chebyt',
    'chebyu',
    'chndtr',
    'chndtridf',
    'chndtrinc',
    'chndtrix',
    'clpmn',
    'comb',
    'cosdg',
    'cosm1',
    'cotdg',
    'dawsn',
    'digamma',
    'diric',
    'ellip_harm',
    'ellip_harm_2',
    'ellip_normal',
    'ellipe',
    'ellipeinc',
    'ellipj',
    'ellipk',
    'ellipkinc',
    'ellipkm1',
    'elliprc',
    'elliprd',
    'elliprf',
    'elliprg',
    'elliprj',
    'entr',
    'erf',
    'erf_zeros',
    'erfc',
    'erfcinv',
    'erfcx',
    'erfi',
    'erfinv',
    'errstate',
    'euler',
    'eval_chebyc',
    'eval_chebys',
    'eval_chebyt',
    'eval_chebyu',
    'eval_gegenbauer',
    'eval_genlaguerre',
    'eval_hermite',
    'eval_hermitenorm',
    'eval_jacobi',
    'eval_laguerre',
    'eval_legendre',
    'eval_sh_chebyt',
    'eval_sh_chebyu',
    'eval_sh_jacobi',
    'eval_sh_legendre',
    'exp1',
    'exp10',
    'exp2',
    'expi',
    'expit',
    'expm1',
    'expn',
    'exprel',
    'factorial',
    'factorial2',
    'factorialk',
    'fdtr',
    'fdtrc',
    'fdtri',
    'fdtridfd',
    'fresnel',
    'fresnel_zeros',
    'fresnelc_zeros',
    'fresnels_zeros',
    'gamma',
    'gammainc',
    'gammaincc',
    'gammainccinv',
    'gammaincinv',
    'gammaln',
    'gammasgn',
    'gdtr',
    'gdtrc',
    'gdtria',
    'gdtrib',
    'gdtrix',
    'gegenbauer',
    'genlaguerre',
    'geterr',
    'h1vp',
    'h2vp',
    'h_roots',
    'hankel1',
    'hankel1e',
    'hankel2',
    'hankel2e',
    'he_roots',
    'hermite',
    'hermitenorm',
    'huber',
    'hyp0f1',
    'hyp1f1',
    'hyp2f1',
    'hyperu',
    'i0',
    'i0e',
    'i1',
    'i1e',
    'inv_boxcox',
    'inv_boxcox1p',
    'it2i0k0',
    'it2j0y0',
    'it2struve0',
    'itairy',
    'iti0k0',
    'itj0y0',
    'itmodstruve0',
    'itstruve0',
    'iv',
    'ive',
    'ivp',
    'j0',
    'j1',
    'j_roots',
    'jacobi',
    'jn',
    'jn_zeros',
    'jnjnp_zeros',
    'jnp_zeros',
    'jnyn_zeros',
    'js_roots',
    'jv',
    'jve',
    'jvp',
    'k0',
    'k0e',
    'k1',
    'k1e',
    'kei',
    'kei_zeros',
    'keip',
    'keip_zeros',
    'kelvin',
    'kelvin_zeros',
    'ker',
    'ker_zeros',
    'kerp',
    'kerp_zeros',
    'kl_div',
    'kn',
    'kolmogi',
    'kolmogorov',
    'kv',
    'kve',
    'kvp',
    'l_roots',
    'la_roots',
    'laguerre',
    'lambertw',
    'legendre',
    'lmbda',
    'log1p',
    'log_expit',
    'log_ndtr',
    'log_softmax',
    'loggamma',
    'logit',
    'logsumexp',
    'lpmn',
    'lpmv',
    'lpn',
    'lqmn',
    'lqn',
    'mathieu_a',
    'mathieu_b',
    'mathieu_cem',
    'mathieu_even_coef',
    'mathieu_modcem1',
    'mathieu_modcem2',
    'mathieu_modsem1',
    'mathieu_modsem2',
    'mathieu_odd_coef',
    'mathieu_sem',
    'modfresnelm',
    'modfresnelp',
    'modstruve',
    'multigammaln',
    'nbdtr',
    'nbdtrc',
    'nbdtri',
    'nbdtrik',
    'nbdtrin',
    'ncfdtr',
    'ncfdtri',
    'ncfdtridfd',
    'ncfdtridfn',
    'ncfdtrinc',
    'nctdtr',
    'nctdtridf',
    'nctdtrinc',
    'nctdtrit',
    'ndtr',
    'ndtri',
    'ndtri_exp',
    'nrdtrimn',
    'nrdtrisd',
    'obl_ang1',
    'obl_ang1_cv',
    'obl_cv',
    'obl_cv_seq',
    'obl_rad1',
    'obl_rad1_cv',
    'obl_rad2',
    'obl_rad2_cv',
    'owens_t',
    'p_roots',
    'pbdn_seq',
    'pbdv',
    'pbdv_seq',
    'pbvv',
    'pbvv_seq',
    'pbwa',
    'pdtr',
    'pdtrc',
    'pdtri',
    'pdtrik',
    'perm',
    'poch',
    'polygamma',
    'pro_ang1',
    'pro_ang1_cv',
    'pro_cv',
    'pro_cv_seq',
    'pro_rad1',
    'pro_rad1_cv',
    'pro_rad2',
    'pro_rad2_cv',
    'ps_roots',
    'pseudo_huber',
    'psi',
    'radian',
    'rel_entr',
    'rgamma',
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
    'round',
    's_roots',
    'seterr',
    'sh_chebyt',
    'sh_chebyu',
    'sh_jacobi',
    'sh_legendre',
    'shichi',
    'sici',
    'sinc',
    'sindg',
    'smirnov',
    'smirnovi',
    'softmax',
    'spence',
    'sph_harm',
    'spherical_in',
    'spherical_jn',
    'spherical_kn',
    'spherical_yn',
    'stdtr',
    'stdtridf',
    'stdtrit',
    'struve',
    't_roots',
    'tandg',
    'tklmbda',
    'ts_roots',
    'u_roots',
    'us_roots',
    'voigt_profile',
    'wofz',
    'wright_bessel',
    'wrightomega',
    'xlog1py',
    'xlogy',
    'y0',
    'y0_zeros',
    'y1',
    'y1_zeros',
    'y1p_zeros',
    'yn',
    'yn_zeros',
    'ynp_zeros',
    'yv',
    'yve',
    'yvp',
    'zeta',
    'zetac',
]

_create_special = functools.partial(
    create_multimethod, domain="numpy.scipy.special"
)


_mark_scalar_tuple_callable_array = functools.partial(
                            Dispatchable,
                            dispatch_type=scalar_tuple_callable_array,
                            coercible=True
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
def comb(N, k, exact=False, repetition=False):
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


@_create_special(_m_n_replacer)
@all_of_type(ndarray)
@_get_docs
def obl_cv_seq(m, n, c):
    return (n,)


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def roots_legendre(n, mu=False):
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


@_create_special(_m_n_replacer)
@all_of_type(ndarray)
@_get_docs
def pro_cv_seq(m, n, c):
    return (n,)


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
def riccati_yn():
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def roots_chebyc(n, mu=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def roots_chebys(n, mu=False):
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
def roots_gegenbauer(n, alpha, mu=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def roots_genlaguerre(n, alpha, mu=False):
    return ()


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
def roots_jacobi(n, alpha, beta, mu=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def roots_laguerre(n, mu=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def roots_legendre(n, mu=False):
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
def roots_sh_jacobi(n, p1, q1, mu=False):
    return ()


@_create_special(_identity_replacer)
@all_of_type(ndarray)
@_get_docs
def roots_sh_legendre(n, mu=False):
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
