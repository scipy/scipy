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
    'cbrt',
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
    'h1vp',
    'h2vp',
    'hankel1',
    'hankel1e',
    'hankel2',
    'hankel2e',
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
    'jacobi',
    'jn',
    'jn_zeros',
    'jnjnp_zeros',
    'jnp_zeros',
    'jnyn_zeros',
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
    'tandg',
    'tklmbda',
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

_mark_output = functools.partial(Dispatchable,
                                 dispatch_type=ndarray,
                                 coercible=False)

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
        if out:
            kw_out["out"] = dispatchables[1:]
        return (dispatchables[0],) + args, kw_out

    return self_method(*args, **kwargs)

def _twoargs_replacer(args, kwargs, dispatchables):
    def self_method(a, b, out=None, *args, **kwargs):
        kw_out = kwargs.copy()
        if out:
            kw_out["out"] = dispatchables[2:]
        return (dispatchables[0], dispatchables[1]) + args, kw_out

    return self_method(*args, **kwargs)


def _threeargs_replacer(args, kwargs, dispatchables):
    def self_method(a, b, c, out=None, *args, **kwargs):
        kw_out = kwargs.copy()
        if out:
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
        if out:
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
        if out:
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


##################################### ufuncs ###################################


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def agm(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def airy(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def airye(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def bdtr(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def bdtrc(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def bdtri(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def bdtrik(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def bdtrin(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def bei(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def beip(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def ber(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def berp(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def besselpoly(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def beta(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def betainc(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def betaincinv(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def betaln(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def binom(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def boxcox(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def boxcox1p(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def btdtr(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def btdtri(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def btdtria(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def btdtrib(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def cbrt(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def chdtr(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def chdtrc(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def chdtri(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def chdtriv(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def chndtr(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def chndtridf(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def chndtrinc(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def chndtrix(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def cosdg(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def cosm1(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def cotdg(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def dawsn(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def ellipe(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def ellipeinc(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def ellipj(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def ellipk(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def ellipkinc(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def ellipkm1(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def elliprc(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def elliprd(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def elliprf(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def elliprg(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_fourargs_replacer)
@all_of_type(ndarray)
@_get_docs
def elliprj(a, b, c, d, out=None, **kwargs):
    return a, b, c, d, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def entr(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def erf(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def erfc(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def erfcinv(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def erfcx(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def erfi(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def erfinv(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def eval_chebyc(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def eval_chebys(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def eval_chebyt(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def eval_chebyu(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def eval_gegenbauer(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def eval_genlaguerre(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def eval_hermite(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def eval_hermitenorm(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_fourargs_replacer)
@all_of_type(ndarray)
@_get_docs
def eval_jacobi(a, b, c, d, out=None, **kwargs):
    return a, b, c, d, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def eval_laguerre(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def eval_legendre(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def eval_sh_chebyt(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def eval_sh_chebyu(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_fourargs_replacer)
@all_of_type(ndarray)
@_get_docs
def eval_sh_jacobi(a, b, c, d, out=None, **kwargs):
    return a, b, c, d, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def eval_sh_legendre(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def exp1(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def exp10(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def exp2(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def expi(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def expit(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def expm1(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def expn(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def exprel(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def fdtr(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def fdtrc(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def fdtri(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def fdtridfd(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def fresnel(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def gamma(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def gammainc(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def gammaincc(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def gammainccinv(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def gammaincinv(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def gammaln(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def gammasgn(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def gdtr(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def gdtrc(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def gdtria(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def gdtrib(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def gdtrix(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def hankel1(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def hankel1e(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def hankel2(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def hankel2e(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def huber(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def hyp0f1(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def hyp1f1(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_fourargs_replacer)
@all_of_type(ndarray)
@_get_docs
def hyp2f1(a, b, c, d, out=None, **kwargs):
    return a, b, c, d, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def hyperu(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def i0(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def i0e(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def i1(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def i1e(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def inv_boxcox(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def inv_boxcox1p(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def it2i0k0(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def it2j0y0(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def it2struve0(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def itairy(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def iti0k0(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def itj0y0(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def itmodstruve0(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def itstruve0(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def iv(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def ive(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def j0(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def j1(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def jv(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def jn(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def jve(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def k0(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def k0e(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def k1(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def k1e(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def kei(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def keip(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def kelvin(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def ker(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def kerp(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def kl_div(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def kn(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def kolmogi(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def kolmogorov(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def kv(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def kve(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def log1p(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def log_expit(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def log_ndtr(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def loggamma(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def logit(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def lpmv(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def mathieu_a(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def mathieu_b(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def mathieu_cem(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def mathieu_modcem1(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def mathieu_modcem2(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def mathieu_modsem1(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def mathieu_modsem2(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def mathieu_sem(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def modfresnelm(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def modfresnelp(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def modstruve(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def nbdtr(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def nbdtrc(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def nbdtri(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def nbdtrik(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def nbdtrin(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_fourargs_replacer)
@all_of_type(ndarray)
@_get_docs
def ncfdtr(a, b, c, d, out=None, **kwargs):
    return a, b, c, d, _mark_output(out)


@_create_special(_fourargs_replacer)
@all_of_type(ndarray)
@_get_docs
def ncfdtri(a, b, c, d, out=None, **kwargs):
    return a, b, c, d, _mark_output(out)


@_create_special(_fourargs_replacer)
@all_of_type(ndarray)
@_get_docs
def ncfdtridfd(a, b, c, d, out=None, **kwargs):
    return a, b, c, d, _mark_output(out)


@_create_special(_fourargs_replacer)
@all_of_type(ndarray)
@_get_docs
def ncfdtridfn(a, b, c, d, out=None, **kwargs):
    return a, b, c, d, _mark_output(out)


@_create_special(_fourargs_replacer)
@all_of_type(ndarray)
@_get_docs
def ncfdtrinc(a, b, c, d, out=None, **kwargs):
    return a, b, c, d, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def nctdtr(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def nctdtridf(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def nctdtrinc(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def nctdtrit(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def ndtr(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def ndtri(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def ndtri_exp(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def nrdtrimn(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def nrdtrisd(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_fourargs_replacer)
@all_of_type(ndarray)
@_get_docs
def obl_ang1(a, b, c, d, out=None, **kwargs):
    return a, b, c, d, _mark_output(out)


@_create_special(_fiveargs_replacer)
@all_of_type(ndarray)
@_get_docs
def obl_ang1_cv(a, b, c, d, e, out=None, **kwargs):
    return a, b, c, d, e, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def obl_cv(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_fourargs_replacer)
@all_of_type(ndarray)
@_get_docs
def obl_rad1(a, b, c, d, out=None, **kwargs):
    return a, b, c, d, _mark_output(out)


@_create_special(_fiveargs_replacer)
@all_of_type(ndarray)
@_get_docs
def obl_rad1_cv(a, b, c, d, e, out=None, **kwargs):
    return a, b, c, d, e, _mark_output(out)


@_create_special(_fourargs_replacer)
@all_of_type(ndarray)
@_get_docs
def obl_rad2(a, b, c, d, out=None, **kwargs):
    return a, b, c, d, _mark_output(out)


@_create_special(_fiveargs_replacer)
@all_of_type(ndarray)
@_get_docs
def obl_rad2_cv(a, b, c, d, e, out=None, **kwargs):
    return a, b, c, d, e, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def owens_t(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def pbdv(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def pbvv(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def pbwa(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def pdtr(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def pdtrc(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def pdtri(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def pdtrik(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def poch(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_fourargs_replacer)
@all_of_type(ndarray)
@_get_docs
def pro_ang1(a, b, c, d, out=None, **kwargs):
    return a, b, c, d, _mark_output(out)


@_create_special(_fiveargs_replacer)
@all_of_type(ndarray)
@_get_docs
def pro_ang1_cv(a, b, c, d, e, out=None, **kwargs):
    return a, b, c, d, e, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def pro_cv(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_fourargs_replacer)
@all_of_type(ndarray)
@_get_docs
def pro_rad1(a, b, c, d, out=None, **kwargs):
    return a, b, c, d, _mark_output(out)


@_create_special(_fiveargs_replacer)
@all_of_type(ndarray)
@_get_docs
def pro_rad1_cv(a, b, c, d, e, out=None, **kwargs):
    return a, b, c, d, e, _mark_output(out)


@_create_special(_fourargs_replacer)
@all_of_type(ndarray)
@_get_docs
def pro_rad2(a, b, c, d, out=None, **kwargs):
    return a, b, c, d, _mark_output(out)


@_create_special(_fiveargs_replacer)
@all_of_type(ndarray)
@_get_docs
def pro_rad2_cv(a, b, c, d, e, out=None, **kwargs):
    return a, b, c, d, e, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def pseudo_huber(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def psi(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def radian(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def rel_entr(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def rgamma(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def round(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def shichi(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def sici(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def sindg(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def smirnov(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def smirnovi(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def spence(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_fourargs_replacer)
@all_of_type(ndarray)
@_get_docs
def sph_harm(a, b, c, d, out=None, **kwargs):
    return a, b, c, d, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def stdtr(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def stdtridf(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def stdtrit(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def struve(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def tandg(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def tklmbda(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def voigt_profile(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def wofz(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_threeargs_replacer)
@all_of_type(ndarray)
@_get_docs
def wright_bessel(a, b, c, out=None, **kwargs):
    return a, b, c, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def wrightomega(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def xlog1py(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def xlogy(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def y0(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def y1(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def yn(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def yv(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def yve(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def zetac(a, out=None, **kwargs):
    return a, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def digamma(a, out=None, **kwargs):
    return a, _mark_output(out)
