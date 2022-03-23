import functools
import numpy as np
from numpy import dtype, ndarray
from scipy._lib.uarray import all_of_type, Dispatchable
from scipy.special import _api
from ._multimethods import _get_docs, _create_special


__all__ = [
    'agm',
    'airy',
    'airye',
    'bdtr',
    'bdtrc',
    'bdtri',
    'bdtrik',
    'bdtrin',
    'bei',
    'beip',
    'ber',
    'berp',
    'besselpoly',
    'beta',
    'betainc',
    'betaincinv',
    'betaln',
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
    'chndtr',
    'chndtridf',
    'chndtrinc',
    'chndtrix',
    'cosdg',
    'cosm1',
    'cotdg',
    'dawsn',
    'digamma',
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
    'erfc',
    'erfcinv',
    'erfcx',
    'erfi',
    'erfinv',
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
    'fdtr',
    'fdtrc',
    'fdtri',
    'fdtridfd',
    'fresnel',
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
    'hankel1',
    'hankel1e',
    'hankel2',
    'hankel2e',
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
    'j0',
    'j1',
    'jn',
    'jv',
    'jve',
    'k0',
    'k0e',
    'k1',
    'k1e',
    'kei',
    'keip',
    'kelvin',
    'ker',
    'kerp',
    'kl_div',
    'kn',
    'kolmogi',
    'kolmogorov',
    'kv',
    'kve',
    'log1p',
    'log_expit',
    'log_ndtr',
    'loggamma',
    'logit',
    'lpmv',
    'mathieu_a',
    'mathieu_b',
    'mathieu_cem',
    'mathieu_modcem1',
    'mathieu_modcem2',
    'mathieu_modsem1',
    'mathieu_modsem2',
    'mathieu_sem',
    'modfresnelm',
    'modfresnelp',
    'modstruve',
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
    'obl_rad1',
    'obl_rad1_cv',
    'obl_rad2',
    'obl_rad2_cv',
    'owens_t',
    'pbdv',
    'pbvv',
    'pbwa',
    'pdtr',
    'pdtrc',
    'pdtri',
    'pdtrik',
    'poch',
    'pro_ang1',
    'pro_ang1_cv',
    'pro_cv',
    'pro_rad1',
    'pro_rad1_cv',
    'pro_rad2',
    'pro_rad2_cv',
    'pseudo_huber',
    'psi',
    'radian',
    'rel_entr',
    'rgamma',
    'round',
    'shichi',
    'sici',
    'sindg',
    'smirnov',
    'smirnovi',
    'spence',
    'sph_harm',
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
    'y1',
    'yn',
    'yv',
    'yve',
    'zetac',

]


_mark_output = functools.partial(Dispatchable,
                                 dispatch_type=ndarray,
                                 coercible=False)


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


########################### auto generated ufuncs #############################


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


########################## manually added ufuncs ##############################


@_create_special(_twoargs_replacer)
@all_of_type(ndarray)
@_get_docs
def jn(a, b, out=None, **kwargs):
    return a, b, _mark_output(out)


@_create_special(_onearg_replacer)
@all_of_type(ndarray)
@_get_docs
def digamma(a, out=None, **kwargs):
    return a, _mark_output(out)
