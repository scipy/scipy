# Type declarations for functions generated dynamically in
# scipy.special._support_alternative_backends.


from collections.abc import Callable, Mapping
from dataclasses import dataclass
from types import ModuleType
from typing import Any

from ._ufunc_tools import _UFuncWrapper


class _CallableUFuncWrapper(_UFuncWrapper):
    """The callable wrapper dynamically constructed from _UFuncWrapper."""

    def __call__(self, *args: Any, **kwargs: Any) -> Any: ...


@dataclass
class _FuncInfo:
    func: Callable[..., Any]
    arg_names: list[str]
    xp_capabilities: Callable[..., Any] | None = ...
    generic_impl: Callable[..., Any] | None = ...
    alt_names_map: dict[str, str] | None = ...
    alt_arg_orders_map: dict[str, tuple[int, ...]] | None = ...
    int_only: tuple[bool, ...] | dict[str, tuple[bool, ...]] | None = ...
    positive_only: (
        bool | tuple[bool, ...] | Mapping[str, bool | tuple[bool, ...]]
    ) = ...
    is_ufunc: bool = ...
    is_elementwise: bool = ...
    shape_mapper: Callable[..., Any] | None = ...
    python_int_only: (
        dict[str, tuple[bool, ...]] | tuple[bool, ...] | None
    ) = ...
    scalar_or_0d_only: (
        dict[str, tuple[bool, ...]] | tuple[bool, ...] | None
    ) = ...
    test_large_ints: bool = ...
    produces_0d: bool = ...
    torch_native: bool = ...
    backends_with_func_in_xp: tuple[str, ...] = ...
    public: bool = ...

    @property
    def name(self) -> str: ...
    @property
    def n_args(self) -> int: ...
    @property
    def wrapper(self) -> _CallableUFuncWrapper | Callable[..., Any]: ...
    def _wrapper_for(self, xp: ModuleType) -> Callable[..., Any]: ...


_special_funcs: tuple[_FuncInfo, ...]
_FUNC_INFO: dict[str, _FuncInfo]
_aliases: dict[str, str]
registered: set[str]


def get_function(name: str) -> _CallableUFuncWrapper | Callable[..., Any]: ...
def __dir__() -> list[str]: ...


# Non-exported helpers.
def _special_namespace_for(*args: Any, **kwargs: Any) -> Any: ...
def _ufunc_kwargs_extra_note(*args: Any, **kwargs: Any) -> Any: ...
def _get_native_func(*args: Any, **kwargs: Any) -> Any: ...
def _rel_entr(*args: Any, **kwargs: Any) -> Any: ...
def _xlogy(*args: Any, **kwargs: Any) -> Any: ...
def _chdtr(*args: Any, **kwargs: Any) -> Any: ...
def _chdtrc(*args: Any, **kwargs: Any) -> Any: ...
def _betaincc(*args: Any, **kwargs: Any) -> Any: ...
def _stdtr(*args: Any, **kwargs: Any) -> Any: ...
def _stdtrit(*args: Any, **kwargs: Any) -> Any: ...
def _ndtri_exp(*args: Any, **kwargs: Any) -> Any: ...
def __riemann_zeta(*args: Any, **kwargs: Any) -> Any: ...
def _poisson_binom_cdf_shape_mapper(*args: Any, **kwargs: Any) -> Any: ...
def _needs_betainc(*args: Any, **kwargs: Any) -> Any: ...


# Public functions re-exported by scipy.special.__init__.py.
agm: _CallableUFuncWrapper
airy: _CallableUFuncWrapper
airye: _CallableUFuncWrapper
bdtr: _CallableUFuncWrapper
bdtrc: _CallableUFuncWrapper
bdtri: _CallableUFuncWrapper
bdtrik: _CallableUFuncWrapper
bdtrin: _CallableUFuncWrapper
bei: _CallableUFuncWrapper
beip: _CallableUFuncWrapper
ber: _CallableUFuncWrapper
berp: _CallableUFuncWrapper
besselpoly: _CallableUFuncWrapper
beta: _CallableUFuncWrapper
betainc: _CallableUFuncWrapper
betaincc: _CallableUFuncWrapper
betainccinv: _CallableUFuncWrapper
betaincinv: _CallableUFuncWrapper
betaln: _CallableUFuncWrapper
binom: _CallableUFuncWrapper
boxcox: _CallableUFuncWrapper
boxcox1p: _CallableUFuncWrapper
btdtria: _CallableUFuncWrapper
btdtrib: _CallableUFuncWrapper
cbrt: _CallableUFuncWrapper
chdtr: _CallableUFuncWrapper
chdtrc: _CallableUFuncWrapper
chdtri: _CallableUFuncWrapper
chdtriv: _CallableUFuncWrapper
chndtr: _CallableUFuncWrapper
chndtridf: _CallableUFuncWrapper
chndtrinc: _CallableUFuncWrapper
chndtrix: _CallableUFuncWrapper
cosdg: _CallableUFuncWrapper
cosm1: _CallableUFuncWrapper
cotdg: _CallableUFuncWrapper
dawsn: _CallableUFuncWrapper
digamma: _CallableUFuncWrapper
digammainv: _CallableUFuncWrapper
ellip_harm: Callable[..., Any]
ellip_harm_2: Callable[..., Any]
ellip_normal: Callable[..., Any]
ellipe: _CallableUFuncWrapper
ellipeinc: _CallableUFuncWrapper
ellipj: _CallableUFuncWrapper
ellipk: _CallableUFuncWrapper
ellipkinc: _CallableUFuncWrapper
ellipkm1: _CallableUFuncWrapper
elliprc: _CallableUFuncWrapper
elliprd: _CallableUFuncWrapper
elliprf: _CallableUFuncWrapper
elliprg: _CallableUFuncWrapper
elliprj: _CallableUFuncWrapper
entr: _CallableUFuncWrapper
erf: _CallableUFuncWrapper
erfc: _CallableUFuncWrapper
erfcinv: _CallableUFuncWrapper
erfcx: _CallableUFuncWrapper
erfi: _CallableUFuncWrapper
erfinv: _CallableUFuncWrapper
eval_chebyc: _CallableUFuncWrapper
eval_chebys: _CallableUFuncWrapper
eval_chebyt: _CallableUFuncWrapper
eval_chebyu: _CallableUFuncWrapper
eval_gegenbauer: _CallableUFuncWrapper
eval_genlaguerre: _CallableUFuncWrapper
eval_hermite: _CallableUFuncWrapper
eval_hermitenorm: _CallableUFuncWrapper
eval_jacobi: _CallableUFuncWrapper
eval_laguerre: _CallableUFuncWrapper
eval_legendre: _CallableUFuncWrapper
eval_sh_chebyt: _CallableUFuncWrapper
eval_sh_chebyu: _CallableUFuncWrapper
eval_sh_jacobi: _CallableUFuncWrapper
eval_sh_legendre: _CallableUFuncWrapper
exp1: _CallableUFuncWrapper
exp10: _CallableUFuncWrapper
exp2: _CallableUFuncWrapper
expi: _CallableUFuncWrapper
expit: _CallableUFuncWrapper
expm1: _CallableUFuncWrapper
expn: _CallableUFuncWrapper
exprel: _CallableUFuncWrapper
fdtr: _CallableUFuncWrapper
fdtrc: _CallableUFuncWrapper
fdtri: _CallableUFuncWrapper
fdtridfd: _CallableUFuncWrapper
fresnel: _CallableUFuncWrapper
gamma: _CallableUFuncWrapper
gammainc: _CallableUFuncWrapper
gammaincc: _CallableUFuncWrapper
gammainccinv: _CallableUFuncWrapper
gammaincinv: _CallableUFuncWrapper
gammaln: _CallableUFuncWrapper
gammasgn: _CallableUFuncWrapper
gdtr: _CallableUFuncWrapper
gdtrc: _CallableUFuncWrapper
gdtria: _CallableUFuncWrapper
gdtrib: _CallableUFuncWrapper
gdtrix: _CallableUFuncWrapper
hankel1: _CallableUFuncWrapper
hankel1e: _CallableUFuncWrapper
hankel2: _CallableUFuncWrapper
hankel2e: _CallableUFuncWrapper
huber: _CallableUFuncWrapper
hyp0f1: _CallableUFuncWrapper
hyp1f1: _CallableUFuncWrapper
hyp2f1: _CallableUFuncWrapper
hyperu: _CallableUFuncWrapper
i0: _CallableUFuncWrapper
i0e: _CallableUFuncWrapper
i1: _CallableUFuncWrapper
i1e: _CallableUFuncWrapper
inv_boxcox: _CallableUFuncWrapper
inv_boxcox1p: _CallableUFuncWrapper
it2i0k0: _CallableUFuncWrapper
it2j0y0: _CallableUFuncWrapper
it2struve0: _CallableUFuncWrapper
itairy: _CallableUFuncWrapper
iti0k0: _CallableUFuncWrapper
itj0y0: _CallableUFuncWrapper
itmodstruve0: _CallableUFuncWrapper
itstruve0: _CallableUFuncWrapper
iv: _CallableUFuncWrapper
ive: _CallableUFuncWrapper
j0: _CallableUFuncWrapper
j1: _CallableUFuncWrapper
jn: _CallableUFuncWrapper
jv: _CallableUFuncWrapper
jve: _CallableUFuncWrapper
k0: _CallableUFuncWrapper
k0e: _CallableUFuncWrapper
k1: _CallableUFuncWrapper
k1e: _CallableUFuncWrapper
kei: _CallableUFuncWrapper
keip: _CallableUFuncWrapper
kelvin: _CallableUFuncWrapper
ker: _CallableUFuncWrapper
kerp: _CallableUFuncWrapper
kl_div: _CallableUFuncWrapper
kn: _CallableUFuncWrapper
kolmogi: _CallableUFuncWrapper
kolmogorov: _CallableUFuncWrapper
kv: _CallableUFuncWrapper
kve: _CallableUFuncWrapper
lambertw: Callable[..., Any]
log1p: _CallableUFuncWrapper
log_expit: _CallableUFuncWrapper
log_gammainc: _CallableUFuncWrapper
log_gammaincc: _CallableUFuncWrapper
log_ndtr: _CallableUFuncWrapper
log_wright_bessel: _CallableUFuncWrapper
loggamma: _CallableUFuncWrapper
logit: _CallableUFuncWrapper
lpmv: _CallableUFuncWrapper
mathieu_a: _CallableUFuncWrapper
mathieu_b: _CallableUFuncWrapper
mathieu_ce: _CallableUFuncWrapper
mathieu_cem: _CallableUFuncWrapper
mathieu_modcem1: _CallableUFuncWrapper
mathieu_modcem2: _CallableUFuncWrapper
mathieu_modsem1: _CallableUFuncWrapper
mathieu_modsem2: _CallableUFuncWrapper
mathieu_se: _CallableUFuncWrapper
mathieu_sem: _CallableUFuncWrapper
modfresnelm: _CallableUFuncWrapper
modfresnelp: _CallableUFuncWrapper
modstruve: _CallableUFuncWrapper
multigammaln: Callable[..., Any]
nbdtr: _CallableUFuncWrapper
nbdtrc: _CallableUFuncWrapper
nbdtri: _CallableUFuncWrapper
nbdtrik: _CallableUFuncWrapper
nbdtrin: _CallableUFuncWrapper
ncfdtr: _CallableUFuncWrapper
ncfdtri: _CallableUFuncWrapper
ncfdtridfd: _CallableUFuncWrapper
ncfdtridfn: _CallableUFuncWrapper
ncfdtrinc: _CallableUFuncWrapper
nctdtr: _CallableUFuncWrapper
nctdtridf: _CallableUFuncWrapper
nctdtrinc: _CallableUFuncWrapper
nctdtrit: _CallableUFuncWrapper
ndtr: _CallableUFuncWrapper
ndtri: _CallableUFuncWrapper
ndtri_exp: _CallableUFuncWrapper
nrdtrimn: _CallableUFuncWrapper
nrdtrisd: _CallableUFuncWrapper
obl_ang1: _CallableUFuncWrapper
obl_ang1_cv: _CallableUFuncWrapper
obl_cv: _CallableUFuncWrapper
obl_rad1: _CallableUFuncWrapper
obl_rad1_cv: _CallableUFuncWrapper
obl_rad2: _CallableUFuncWrapper
obl_rad2_cv: _CallableUFuncWrapper
owens_t: _CallableUFuncWrapper
pbdv: _CallableUFuncWrapper
pbvv: _CallableUFuncWrapper
pbwa: _CallableUFuncWrapper
pdtr: _CallableUFuncWrapper
pdtrc: _CallableUFuncWrapper
pdtri: _CallableUFuncWrapper
pdtrik: _CallableUFuncWrapper
poch: _CallableUFuncWrapper
poisson_binom_cdf: _CallableUFuncWrapper
polygamma: Callable[..., Any]
powm1: _CallableUFuncWrapper
pro_ang1: _CallableUFuncWrapper
pro_ang1_cv: _CallableUFuncWrapper
pro_cv: _CallableUFuncWrapper
pro_rad1: _CallableUFuncWrapper
pro_rad1_cv: _CallableUFuncWrapper
pro_rad2: _CallableUFuncWrapper
pro_rad2_cv: _CallableUFuncWrapper
pseudo_huber: _CallableUFuncWrapper
psi: _CallableUFuncWrapper
radian: _CallableUFuncWrapper
rel_entr: _CallableUFuncWrapper
rgamma: _CallableUFuncWrapper
round: _CallableUFuncWrapper
shichi: _CallableUFuncWrapper
sici: _CallableUFuncWrapper
sinc: Callable[..., Any]
sindg: _CallableUFuncWrapper
smirnov: _CallableUFuncWrapper
smirnovi: _CallableUFuncWrapper
spence: _CallableUFuncWrapper
stdtr: _CallableUFuncWrapper
stdtridf: _CallableUFuncWrapper
stdtrit: _CallableUFuncWrapper
struve: _CallableUFuncWrapper
tandg: _CallableUFuncWrapper
tklmbda: _CallableUFuncWrapper
voigt_profile: _CallableUFuncWrapper
wofz: _CallableUFuncWrapper
wright_bessel: _CallableUFuncWrapper
wrightomega: _CallableUFuncWrapper
xlog1py: _CallableUFuncWrapper
xlogy: _CallableUFuncWrapper
y0: _CallableUFuncWrapper
y1: _CallableUFuncWrapper
yn: _CallableUFuncWrapper
yv: _CallableUFuncWrapper
yve: _CallableUFuncWrapper
zetac: _CallableUFuncWrapper

# Private wrappers consumed by scipy.stats and internal special code.
_beta_pdf: _CallableUFuncWrapper
_binom_cdf: _CallableUFuncWrapper
_binom_isf: _CallableUFuncWrapper
_binom_pmf: _CallableUFuncWrapper
_binom_ppf: _CallableUFuncWrapper
_binom_sf: _CallableUFuncWrapper
_bivariate_normal_cdf: _CallableUFuncWrapper
_cauchy_isf: _CallableUFuncWrapper
_cauchy_ppf: _CallableUFuncWrapper
_cosine_cdf: _CallableUFuncWrapper
_cosine_invcdf: _CallableUFuncWrapper
_cospi: _CallableUFuncWrapper
_ellip_harm: _CallableUFuncWrapper
_gen_harmonic: _CallableUFuncWrapper
_hypergeom_cdf: _CallableUFuncWrapper
_hypergeom_mean: _CallableUFuncWrapper
_hypergeom_pmf: _CallableUFuncWrapper
_hypergeom_sf: _CallableUFuncWrapper
_hypergeom_skewness: _CallableUFuncWrapper
_hypergeom_variance: _CallableUFuncWrapper
_igam_fac: _CallableUFuncWrapper
_invgauss_isf: _CallableUFuncWrapper
_invgauss_ppf: _CallableUFuncWrapper
_iv_ratio: _CallableUFuncWrapper
_iv_ratio_c: _CallableUFuncWrapper
_iv_ratioinv: _CallableUFuncWrapper
_kolmogc: _CallableUFuncWrapper
_kolmogci: _CallableUFuncWrapper
_kolmogp: _CallableUFuncWrapper
_lambertw: _CallableUFuncWrapper
_landau_cdf: _CallableUFuncWrapper
_landau_isf: _CallableUFuncWrapper
_landau_pdf: _CallableUFuncWrapper
_landau_ppf: _CallableUFuncWrapper
_landau_sf: _CallableUFuncWrapper
_lgam1p: _CallableUFuncWrapper
_log1mexp: _CallableUFuncWrapper
_log1pmx: _CallableUFuncWrapper
_mathieu_ce: _CallableUFuncWrapper
_mathieu_cem: _CallableUFuncWrapper
_mathieu_se: _CallableUFuncWrapper
_mathieu_sem: _CallableUFuncWrapper
_nbinom_cdf: _CallableUFuncWrapper
_nbinom_isf: _CallableUFuncWrapper
_nbinom_kurtosis_excess: _CallableUFuncWrapper
_nbinom_mean: _CallableUFuncWrapper
_nbinom_pmf: _CallableUFuncWrapper
_nbinom_ppf: _CallableUFuncWrapper
_nbinom_sf: _CallableUFuncWrapper
_nbinom_skewness: _CallableUFuncWrapper
_nbinom_variance: _CallableUFuncWrapper
_ncf_isf: _CallableUFuncWrapper
_ncf_kurtosis_excess: _CallableUFuncWrapper
_ncf_mean: _CallableUFuncWrapper
_ncf_pdf: _CallableUFuncWrapper
_ncf_sf: _CallableUFuncWrapper
_ncf_skewness: _CallableUFuncWrapper
_ncf_variance: _CallableUFuncWrapper
_nct_isf: _CallableUFuncWrapper
_nct_kurtosis_excess: _CallableUFuncWrapper
_nct_mean: _CallableUFuncWrapper
_nct_pdf: _CallableUFuncWrapper
_nct_sf: _CallableUFuncWrapper
_nct_skewness: _CallableUFuncWrapper
_nct_variance: _CallableUFuncWrapper
_ncx2_isf: _CallableUFuncWrapper
_ncx2_pdf: _CallableUFuncWrapper
_ncx2_sf: _CallableUFuncWrapper
_normalized_gen_harmonic: _CallableUFuncWrapper
_riemann_zeta: _CallableUFuncWrapper
_scaled_exp1: _CallableUFuncWrapper
_sf_error_test_function: _CallableUFuncWrapper
_sinpi: _CallableUFuncWrapper
_skewnorm_cdf: _CallableUFuncWrapper
_skewnorm_isf: _CallableUFuncWrapper
_skewnorm_ppf: _CallableUFuncWrapper
_smirnovc: _CallableUFuncWrapper
_smirnovci: _CallableUFuncWrapper
_smirnovp: _CallableUFuncWrapper
_spherical_in: _CallableUFuncWrapper
_spherical_in_d: _CallableUFuncWrapper
_spherical_jn: _CallableUFuncWrapper
_spherical_jn_d: _CallableUFuncWrapper
_spherical_kn: _CallableUFuncWrapper
_spherical_kn_d: _CallableUFuncWrapper
_spherical_yn: _CallableUFuncWrapper
_spherical_yn_d: _CallableUFuncWrapper
_stirling2_inexact: _CallableUFuncWrapper
_trivariate_normal_cdf: _CallableUFuncWrapper
_von_mises_cdf: _CallableUFuncWrapper
_zeta: _CallableUFuncWrapper

# Keep the public star-import surface exact.
__all__ = [
    "agm",
    "airy",
    "airye",
    "bdtr",
    "bdtrc",
    "bdtri",
    "bdtrik",
    "bdtrin",
    "bei",
    "beip",
    "ber",
    "berp",
    "besselpoly",
    "beta",
    "betainc",
    "betaincc",
    "betainccinv",
    "betaincinv",
    "betaln",
    "binom",
    "boxcox",
    "boxcox1p",
    "btdtria",
    "btdtrib",
    "cbrt",
    "chdtr",
    "chdtrc",
    "chdtri",
    "chdtriv",
    "chndtr",
    "chndtridf",
    "chndtrinc",
    "chndtrix",
    "cosdg",
    "cosm1",
    "cotdg",
    "dawsn",
    "digamma",
    "digammainv",
    "ellip_harm",
    "ellip_harm_2",
    "ellip_normal",
    "ellipe",
    "ellipeinc",
    "ellipj",
    "ellipk",
    "ellipkinc",
    "ellipkm1",
    "elliprc",
    "elliprd",
    "elliprf",
    "elliprg",
    "elliprj",
    "entr",
    "erf",
    "erfc",
    "erfcinv",
    "erfcx",
    "erfi",
    "erfinv",
    "eval_chebyc",
    "eval_chebys",
    "eval_chebyt",
    "eval_chebyu",
    "eval_gegenbauer",
    "eval_genlaguerre",
    "eval_hermite",
    "eval_hermitenorm",
    "eval_jacobi",
    "eval_laguerre",
    "eval_legendre",
    "eval_sh_chebyt",
    "eval_sh_chebyu",
    "eval_sh_jacobi",
    "eval_sh_legendre",
    "exp1",
    "exp10",
    "exp2",
    "expi",
    "expit",
    "expm1",
    "expn",
    "exprel",
    "fdtr",
    "fdtrc",
    "fdtri",
    "fdtridfd",
    "fresnel",
    "gamma",
    "gammainc",
    "gammaincc",
    "gammainccinv",
    "gammaincinv",
    "gammaln",
    "gammasgn",
    "gdtr",
    "gdtrc",
    "gdtria",
    "gdtrib",
    "gdtrix",
    "hankel1",
    "hankel1e",
    "hankel2",
    "hankel2e",
    "huber",
    "hyp0f1",
    "hyp1f1",
    "hyp2f1",
    "hyperu",
    "i0",
    "i0e",
    "i1",
    "i1e",
    "inv_boxcox",
    "inv_boxcox1p",
    "it2i0k0",
    "it2j0y0",
    "it2struve0",
    "itairy",
    "iti0k0",
    "itj0y0",
    "itmodstruve0",
    "itstruve0",
    "iv",
    "ive",
    "j0",
    "j1",
    "jn",
    "jv",
    "jve",
    "k0",
    "k0e",
    "k1",
    "k1e",
    "kei",
    "keip",
    "kelvin",
    "ker",
    "kerp",
    "kl_div",
    "kn",
    "kolmogi",
    "kolmogorov",
    "kv",
    "kve",
    "lambertw",
    "log1p",
    "log_expit",
    "log_gammainc",
    "log_gammaincc",
    "log_ndtr",
    "log_wright_bessel",
    "loggamma",
    "logit",
    "lpmv",
    "mathieu_a",
    "mathieu_b",
    "mathieu_ce",
    "mathieu_cem",
    "mathieu_modcem1",
    "mathieu_modcem2",
    "mathieu_modsem1",
    "mathieu_modsem2",
    "mathieu_se",
    "mathieu_sem",
    "modfresnelm",
    "modfresnelp",
    "modstruve",
    "multigammaln",
    "nbdtr",
    "nbdtrc",
    "nbdtri",
    "nbdtrik",
    "nbdtrin",
    "ncfdtr",
    "ncfdtri",
    "ncfdtridfd",
    "ncfdtridfn",
    "ncfdtrinc",
    "nctdtr",
    "nctdtridf",
    "nctdtrinc",
    "nctdtrit",
    "ndtr",
    "ndtri",
    "ndtri_exp",
    "nrdtrimn",
    "nrdtrisd",
    "obl_ang1",
    "obl_ang1_cv",
    "obl_cv",
    "obl_rad1",
    "obl_rad1_cv",
    "obl_rad2",
    "obl_rad2_cv",
    "owens_t",
    "pbdv",
    "pbvv",
    "pbwa",
    "pdtr",
    "pdtrc",
    "pdtri",
    "pdtrik",
    "poch",
    "poisson_binom_cdf",
    "polygamma",
    "powm1",
    "pro_ang1",
    "pro_ang1_cv",
    "pro_cv",
    "pro_rad1",
    "pro_rad1_cv",
    "pro_rad2",
    "pro_rad2_cv",
    "pseudo_huber",
    "psi",
    "radian",
    "rel_entr",
    "rgamma",
    "round",
    "shichi",
    "sici",
    "sinc",
    "sindg",
    "smirnov",
    "smirnovi",
    "spence",
    "stdtr",
    "stdtridf",
    "stdtrit",
    "struve",
    "tandg",
    "tklmbda",
    "voigt_profile",
    "wofz",
    "wright_bessel",
    "wrightomega",
    "xlog1py",
    "xlogy",
    "y0",
    "y1",
    "yn",
    "yv",
    "yve",
    "zetac",
]
