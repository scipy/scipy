#!/usr/bin/env python
import sys
import importlib
import types
import warnings

from numpydoc.validate import validate

from scipy._lib._public_api import PUBLIC_MODULES


skip_errors = [
    "GL01",  # inconsistent standards; see gh-24348
    "GL02",  # inconsistent standards; see gh-24348
    "GL03",  # overlaps with GL02; see gh-24348
    "GL09",
    "SS05",  # inconsistent standards; see gh-24348
    "SS06",
    "ES01",
    "PR01",
    "PR08",
    "PR09",
    "RT02",  # questionable rule; see gh-24348
    "RT04",
    "RT05",
    "SA01",  # questionable rule; see gh-24348
    "SA02",
    "SA03",
    "SA04",
    "EX01",  # remove when gh-7168 is resolved
]


compiled_code_skips = {  # compiled code ignores "numpydoc ignore=" comments"
    "scipy.spatial.cKDTree" : ["SS02"],
    # TODO: the following ufuncs do not have a `__signature__`, this is worth fixing
    # for various typing and introspection reasons. One possible fix is to follow the
    # pattern used in https://github.com/numpy/numpy/pull/30211
    "scipy.special.agm" : ['PR02'],
    "scipy.special.bdtr" : ['PR02'],
    "scipy.special.bdtrc" : ['PR02'],
    "scipy.special.bdtri" : ['PR02'],
    "scipy.special.bdtrik" : ['PR02'],
    "scipy.special.bdtrin" : ['PR02'],
    "scipy.special.betainc" : ['PR02'],
    "scipy.special.betaincc" : ['PR02'],
    "scipy.special.betainccinv" : ['PR02'],
    "scipy.special.betaincinv" : ['PR02'],
    "scipy.special.boxcox" : ['PR02'],
    "scipy.special.boxcox1p" : ['PR02'],
    "scipy.special.btdtria" : ['PR02'],
    "scipy.special.btdtrib" : ['PR02'],
    "scipy.special.chdtr" : ['PR02'],
    "scipy.special.chdtrc" : ['PR02'],
    "scipy.special.chdtri" : ['PR02'],
    "scipy.special.chdtriv" : ['PR02'],
    "scipy.special.chndtr" : ['PR02'],
    "scipy.special.chndtridf" : ['PR02'],
    "scipy.special.chndtrinc" : ['PR02'],
    "scipy.special.chndtrix" : ['PR02'],
    "scipy.special.elliprc" : ['PR02'],
    "scipy.special.elliprd" : ['PR02'],
    "scipy.special.elliprf" : ['PR02'],
    "scipy.special.elliprg" : ['PR02'],
    "scipy.special.elliprj" : ['PR02'],
    "scipy.special.entr" : ['PR02'],
    "scipy.special.erfcinv" : ['PR02'],
    "scipy.special.erfinv" : ['PR02'],
    "scipy.special.eval_chebyc" : ['PR02'],
    "scipy.special.eval_chebys" : ['PR02'],
    "scipy.special.eval_chebyt" : ['PR02'],
    "scipy.special.eval_chebyu" : ['PR02'],
    "scipy.special.eval_gegenbauer" : ['PR02'],
    "scipy.special.eval_genlaguerre" : ['PR02'],
    "scipy.special.eval_hermite" : ['PR02'],
    "scipy.special.eval_hermitenorm" : ['PR02'],
    "scipy.special.eval_jacobi" : ['PR02'],
    "scipy.special.eval_laguerre" : ['PR02'],
    "scipy.special.eval_legendre" : ['PR02'],
    "scipy.special.eval_sh_chebyt" : ['PR02'],
    "scipy.special.eval_sh_chebyu" : ['PR02'],
    "scipy.special.eval_sh_jacobi" : ['PR02'],
    "scipy.special.eval_sh_legendre" : ['PR02'],
    "scipy.special.expn" : ['PR02'],
    "scipy.special.fdtr" : ['PR02'],
    "scipy.special.fdtrc" : ['PR02'],
    "scipy.special.fdtri" : ['PR02'],
    "scipy.special.fdtridfd" : ['PR02'],
    "scipy.special.gdtr" : ['PR02'],
    "scipy.special.gdtrc" : ['PR02'],
    "scipy.special.gdtria" : ['PR02'],
    "scipy.special.gdtrib" : ['PR02'],
    "scipy.special.gdtrix" : ['PR02'],
    "scipy.special.huber" : ['PR02'],
    "scipy.special.hyp0f1" : ['PR02'],
    "scipy.special.hyp1f1" : ['PR02'],
    "scipy.special.hyperu" : ['PR02'],
    "scipy.special.inv_boxcox" : ['PR02'],
    "scipy.special.inv_boxcox1p" : ['PR02'],
    "scipy.special.kl_div" : ['PR02'],
    "scipy.special.kn" : ['PR02'],
    "scipy.special.kolmogi" : ['PR02'],
    "scipy.special.kolmogorov" : ['PR02'],
    "scipy.special.lpmv" : ['PR02'],
    "scipy.special.nbdtr" : ['PR02'],
    "scipy.special.nbdtrc" : ['PR02'],
    "scipy.special.nbdtri" : ['PR02'],
    "scipy.special.nbdtrik" : ['PR02'],
    "scipy.special.nbdtrin" : ['PR02'],
    "scipy.special.ncfdtr" : ['PR02'],
    "scipy.special.ncfdtri" : ['PR02'],
    "scipy.special.ncfdtridfd" : ['PR02'],
    "scipy.special.ncfdtridfn" : ['PR02'],
    "scipy.special.ncfdtrinc" : ['PR02'],
    "scipy.special.nctdtr" : ['PR02'],
    "scipy.special.nctdtridf" : ['PR02'],
    "scipy.special.nctdtrinc" : ['PR02'],
    "scipy.special.nctdtrit" : ['PR02'],
    "scipy.special.ndtri" : ['PR02'],
    "scipy.special.ndtri_exp" : ['PR02'],
    "scipy.special.nrdtrimn" : ['PR02'],
    "scipy.special.nrdtrisd" : ['PR02'],
    "scipy.special.owens_t" : ['PR02'],
    "scipy.special.pdtr" : ['PR02'],
    "scipy.special.pdtrc" : ['PR02'],
    "scipy.special.pdtri" : ['PR02'],
    "scipy.special.pdtrik" : ['PR02'],
    "scipy.special.poch" : ['PR02'],
    "scipy.special.powm1" : ['PR02'],
    "scipy.special.pseudo_huber" : ['PR02'],
    "scipy.special.rel_entr" : ['PR02'],
    "scipy.special.round" : ['PR02'],
    "scipy.special.shichi" : ['PR02'],
    "scipy.special.sici" : ['PR02'],
    "scipy.special.smirnov" : ['PR02'],
    "scipy.special.smirnovi" : ['PR02'],
    "scipy.special.spence" : ['PR02'],
    "scipy.special.stdtr" : ['PR02'],
    "scipy.special.stdtridf" : ['PR02'],
    "scipy.special.stdtrit" : ['PR02'],
    "scipy.special.tklmbda" : ['PR02'],
    "scipy.special.wrightomega" : ['PR02'],
    "scipy.special.yn" : ['PR02'],
    "scipy.special.jn" : ['PR02'],
    "scipy.special.airy" : ['PR02'],
    "scipy.special.airye" : ['PR02'],
    "scipy.special.bei" : ['PR02'],
    "scipy.special.beip" : ['PR02'],
    "scipy.special.ber" : ['PR02'],
    "scipy.special.berp" : ['PR02'],
    "scipy.special.binom" : ['PR02'],
    "scipy.special.exp1" : ['PR02'],
    "scipy.special.expi" : ['PR02'],
    "scipy.special.expit" : ['PR02'],
    "scipy.special.exprel" : ['PR02'],
    "scipy.special.gamma" : ['PR02'],
    "scipy.special.gammaln" : ['PR02'],
    "scipy.special.hankel1" : ['PR02'],
    "scipy.special.hankel1e" : ['PR02'],
    "scipy.special.hankel2" : ['PR02'],
    "scipy.special.hankel2e" : ['PR02'],
    "scipy.special.hyp2f1" : ['PR02'],
    "scipy.special.it2i0k0" : ['PR02'],
    "scipy.special.it2j0y0" : ['PR02'],
    "scipy.special.it2struve0" : ['PR02'],
    "scipy.special.itairy" : ['PR02'],
    "scipy.special.iti0k0" : ['PR02'],
    "scipy.special.itj0y0" : ['PR02'],
    "scipy.special.itmodstruve0" : ['PR02'],
    "scipy.special.itstruve0" : ['PR02'],
    "scipy.special.iv" : ['PR02'],
    "scipy.special.ive" : ['PR02'],
    "scipy.special.jv" : ['PR02'],
    "scipy.special.jve" : ['PR02'],
    "scipy.special.kei" : ['PR02'],
    "scipy.special.keip" : ['PR02'],
    "scipy.special.kelvin" : ['PR02'],
    "scipy.special.ker" : ['PR02'],
    "scipy.special.kerp" : ['PR02'],
    "scipy.special.kv" : ['PR02'],
    "scipy.special.kve" : ['PR02'],
    "scipy.special.log_expit" : ['PR02'],
    "scipy.special.log_wright_bessel" : ['PR02'],
    "scipy.special.loggamma" : ['PR02'],
    "scipy.special.logit" : ['PR02'],
    "scipy.special.mathieu_a" : ['PR02'],
    "scipy.special.mathieu_b" : ['PR02'],
    "scipy.special.mathieu_cem" : ['PR02'],
    "scipy.special.mathieu_modcem1" : ['PR02'],
    "scipy.special.mathieu_modcem2" : ['PR02'],
    "scipy.special.mathieu_modsem1" : ['PR02'],
    "scipy.special.mathieu_modsem2" : ['PR02'],
    "scipy.special.mathieu_sem" : ['PR02'],
    "scipy.special.modfresnelm" : ['PR02'],
    "scipy.special.modfresnelp" : ['PR02'],
    "scipy.special.obl_ang1" : ['PR02'],
    "scipy.special.obl_ang1_cv" : ['PR02'],
    "scipy.special.obl_cv" : ['PR02'],
    "scipy.special.obl_rad1" : ['PR02'],
    "scipy.special.obl_rad1_cv" : ['PR02'],
    "scipy.special.obl_rad2" : ['PR02'],
    "scipy.special.obl_rad2_cv" : ['PR02'],
    "scipy.special.pbdv" : ['PR02'],
    "scipy.special.pbvv" : ['PR02'],
    "scipy.special.pbwa" : ['PR02'],
    "scipy.special.pro_ang1" : ['PR02'],
    "scipy.special.pro_ang1_cv" : ['PR02'],
    "scipy.special.pro_cv" : ['PR02'],
    "scipy.special.pro_rad1" : ['PR02'],
    "scipy.special.pro_rad1_cv" : ['PR02'],
    "scipy.special.pro_rad2" : ['PR02'],
    "scipy.special.pro_rad2_cv" : ['PR02'],
    "scipy.special.psi" : ['PR02'],
    "scipy.special.rgamma" : ['PR02'],
    "scipy.special.wright_bessel" : ['PR02'],
    "scipy.special.yv" : ['PR02'],
    "scipy.special.yve" : ['PR02'],
    "scipy.special.zetac" : ['PR02'],
    "scipy.special.sindg" : ['PR02'],
    "scipy.special.cosdg" : ['PR02'],
    "scipy.special.tandg" : ['PR02'],
    "scipy.special.cotdg" : ['PR02'],
    "scipy.special.i0" : ['PR02'],
    "scipy.special.i0e" : ['PR02'],
    "scipy.special.i1" : ['PR02'],
    "scipy.special.i1e" : ['PR02'],
    "scipy.special.k0" : ['PR02'],
    "scipy.special.k0e" : ['PR02'],
    "scipy.special.k1" : ['PR02'],
    "scipy.special.k1e" : ['PR02'],
    "scipy.special.y0" : ['PR02'],
    "scipy.special.y1" : ['PR02'],
    "scipy.special.j0" : ['PR02'],
    "scipy.special.j1" : ['PR02'],
    "scipy.special.struve" : ['PR02'],
    "scipy.special.modstruve" : ['PR02'],
    "scipy.special.beta" : ['PR02'],
    "scipy.special.betaln" : ['PR02'],
    "scipy.special.besselpoly" : ['PR02'],
    "scipy.special.gammasgn" : ['PR02'],
    "scipy.special.cbrt" : ['PR02'],
    "scipy.special.radian" : ['PR02'],
    "scipy.special.cosm1" : ['PR02'],
    "scipy.special.gammainc" : ['PR02'],
    "scipy.special.gammaincinv" : ['PR02'],
    "scipy.special.gammaincc" : ['PR02'],
    "scipy.special.gammainccinv" : ['PR02'],
    "scipy.special.fresnel" : ['PR02'],
    "scipy.special.ellipe" : ['PR02'],
    "scipy.special.ellipeinc" : ['PR02'],
    "scipy.special.ellipk" : ['PR02'],
    "scipy.special.ellipkinc" : ['PR02'],
    "scipy.special.ellipkm1" : ['PR02'],
    "scipy.special.ellipj" : ['PR02'],
    "scipy.special.erf" : ['PR02'],
    "scipy.special.erfc" : ['PR02'],
    "scipy.special.erfcx" : ['PR02'],
    "scipy.special.erfi" : ['PR02'],
    "scipy.special.voigt_profile" : ['PR02'],
    "scipy.special.wofz" : ['PR02'],
    "scipy.special.dawsn" : ['PR02'],
    "scipy.special.ndtr" : ['PR02'],
    "scipy.special.log_ndtr" : ['PR02'],
    "scipy.special.exp2" : ['PR02'],
    "scipy.special.exp10" : ['PR02'],
    "scipy.special.expm1" : ['PR02'],
    "scipy.special.log1p" : ['PR02'],
    "scipy.special.xlogy" : ['PR02'],
    "scipy.special.xlog1py" : ['PR02'],
    "scipy.special.digamma" : ['PR02'],
    "scipy.special.assoc_legendre_p" : ['PR02'],
    "scipy.special.legendre_p" : ['PR02'],
    "scipy.special.sph_harm_y" : ['PR02'],
    "scipy.special.sph_legendre_p" : ['PR02'],
}

legacy_functions = [
    "scipy.integrate.complex_ode",
    "scipy.integrate.ode",
    "scipy.stats.rv_histogram",
    "scipy.stats.distributions.rv_histogram",
    "scipy.stats.rv_continuous",
    "scipy.stats.distributions.rv_continuous",
    "scipy.stats.rv_discrete",
    "scipy.stats.distributions.rv_discrete",
    "scipy.interpolate.InterpolatedUnivariateSpline",
    "scipy.interpolate.LSQUnivariateSpline",
    "scipy.interpolate.UnivariateSpline",
    "scipy.interpolate.splder",
    "scipy.interpolate.Rbf",
    "scipy.sparse.lil_matrix",
    "scipy.sparse.dok_matrix",
    "scipy.sparse.dia_matrix",
    "scipy.sparse.csc_matrix",
    "scipy.sparse.csr_matrix",
    "scipy.sparse.coo_matrix",
    "scipy.sparse.bsr_matrix",
    "scipy.sparse.isspmatrix_lil",
    "scipy.sparse.isspmatrix_dok",
    "scipy.sparse.isspmatrix_dia",
    "scipy.sparse.isspmatrix_csr",
    "scipy.sparse.isspmatrix_csc",
    "scipy.sparse.isspmatrix_coo",
    "scipy.sparse.isspmatrix_bsr",
    "scipy.sparse.isspmatrix",
    "scipy.sparse.spmatrix",
    "scipy.optimize.BroydenFirst",
    "scipy.optimize.KrylovJacobian",
    "scipy.optimize.newton_krylov",
    "scipy.linalg.LinAlgError",  # this is from numpy
    "scipy.optimize.fmin_bfgs",
    "scipy.optimize.fmin_tnc",
    "scipy.optimize.fmin_ncg",
    "scipy.optimize.fmin_cobyla"
]

# the method of these classes have no __doc__, skip for now
false_positives = ["scipy.stats.Uniform",
                   "scipy.stats.Normal",
                   "scipy.stats.Mixture",
                   "scipy.stats.Binomial",
                   "scipy.stats.Logistic"]

skip_modules = [
    "scipy.odr",
    "scipy.fftpack",
    "scipy.stats.mstats",
    "scipy.linalg.cython_lapack",
    "scipy.linalg.cython_blas",
]


def walk_class(module_str, class_, public_api):
    class_str = class_.__name__
    # skip private methods (and dunder methods)
    attrs = {a for a in dir(class_) if not a.startswith("_")}
    for attr in attrs:
        item = getattr(class_, attr)
        if isinstance(item, types.FunctionType):
            public_api.append(f"{module_str}.{class_str}.{attr}")


def walk_module(module_str):
    public_api = []

    module = importlib.import_module(module_str)

    for item_str in module.__all__:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", category=DeprecationWarning)
            item = getattr(module, item_str)
            if w:
                continue

        if isinstance(item, float | int | dict):  # ignore constants
            continue
        elif isinstance(item, types.ModuleType):
            continue
        elif isinstance(item, type):  # classes
            public_api.append(f"{module_str}.{item_str}")
            walk_class(module_str, item, public_api)  # methods
        else:  # functions
            public_api.append(f"{module_str}.{item_str}")
    return public_api


def main():
    public_api = []

    # get a list of all public objects
    for module in PUBLIC_MODULES:
        if module in skip_modules:
            # deprecated / legacy modules
            continue
        public_api += walk_module(module)

    errors = 0
    for item in public_api:
        if (any(func in item for func in legacy_functions) or
            any(func in item for func in false_positives)):
            continue
        try:
            res = validate(item)
        except AttributeError:
            continue
        for err in res["errors"]:
            if (err[0] not in skip_errors and
                err[0] not in compiled_code_skips.get(item, [])):
                print(f"{item}: {err}")
                errors += 1
    sys.exit(errors)


if __name__ == '__main__':
    main()
