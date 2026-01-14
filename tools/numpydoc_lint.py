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
    "SS03",
    "SS05",  # inconsistent standards; see gh-24348
    "SS06",
    "ES01",
    "PR01",
    "PR02",
    "PR03",
    "PR04",
    "PR06",
    "PR07",
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
    "scipy.spatial.cKDTree" : ["SS02"]
}

legacy_functions = [
    "scipy.integrate.complex_ode",
    "scipy.integrate.ode",
    "rv_histogram",
    "rv_continuous",
    "rv_discrete",
    "scipy.interpolate.InterpolatedUnivariateSpline",
    "scipy.interpolate.LSQUnivariateSpline",
    "scipy.interpolate.UnivariateSpline",
    "scipy.sparse.lil_matrix",
    "scipy.sparse.dok_matrix",
    "scipy.sparse.dia_matrix",
    "scipy.sparse.csc_matrix",
    "scipy.sparse.csr_matrix",
    "scipy.sparse.coo_matrix",
    "scipy.sparse.bsr_matrix",
    "scipy.sparse.spmatrix",
    "scipy.optimize.BroydenFirst",
    "scipy.optimize.KrylovJacobian"
]


def walk_class(module_str, class_, public_api):
    class_str = class_.__name__
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
        if ("mstats" in module or "odr" in module or "fftpack" in module or
            "cython" in module):
            # deprecated / legacy modules
            continue
        public_api += walk_module(module)

    errors = 0
    for item in public_api:
        if any(func in item for func in legacy_functions):
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
