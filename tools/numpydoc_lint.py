#!/usr/bin/env python
import sys
import importlib

from numpydoc.validate import validate

from scipy._lib._public_api import PUBLIC_MODULES


skip_errors = [
    "GL01",  # inconsistent standards; see gh-24348
    "GL02",  # inconsistent standards; see gh-24348
    "GL03",  # overlaps with GL02; see gh-24348
    "GL09",
    "SS02",
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


def main():
    # load in numpydoc config
    to_check = []

    # get a list of all public objects
    for module in PUBLIC_MODULES:
        if "odr" in module:
            # skip ODR as deprecated
            continue
        mod = importlib.import_module(module)
        try:
            to_check.extend(obj for f in mod.__all__
                            if (obj := f"{module}.{f}") not in PUBLIC_MODULES)
        except AttributeError:
            # needed for some deprecated modules
            continue

    errors = 0
    for item in to_check:
        try:
            res = validate(item)
        except AttributeError:
            continue
        if res["type"] in ("module", "float", "int", "dict"):
            continue
        for err in res["errors"]:
            if err[0] not in skip_errors:
                print(f"{item}: {err}")
                errors += 1
    sys.exit(errors)


if __name__ == '__main__':
    main()
