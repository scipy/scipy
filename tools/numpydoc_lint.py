#!/usr/bin/env python
import sys
import importlib
import os

from numpydoc.validate import validate
from numpydoc.hooks.validate_docstrings import parse_config


from scipy._lib._public_api import PUBLIC_MODULES

def main():
    root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

    # load in numpydoc config
    config = parse_config(root_path)
    to_check = []

    # get a list of all public objects
    for module in PUBLIC_MODULES:
        mod = importlib.import_module(module)
        try:
            to_check.extend([obj for f in mod.__all__
                             if (obj := f"{module}.{f}") not in PUBLIC_MODULES])
        except AttributeError:
            # needed for some deprecated modules
            continue

    errors = 0
    for item in to_check:
        try:
            res = validate(item)
        except AttributeError:
            continue
        if res["type"] == "module":
            continue
        for err in res["errors"]:
            if err[0] in config["checks"]:
                print(f"{item}: {err}")
                errors += 1
    sys.exit(errors)


if __name__ == '__main__':
    main()
