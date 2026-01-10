#!/usr/bin/env python
import sys
import importlib
import os

from numpydoc.validate import validate
from numpydoc.hooks.validate_docstrings import parse_config


from scipy._lib._public_api import PUBLIC_MODULES

def main():
    root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    config = parse_config(root_path)
    error = 0
    too_check = []
    for module in PUBLIC_MODULES:
        mod = importlib.import_module(module)
        try:
            too_check.extend([obj for f in mod.__all__ 
                              if (obj := f"{module}.{f}") not in PUBLIC_MODULES])
        except AttributeError:
            continue
    for item in too_check:
        try:
            res = validate(item)
        except AttributeError:
            continue
        if res["type"] == "module":
            continue
        for err in res["errors"]:
            if err[0] in config["checks"]:
                print(f"{res["file"]}: {res["file_line"]}: {err}")
                error += 1
    sys.exit(error)


if __name__ == '__main__':
    main()
