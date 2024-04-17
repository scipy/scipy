1#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Check that Python.h is included before any stdlib headers.

May be a bit overzealous, but it should get the job done.
"""
import argparse
import io
import re
import sys

HEADER_PATTERN = re.compile(r'^\s*#\s*include\s*[<"]((?:\w+/)*\w+(?:\.h[hp+]{0,2})?)[>"]\s*$')

PYTHON_INCLUDING_HEADERS = (
  "Python.h",
  "numpy/arrayobject.h",
  "numpy/ndarrayobject.h",
  "numpy/npy_common.h",
  "numpy/npy_math.h",
  "numpy/random/distributions.h",
  "pybind11/pybind11.h"
)

PARSER = argparse.ArgumentParser(description=__doc__)
PARSER.add_argument("file_list", nargs="+", type=str)


def check_python_h_included_first(name_to_check: str) -> int:
    """Check that the passed file includes Python.h first if it does at all.

    Perhaps overzealous, but that should work around concerns with
    recursion.

    Parameters
    ----------
    name_to_check : str
        The name of the file to check.

    Returns
    -------
    int
        The number of headers before Python.h
    """
    included_python = False
    included_other = []
    with open(name_to_check, "r") as in_file:
        for i, line in enumerate(in_file, 1):
            match = re.match(HEADER_PATTERN, line)
            if match:
                if match.group(1) in PYTHON_INCLUDING_HEADERS:
                    if included_other and not included_python:
                        print(
                            f"Header before Python.h in file {name_to_check:s}\n"
                            f"Python.h on line {i:d}, other header(s) on line(s) {included_other}",
                            file=sys.stderr
                        )
                    included_python = True
                elif not included_python and ("numpy" in match.group(1) and match.group(1) != "numpy/utils.h"):
                        print(
                            f"Python.h not included before python-including header "
                            f"in file {name_to_check:s}\n"
                            f"pybind11/pybind11.h on line {i:d}",
                            file=sys.stderr
                        )
                elif not included_python:
                    included_other.append(i)
            elif not included_python and ("py::" in line or "PYBIND11_" in line):
                print(
                    "Python-including header not used before python constructs in file "
                    f"{name_to_check:s}\nConstruct on line {i:d}",
                    file=sys.stderr
                )
    return included_python and len(included_other)


if __name__ == "__main__":
    args = PARSER.parse_args()
    n_out_of_order = 0
    for name_to_check in args.file_list:
        try:
            n_out_of_order += check_python_h_included_first(name_to_check)
        except UnicodeDecodeError as err:
            print(f"File {name_to_check:s} not utf-8", sys.stdout)
    sys.exit(n_out_of_order)
