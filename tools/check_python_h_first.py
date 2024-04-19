1#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Check that Python.h is included before any stdlib headers.

May be a bit overzealous, but it should get the job done.
"""
import argparse
import io
import os.path
import re
import sys

HEADER_PATTERN = re.compile(r'^\s*#\s*include\s*[<"]((?:\w+/)*\w+(?:\.h[hp+]{0,2})?)[>"]\s*$')

PYTHON_INCLUDING_HEADERS = [
  "Python.h",
  "numpy/arrayobject.h",
  "numpy/ndarrayobject.h",
  "numpy/npy_common.h",
  "numpy/npy_math.h",
  "numpy/random/distributions.h",
  "pybind11/pybind11.h",
  # Python-including headers the sort doesn't pick up
  "ni_support.h",
]

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
    included_non_python_headers = []
    included_headers = []
    duplicate_header_count = 0
    warned_python_construct = False
    basename_to_check = os.path.basename(name_to_check)
    in_comment = False
    with open(name_to_check, "r") as in_file:
        for i, line in enumerate(in_file, 1):
            # Very basic comment parsing
            # Assumes /*...*/ comments are on their own lines
            if "/*" in line:
                if "*/" not in line:
                    in_comment = True
                # else-branch could use regex to remove comment and continue
                continue
            if in_comment:
                if "*/" in line:
                    in_comment = False
                continue
            match = re.match(HEADER_PATTERN, line)
            if match:
                this_header = match.group(1)
                if this_header in PYTHON_INCLUDING_HEADERS:
                    if included_non_python_headers and not included_python:
                        print(
                            f"Header before Python.h in file {name_to_check:s}\n"
                            f"Python.h on line {i:d}, other header(s) on line(s) {included_non_python_headers}",
                            file=sys.stderr
                        )
                    included_python = True
                    PYTHON_INCLUDING_HEADERS.append(basename_to_check)
                elif not included_python and ("numpy" in this_header and this_header != "numpy/utils.h"):
                        print(
                            f"Python.h not included before python-including header "
                            f"in file {name_to_check:s}\n"
                            f"{this_header:s} on line {i:d}",
                            file=sys.stderr
                        )
                elif not included_python:
                    included_non_python_headers.append(i)

                if this_header in included_headers:
                    print(
                        f"Header {this_header:s} included twice in file {name_to_check:s}, second on line {i:d}",
                        file=sys.stderr
                    )
                    duplicate_header_count += 1
                included_headers.append(this_header)
            elif (not included_python and not warned_python_construct and ".h" not in basename_to_check) and (
                "py::" in line or "PYBIND11_" in line or "npy_" in line
            ):
                print(
                    "Python-including header not used before python constructs in file "
                    f"{name_to_check:s}\nConstruct on line {i:d}",
                    file=sys.stderr
                )
                warned_python_construct = True
    return (included_python and len(included_non_python_headers)) + duplicate_header_count


if __name__ == "__main__":
    args = PARSER.parse_args()
    n_out_of_order = 0
    # See which of the headers include Python.h and add them to the list
    for name_to_check in sorted(args.file_list, key=lambda name: "h" not in name.lower()):
        try:
            n_out_of_order += check_python_h_included_first(name_to_check)
        except UnicodeDecodeError as err:
            print(f"File {name_to_check:s} not utf-8", sys.stdout)
    sys.exit(n_out_of_order)
