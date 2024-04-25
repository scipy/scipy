#!/usr/bin/env python
"""Check that Python.h is included before any stdlib headers.

May be a bit overzealous, but it should get the job done.
"""
import argparse
import fnmatch
import os.path
import re
import subprocess
import sys

HEADER_PATTERN = re.compile(
    r'^\s*#\s*include\s*[<"]((?:\w+/)*\w+(?:\.h[hp+]{0,2})?)[>"]\s*$'
)

PYTHON_INCLUDING_HEADERS = [
    "Python.h",
    # This isn't all of Python.h, but it is the visibility macros
    "pyconfig.h",
    "numpy/arrayobject.h",
    "numpy/ndarrayobject.h",
    "numpy/npy_common.h",
    "numpy/npy_math.h",
    "numpy/random/distributions.h",
    "pybind11/pybind11.h",
    # Boost::Python
    "boost/python.hpp",
    "boost/python/args.hpp",
    "boost/python/detail/prefix.hpp",
    "boost/python/detail/wrap_python.hpp",
    "boost/python/ssize_t.hpp",
    "boost/python/object.hpp",
    "boost/mpi/python.hpp",
    # Pythran
    "pythonic/core.hpp",
    # Python-including headers the sort doesn't pick up
    "ni_support.h",
]
LEAF_HEADERS = []

C_CPP_EXTENSIONS = (".c", ".h", ".cpp", ".hpp", ".cc", ".hh", ".cxx", ".hxx")
# check against list in diff_files

PARSER = argparse.ArgumentParser(description=__doc__)
PARSER.add_argument(
    "--diff-against",
    dest="branch",
    type=str,
    default=None,
    help="Diff against "
    "this branch and lint modified files. Use either "
    "`--diff-against` or `--files`, but not both. "
    "Likely to produce false positives.",
)
PARSER.add_argument(
    "files",
    nargs="*",
    help="Lint these files or directories; " "use **/*.py to lint all files",
)


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
    included_non_python_header = []
    warned_python_construct = False
    basename_to_check = os.path.basename(name_to_check)
    in_comment = False
    includes_headers = False
    with open(name_to_check) as in_file:
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
            match = HEADER_PATTERN.match(line)
            if match:
                includes_headers = True
                this_header = match.group(1)
                if this_header in PYTHON_INCLUDING_HEADERS:
                    if included_non_python_header and not included_python:
                        print(
                            f"Header before Python.h in file {name_to_check:s}\n"
                            f"Python.h on line {i:d}, other header(s) on line(s)"
                            f" {included_non_python_header}",
                            file=sys.stderr,
                        )
                    included_python = True
                    PYTHON_INCLUDING_HEADERS.append(basename_to_check)
                elif not included_python and (
                    "numpy" in this_header
                    and this_header != "numpy/utils.h"
                    or "python" in this_header
                ):
                    print(
                        f"Python.h not included before python-including header "
                        f"in file {name_to_check:s}\n"
                        f"{this_header:s} on line {i:d}",
                        file=sys.stderr,
                    )
                elif not included_python and this_header not in LEAF_HEADERS:
                    included_non_python_header.append(i)
            elif (
                not included_python
                and not warned_python_construct
                and ".h" not in basename_to_check
            ) and ("py::" in line or "PYBIND11_" in line or "npy_" in line):
                print(
                    "Python-including header not used before python constructs "
                    f"in file {name_to_check:s}\nConstruct on line {i:d}",
                    file=sys.stderr,
                )
                warned_python_construct = True
    if includes_headers:
        LEAF_HEADERS.append(this_header)
    return included_python and len(included_non_python_header)


def process_files(file_list: list[str]) -> int:
    n_out_of_order = 0
    for name_to_check in sorted(
        file_list, key=lambda name: "h" not in os.path.splitext(name)[1].lower()
    ):
        try:
            n_out_of_order += check_python_h_included_first(name_to_check)
        except UnicodeDecodeError:
            print(f"File {name_to_check:s} not utf-8", sys.stdout)
    return n_out_of_order


def find_c_cpp_files(root: str) -> list[str]:

    result = []

    for dirpath, dirnames, filenames in os.walk("scipy"):
        # I'm assuming other people have checked boost
        for name in ("build", ".git", "boost"):
            try:
                dirnames.remove(name)
            except ValueError:
                pass
        for name in fnmatch.filter(dirnames, "*.p"):
            dirnames.remove(name)
        result.extend(
            [
                os.path.join(dirpath, name)
                for name in filenames
                if os.path.splitext(name)[1].lower() in C_CPP_EXTENSIONS
            ]
        )
    return result


def diff_files(sha: str) -> list[str]:
    """Find the diff since the given SHA.

    Adapted from lint.py
    """
    res = subprocess.run(
        [
            "git",
            "diff",
            "--name-only",
            "--diff-filter=ACMR",
            "-z",
            sha,
            "--",
            # Check against C_CPP_EXTENSIONS
            "*.[chCH]",
            "*.[ch]pp",
            "*.[ch]xx",
            "*.cc",
            "*.hh",
        ],
        stdout=subprocess.PIPE,
        encoding="utf-8",
    )
    res.check_returncode()
    return [f for f in res.stdout.split("\0") if f]


if __name__ == "__main__":
    from lint import find_branch_point

    args = PARSER.parse_args()

    if not ((len(args.files) == 0) ^ (args.branch is None)):
        files = find_c_cpp_files("scipy")
    elif args.branch:
        branch_point = find_branch_point(args.branch)
        files = diff_files(branch_point)
    else:
        files = args.files

    # See which of the headers include Python.h and add them to the list
    n_out_of_order = process_files(files)
    sys.exit(n_out_of_order)
