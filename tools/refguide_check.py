#!/usr/bin/env python
"""
refguide_check.py [OPTIONS] [-- ARGS]

Check for a Scipy submodule whether the objects in its __all__ dict
correspond to the objects included in the reference guide.

Example of usage::

    $ python refguide_check.py optimize

Note that this is a helper script to be able to check if things are missing;
the output of this script does need to be checked manually.  In some cases
objects are left out of the refguide for a good reason (it's an alias of
another function, or deprecated, or ...)

"""

import sys
import re
import copy
import inspect

from argparse import ArgumentParser, REMAINDER

import scipy
from scipy import (cluster, constants, fftpack, integrate, interpolate, io,
                   linalg, misc, ndimage, odr, optimize, signal, sparse,
                   spatial, special, stats)
# TODO: sparse.csgraph, sparse.linalg, stats.mstats, cluster.vq,
#       cluster.hierarchy


def find_funcnames(module):
    funcnames = set()
    # 3 spaces followed by function name; only function names listed in
    # refguide are indented like this (mostly, there may be some false
    # positives)
    pattern = re.compile("(\s\s\s[a-z_0-9A-Z]+)")
    for line in module.__doc__.splitlines():
        res = re.search(pattern, line)
        if res is not None:
            funcname = res.groups()[0].lstrip()
            funcnames.add(funcname)

    return funcnames


def get_all_dict(module):
    """Return a copy of the __all__ dict with irrelevant items removed."""
    all = copy.deepcopy(module.__all__)
    for name in ['absolute_import', 'division', 'print_function']:
        try:
            all.remove(name)
        except ValueError:
            pass

    # somehow some modules survive the first iteration (?)
    for _ in range(2):
        for name in all:
            if inspect.ismodule(getattr(module, name)):
                all.remove(name)

    return all


def compare(all, funcnames):
    """Return sets of objects only in one of __all__, refguide."""
    only_all = set()
    for name in all:
        if name not in funcnames:
            only_all.add(name)

    only_ref = set()
    for name in funcnames:
        if name not in all:
            only_ref.add(name)

    return only_all, only_ref


def report(all, funcnames, module_name):
    """Print out a report for the module"""
    num_all = len(all)
    num_ref = len(funcnames)
    print("Number of functions in __all__: %i" % num_all)
    print("Number of functions in refguide: %i" % num_ref)

    only_all, only_ref = compare(all, funcnames)
    if len(only_all) == len(only_ref) == 0:
        print("\nAll good!")
    else:
        if len(only_all) > 0:
            print("")
            print("Objects in %s.__all__ but not in refguide:" % module_name)
            print("------------------------------------------")
            for name in only_all:
                print(name)

        if len(only_ref) > 0:
            print("")
            print("Objects in refguide but not in %s.__all__:" % module_name)
            print("------------------------------------------")
            for name in only_ref:
                print(name)


def main(argv):
    parser = ArgumentParser(usage=__doc__.lstrip())
    parser.add_argument("module_name", metavar="ARGS", default=[],
                        nargs=REMAINDER, help="Valid Scipy submodule name")
    args = parser.parse_args(argv)

    module_name = args.module_name[0]
    module = getattr(scipy, module_name)

    funcnames = find_funcnames(module)
    all = get_all_dict(module)
    report(all, funcnames, module_name)


if __name__ == '__main__':
    main(argv=sys.argv[1:])
