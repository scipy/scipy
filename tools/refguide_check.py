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

Another use of this helper script is to check validity of code samples
in docstrings. This is different from doctesting [we do not aim to have
scipy docstrings doctestable!], this is just to make sure that code in
docstrings is valid python::

    $ python refguide_check.py --check_docs optimize

"""
from __future__ import print_function

import sys
import re
import copy
import inspect
import doctest
import warnings

from argparse import ArgumentParser, REMAINDER

import numpy as np

import scipy
from scipy import (cluster, constants, fftpack, integrate, interpolate, io,
                   linalg, misc, ndimage, odr, optimize, signal, sparse,
                   spatial, special, stats)
# TODO: sparse.csgraph, sparse.linalg, stats.mstats, cluster.vq,
#       cluster.hierarchy


def find_funcnames(module):
    funcnames = set()
    # 3 spaces followed by function name, and maybe some spaces, some
    # dashes, and an explanation; only function names listed in
    # refguide are formatted like this (mostly, there may be some false
    # positives)
    pattern = re.compile("^\s\s\s([a-z_0-9A-Z]+)(\s+-+.*)?$")
    for line in module.__doc__.splitlines():
        res = re.search(pattern, line)
        if res is not None:
            funcname = res.groups()[0]
            funcnames.add(funcname)

    return funcnames


def get_all_dict(module):
    """Return a copy of the __all__ dict with irrelevant items removed."""
    all_dict = copy.deepcopy(module.__all__)
    for name in ['absolute_import', 'division', 'print_function']:
        try:
            all_dict.remove(name)
        except ValueError:
            pass

    # FIXME: shouldn't modules be in the refguide, actually? if they're in __all__?
    # somehow some modules survive the first iteration (?)
    for _ in range(2):
        for name in all_dict:
            if inspect.ismodule(getattr(module, name)):
                all_dict.remove(name)

    deprecated = []
    for name in all_dict:
        f = getattr(module, name)
        if callable(f) and is_deprecated(f):
            all_dict.remove(name)
            deprecated.append(name)
            
    return all_dict, deprecated


def compare(all_dict, funcnames):
    """Return sets of objects only in one of __all__, refguide."""
    only_all = set()
    for name in all_dict:
        if name not in funcnames:
            only_all.add(name)

    only_ref = set()
    for name in funcnames:
        if name not in all_dict:
            only_ref.add(name)

    return only_all, only_ref

def is_deprecated(f):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("error")
        try:
            f(**{"not a kwarg":None})
        except DeprecationWarning:
            return True
        except:
            pass
        return False

def report(all_dict, funcnames, deprecated, module_name):
    """Print out a report for the module"""
    num_all = len(all_dict)
    num_ref = len(funcnames)
    print("Number of non-deprecated functions in __all__: %i" % num_all)
    print("Number of functions in refguide: %i" % num_ref)

    only_all, only_ref = compare(all_dict, funcnames)
    dep_in_ref = set(only_ref).intersection(deprecated)
    only_ref = set(only_ref).difference(deprecated)
    if len(only_all) == len(only_ref) == 0:
        print("\nAll good!")
    else:
        if len(only_all) > 0:
            print("")
            print("Functions in %s.__all__ but not in refguide:" % module_name)
            print("------------------------------------------")
            for name in only_all:
                print(name)

        if len(only_ref) > 0:
            print("")
            print("Objects in refguide but not functions in %s.__all__:" % module_name)
            print("------------------------------------------")
            for name in only_ref:
                print(name)

        if len(dep_in_ref) > 0:
            print("")
            print("Deprecated objects in refguide:")
            print("------------------------------------------")
            for name in deprecated:
                print(name)


def check_docstrings(module):
    """Check the code in the docstrings of the module's public symbols.
    """

    class DTRunner(doctest.DocTestRunner):
        def report_failure(self, out, test, example, got):
            # do not complain if output does not match
            pass

    # namespace to run examples in
    ns = {'np': np,
          'assert_allclose': np.testing.assert_allclose,
          'assert_equal': np.testing.assert_equal}

    for name in get_all_dict(module):
        obj = getattr(module, name)

        finder = doctest.DocTestFinder()
        tests = finder.find(obj, name, globs=ns)

        print(name)

        runner = DTRunner()
        for t in tests:
            # do not show MPL figures
            for j, ex in enumerate(t.examples):
                if 'show()' in ex.source:
                    t.examples[j].source = ex.source.replace('show()', 'show')
            runner.run(t)


def main(argv):
    parser = ArgumentParser(usage=__doc__.lstrip())
    parser.add_argument("module_name", metavar="ARGS", default=[],
                        nargs=REMAINDER, help="Valid Scipy submodule name")
    parser.add_argument("--check_docs", action="store_true")
    args = parser.parse_args(argv)

    module_name = args.module_name[0]
    module = getattr(scipy, module_name)

    if args.check_docs:
        check_docstrings(module)
    else:
        funcnames = find_funcnames(module)
        all_dict, deprecated = get_all_dict(module)
        report(all_dict, funcnames, deprecated, module_name)


if __name__ == '__main__':
    main(argv=sys.argv[1:])
