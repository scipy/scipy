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
import os
import re
import copy
import inspect
import warnings
import doctest
import tempfile
import shutil
from doctest import NORMALIZE_WHITESPACE, ELLIPSIS, IGNORE_EXCEPTION_DETAIL

from argparse import ArgumentParser, REMAINDER

import numpy as np


BASE_MODULE = "scipy"

PUBLIC_SUBMODULES = [
    'linalg',
    'cluster',
    'cluster.vq',
    'cluster.hierarchy',
    'fftpack',
    'interpolate',
    'integrate',
    'io',
    'misc',
    'ndimage',
    'odr',
    'optimize',
    'signal',
    'spatial',
    'sparse',
    'sparse.csgraph',
    'sparse.linalg',
    'stats',
    'stats.mstats',
]

# these names are known to fail doctesting and we like to keep it that way
# e.g. sometimes pseudocode is acceptable etc
DOCTEST_SKIPLIST = set([
    'scipy.integrate.quad',
    'scipy.interpolate.UnivariateSpline',
    'scipy.stats.levy_stable'
])


def short_path(path, cwd=None):
    """
    Return relative or absolute path name, whichever is shortest.
    """
    if not isinstance(path, str):
        return path
    if cwd is None:
        cwd = os.getcwd()
    abspath = os.path.abspath(path)
    relpath = os.path.relpath(path, cwd)
    if len(abspath) <= len(relpath):
        return abspath
    return relpath


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
    if hasattr(module, "__all__"):
        all_dict = copy.deepcopy(module.__all__)
    else:
        all_dict = copy.deepcopy(dir(module))
        all_dict = [name for name in all_dict
                    if not name.startswith("__")]
    for name in ['absolute_import', 'division', 'print_function']:
        try:
            all_dict.remove(name)
        except ValueError:
            pass

    # Modules are almost always private; real submodules need a separate
    # run of refguide_check.
    all_dict = [name for name in all_dict
                if not inspect.ismodule(getattr(module, name, None))]

    deprecated = []
    not_deprecated = []
    for name in all_dict:
        f = getattr(module, name, None)
        if callable(f) and is_deprecated(f):
            deprecated.append(name)
        else:
            not_deprecated.append(name)

    return not_deprecated, deprecated


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

    print("\n\n" + "=" * len(module_name))
    print(module_name)
    print("=" * len(module_name) + "\n")

    num_all = len(all_dict)
    num_ref = len(funcnames)
    print("Non-deprecated objects in __all__: %i" % num_all)
    print("Objects in refguide: %i" % num_ref)

    only_all, only_ref = compare(all_dict, funcnames)
    dep_in_ref = set(only_ref).intersection(deprecated)
    only_ref = set(only_ref).difference(deprecated)
    if len(only_all) == len(only_ref) == 0:
        print("\nNo missing or extraneous items!")
    else:
        if len(only_all) > 0:
            print("")
            print("Objects in %s.__all__ but not in refguide::\n" % module_name)
            for name in only_all:
                print("    " + name)

        if len(only_ref) > 0:
            print("")
            print("Objects in refguide but not in %s.__all__::\n" % module_name)
            for name in only_ref:
                print("    " + name)

        if len(dep_in_ref) > 0:
            print("")
            print("Deprecated objects in refguide::\n")
            for name in deprecated:
                print("    " + name)


def check_docstrings(module, verbose):
    """Check code in docstrings of the module's public symbols.
    """
    # the namespace to run examples in
    ns = {'np': np,
          'assert_allclose': np.testing.assert_allclose,
          'assert_equal': np.testing.assert_equal,
          # recognize numpy repr's
          'array': np.array,
          'int64': np.int64,
          'uint64': np.uint64,
          'int8': np.int8,
          'int32': np.int32,
          'float64': np.float64,
          'dtype': np.dtype,
          'nan': np.nan,
          'NaN': np.nan,
          'inf': np.inf,
          'Inf': np.inf,}

    # if MPL is available, use display-less backend
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        have_matplotlib = True
    except ImportError:
        have_matplotlib = False


    def format_item_header(name):
        return "\n\n" + name + "\n" + "-" * len(name)


    class DTRunner(doctest.DocTestRunner):
        stopwords = {'plt.', '.hist', '.show', '.ylim', '.subplot(',
                     'set_title', 'imshow', 'plt.show', 'ax.axis', 'plt.plot(',
                     '.title', '.ylabel', '.xlabel', 'set_ylim', 'set_xlim'}
        rndm_markers = {'# random', '# Random', '#random', '#Random'}
        DIVIDER = "\n"

        def __init__(self, item_name, checker=None, verbose=None, optionflags=0):
            self._item_name = item_name
            doctest.DocTestRunner.__init__(self, checker=checker, verbose=verbose,
                                           optionflags=optionflags)

        def _report_item_name(self, out, new_line=False):
            if self._item_name is not None:
                out(format_item_header(self._item_name))
                if new_line:
                    out("\n")
                self._item_name = None

        def report_success(self, out, test, example, got):
            if self._verbose:
                self._report_item_name(out, new_line=True)
            return doctest.DocTestRunner.report_success(self, out, test, example, got)

        def report_unexpected_exception(self, out, test, example, exc_info):
            self._report_item_name(out)
            return doctest.DocTestRunner.report_unexpected_exception(
                self, out, test, example, exc_info)

        def report_failure(self, out, test, example, got):
            if (any(word in example.source for word in self.stopwords) or
                any(word in example.want for word in self.rndm_markers)):
                # do not complain if output does not match
                pass
            else:
                self._report_item_name(out)
                return doctest.DocTestRunner.report_failure(self, out, test,
                                                            example, got)

    class Checker(doctest.OutputChecker):
        obj_pattern = re.compile('at 0x[0-9a-fA-F]+>')
        vanilla = doctest.OutputChecker()

        def __init__(self, parse_namedtuples=True, atol=1e-8, rtol=1e-2):
            self.parse_namedtuples = parse_namedtuples
            self.atol, self.rtol = atol, rtol

        def check_output(self, want, got, optionflags):

            # cut it short if they are equal
            if want == got:
                return True

            # skip function/object addresses
            if self.obj_pattern.search(got):
                return True

            # ignore comments (e.g. signal.freqresp)
            if want.lstrip().startswith("#"):
                return True

            # try the standard doctest
            try:
                if self.vanilla.check_output(want, got, optionflags):
                    return True
            except Exception:
                pass

            # OK then, convert strings to objects
            try:
                a_want = eval(want, dict(ns))
                a_got = eval(got, dict(ns))
            except:
                if not self.parse_namedtuples:
                    return False
                # suppose that "want"  is a tuple, and "got" is smth like
                # MoodResult(statistic=10, pvalue=0.1).
                # Then convert the latter to the tuple (10, 0.1), 
                # and then compare the tuples.
                try:
                    num = len(a_want)
                    regex = ('[\w\d_]+\(' +
                             ', '.join(['[\w\d_]+=(.+)']*num) +
                             '\)')
                    grp = re.findall(regex, got.replace('\n', ' '))
                    if len(grp) > 1:  # no more than one for now
                        return False
                    # fold it back to a tuple
                    got_again = '(' + ', '.join(grp[0]) + ')'
                    return self.check_output(want, got_again, optionflags)
                except Exception:
                    return False

            # ... and defer to numpy
            try:
                return self._do_check(a_want, a_got)
            except Exception:
                # heterog tuple, eg (1, np.array([1., 2.]))
               try: 
                    return all(self._do_check(w, g) for w, g in zip(a_want, a_got))
               except TypeError:
                    return False

        def _do_check(self, want, got):
            # This should be done exactly as written to correctly handle all of
            # numpy-comparable objects, strings, and heterogenous tuples
            try:
                if want == got:
                    return True
            except Exception:
                pass
            return np.allclose(want, got, atol=self.atol, rtol=self.rtol)

    # Loop over non-deprecated items
    all_success = True
    for name in get_all_dict(module)[0]:
        full_name = module.__name__ + '.' + name

        if full_name in DOCTEST_SKIPLIST:
            continue

        try:
            obj = getattr(module, name)
        except AttributeError:
            import traceback
            print(format_item_header(full_name))
            print("Missing item!")
            print(traceback.format_exc())
            continue

        finder = doctest.DocTestFinder()
        try:
            tests = finder.find(obj, name, globs=dict(ns))
        except:
            import traceback
            print(format_item_header(full_name))
            print("Failed to get doctests:")
            print(traceback.format_exc())
            continue

        flags = NORMALIZE_WHITESPACE | ELLIPSIS | IGNORE_EXCEPTION_DETAIL
        runner = DTRunner(full_name, checker=Checker(), optionflags=flags,
                          verbose=verbose)

        # Run tests, trying to restore global state afterward
        old_printoptions = np.get_printoptions()
        old_errstate = np.seterr()
        cwd = os.getcwd()
        tmpdir = tempfile.mkdtemp()
        try:
            os.chdir(tmpdir)

            for t in tests:
                t.filename = short_path(t.filename, cwd)
                fails, successes = runner.run(t)
                if fails > 0:
                    all_success = False

            if have_matplotlib:
                plt.close('all')
        finally:
            os.chdir(cwd)
            shutil.rmtree(tmpdir)
            np.set_printoptions(**old_printoptions)
            np.seterr(**old_errstate)

    if not verbose and all_success:
        # Print at least a success message if no other output was produced
        print("\nAll doctests pass!")


def main(argv):
    parser = ArgumentParser(usage=__doc__.lstrip())
    parser.add_argument("module_names", metavar="SUBMODULES", default=list(PUBLIC_SUBMODULES),
                        nargs='*', help="Submodules to check (default: all public)")
    parser.add_argument("--doctests", action="store_true", help="Run also doctests")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    for submodule_name in args.module_names:
        module_name = BASE_MODULE + '.' + submodule_name
        __import__(module_name)
        module = sys.modules[module_name]

        funcnames = find_funcnames(module)
        all_dict, deprecated = get_all_dict(module)
        report(all_dict, funcnames, deprecated, module_name)

        if args.doctests:
            check_docstrings(module, args.verbose)


if __name__ == '__main__':
    main(argv=sys.argv[1:])
