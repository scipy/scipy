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
import io
import docutils.core
from docutils.parsers.rst import directives
import shutil
from doctest import NORMALIZE_WHITESPACE, ELLIPSIS, IGNORE_EXCEPTION_DETAIL
from argparse import ArgumentParser, REMAINDER
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'doc', 'sphinxext'))
from numpydoc.docscrape_sphinx import get_doc_object
# Remove sphinx directives that don't run without Sphinx environment
directives._directives.pop('versionadded', None)
directives._directives.pop('moduleauthor', None)
directives._directives.pop('sectionauthor', None)
directives._directives.pop('codeauthor', None)
directives._directives.pop('toctree', None)


BASE_MODULE = "scipy"

PUBLIC_SUBMODULES = [
    'cluster',
    'cluster.hierarchy',
    'cluster.vq',
    'constants',
    'fftpack',
    'fftpack.convolve',
    'integrate',
    'interpolate',
    'io',
    'io.arff',
    'io.wavfile',
    'linalg',
    'linalg.blas',
    'linalg.lapack',
    'linalg.interpolative',
    'misc',
    'ndimage',
    'odr',
    'optimize',
    'signal',
    'sparse',
    'sparse.csgraph',
    'sparse.linalg',
    'spatial',
    'spatial.distance',
    'special',
    'stats',
    'stats.mstats',
]

# Docs for these modules are included in the parent module
OTHER_MODULE_DOCS = {
    'fftpack.convolve': 'fftpack',
    'io.wavfile': 'io',
    'io.arff': 'io',
}

# these names are known to fail doctesting and we like to keep it that way
# e.g. sometimes pseudocode is acceptable etc
DOCTEST_SKIPLIST = set([
    'scipy.stats.kstwobign', # inaccurate cdf or ppf
    'scipy.stats.levy_stable',
    'scipy.special.sinc', # comes from numpy
    'scipy.misc.who', # comes from numpy
])

# these names are not required to be present in ALL despite being in
# autosummary:: listing
REFGUIDE_ALL_SKIPLIST = [
    r'scipy\.sparse\.csgraph',
    r'scipy\.sparse\.linalg',
    r'scipy\.spatial\.distance',
    r'scipy\.linalg\.blas\.[sdczi].*',
    r'scipy\.linalg\.lapack\.[sdczi].*',
]


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


def find_names(module, names_dict):
    # Refguide entries:
    #
    # - 3 spaces followed by function name, and maybe some spaces, some
    #   dashes, and an explanation; only function names listed in
    #   refguide are formatted like this (mostly, there may be some false
    #   positives)
    #
    # - special directives, such as data and function
    #
    # - (scipy.constants only): quoted list
    #
    patterns = [
        r"^\s\s\s([a-z_0-9A-Z]+)(\s+-+.*)?$",
        r"^\.\. (?:data|function)::\s*([a-z_0-9A-Z]+)\s*$"
    ]

    if module.__name__ == 'scipy.constants':
        patterns += ["^``([a-z_0-9A-Z]+)``"]

    patterns = [re.compile(pattern) for pattern in patterns]
    module_name = module.__name__

    for line in module.__doc__.splitlines():
        res = re.search(r"^\s*\.\. (?:currentmodule|module):: ([a-z0-9A-Z_.]+)\s*$", line)
        if res:
            module_name = res.group(1)
            continue

        for pattern in patterns:
            res = re.match(pattern, line)
            if res is not None:
                name = res.group(1)
                entry = '.'.join([module_name, name])
                names_dict.setdefault(module_name, set()).add(name)
                break


def get_all_dict(module):
    """Return a copy of the __all__ dict with irrelevant items removed."""
    if hasattr(module, "__all__"):
        all_dict = copy.deepcopy(module.__all__)
    else:
        all_dict = copy.deepcopy(dir(module))
        all_dict = [name for name in all_dict
                    if not name.startswith("_")]
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

    others = set(dir(module)).difference(set(deprecated)).difference(set(not_deprecated))

    return not_deprecated, deprecated, others


def compare(all_dict, others, names, module_name):
    """Return sets of objects only in __all__, refguide, or completely missing."""
    only_all = set()
    for name in all_dict:
        if name not in names:
            only_all.add(name)

    only_ref = set()
    missing = set()
    for name in names:
        if name not in all_dict:
            for pat in REFGUIDE_ALL_SKIPLIST:
                if re.match(pat, module_name + '.' + name):
                    if name not in others:
                        missing.add(name)
                    break
            else:
                only_ref.add(name)

    return only_all, only_ref, missing

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

def check_items(all_dict, names, deprecated, others, module_name, dots=True):
    num_all = len(all_dict)
    num_ref = len(names)

    output = ""

    output += "Non-deprecated objects in __all__: %i\n" % num_all
    output += "Objects in refguide: %i\n\n" % num_ref

    only_all, only_ref, missing = compare(all_dict, others, names, module_name)
    dep_in_ref = set(only_ref).intersection(deprecated)
    only_ref = set(only_ref).difference(deprecated)

    if len(dep_in_ref) > 0:
        output += "Deprecated objects in refguide::\n\n"
        for name in sorted(deprecated):
            output += "    " + name + "\n"

    if len(only_all) == len(only_ref) == len(missing) == 0:
        if dots:
            output_dot('.')
        return [(None, True, output)]
    else:
        if len(only_all) > 0:
            output += "ERROR: objects in %s.__all__ but not in refguide::\n\n" % module_name
            for name in sorted(only_all):
                output += "    " + name + "\n"

        if len(only_ref) > 0:
            output += "ERROR: objects in refguide but not in %s.__all__::\n\n" % module_name
            for name in sorted(only_ref):
                output += "    " + name + "\n"

        if len(missing) > 0:
            output += "ERROR: missing objects::\n\n"
            for name in sorted(missing):
                output += "    " + name + "\n"

        if dots:
            output_dot('F')
        return [(None, False, output)]


def validate_rst_syntax(text, name, dots=True):
    if text is None:
        if dots:
            output_dot('E')
        return False, "ERROR: %s: no documentation" % (name,)

    ok_unknown_items = set([
        'mod', 'currentmodule', 'autosummary', 'data',
        'obj', 'versionadded', 'module', 'class',
        'ref', 'func', 'toctree', 'moduleauthor',
        'sectionauthor', 'codeauthor',
    ])

    # Run through docutils
    error_stream = io.StringIO()

    def resolve(name, is_label=False):
        return ("http://foo", name)

    token = '<RST-VALIDATE-SYNTAX-CHECK>'

    docutils.core.publish_doctree(
        text, token,
        settings_overrides = dict(halt_level=5,
                                  traceback=True,
                                  default_reference_context='title-reference',
                                  default_role='emphasis',
                                  link_base='',
                                  resolve_name=resolve,
                                  stylesheet_path='',
                                  raw_enabled=0,
                                  file_insertion_enabled=0,
                                  warning_stream=error_stream))

    # Print errors, disregarding unimportant ones
    error_msg = error_stream.getvalue()
    errors = error_msg.split(token)
    success = True
    output = ""

    for error in errors:
        lines = error.splitlines()
        if not lines:
            continue

        m = re.match(r'.*Unknown (?:interpreted text role|directive type) "(.*)".*$', lines[0])
        if m:
            if m.group(1) in ok_unknown_items:
                continue

        output += name + lines[0] + "::\n    " + "\n    ".join(lines[1:]).rstrip() + "\n"
        success = False

    if not success:
        output += "    " + "-"*72 + "\n"
        for lineno, line in enumerate(text.splitlines()):
            output += "    %-4d    %s\n" % (lineno+1, line)
        output += "    " + "-"*72 + "\n\n"

    if dots:
        output_dot('.' if success else 'F')
    return success, output


def output_dot(msg='.', stream=sys.stderr):
    stream.write(msg)
    stream.flush()


def check_rest(module, names, dots=True):
    """
    Check reStructuredText formatting of docstrings

    Returns: [(name, success_flag, output), ...]
    """

    skip_types = (dict, str, unicode, float, int)

    results = []

    if module.__name__[6:] not in OTHER_MODULE_DOCS:
        results += [(module.__name__,) +
                    validate_rst_syntax(inspect.getdoc(module),
                                        module.__name__, dots=dots)]

    for name in names:
        full_name = module.__name__ + '.' + name
        obj = getattr(module, name, None)

        if obj is None:
            results.append((full_name, False, "%s has no docstring" % (full_name,)))
            continue
        elif isinstance(obj, skip_types):
            continue

        if inspect.ismodule(obj):
            text = inspect.getdoc(obj)
        else:
            try:
                text = str(get_doc_object(obj))
            except:
                import traceback
                results.append((full_name, False,
                                "Error in docstring format!\n" +
                                traceback.format_exc()))
                continue

        try:
            src_file = short_path(inspect.getsourcefile(obj))
        except TypeError:
            src_file = None

        if src_file:
            file_full_name = src_file + ':' + full_name
        else:
            file_full_name = full_name

        results.append((full_name,) + validate_rst_syntax(text, file_full_name, dots=dots))

    return results


def check_doctests(module, verbose, dots=True, doctest_warnings=False):
    """Check code in docstrings of the module's public symbols.

    Returns: list of [(item_name, success_flag, output), ...]
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

    class DTRunner(doctest.DocTestRunner):
        DIVIDER = "\n"

        def __init__(self, item_name, checker=None, verbose=None, optionflags=0):
            self._item_name = item_name
            doctest.DocTestRunner.__init__(self, checker=checker, verbose=verbose,
                                           optionflags=optionflags)

        def _report_item_name(self, out, new_line=False):
            if self._item_name is not None:
                if new_line:
                    out("\n")
                self._item_name = None

        def report_start(self, out, test, example):
            self._checker._source = example.source
            return doctest.DocTestRunner.report_start(self, out, test, example)

        def report_success(self, out, test, example, got):
            if self._verbose:
                self._report_item_name(out, new_line=True)
            return doctest.DocTestRunner.report_success(self, out, test, example, got)

        def report_unexpected_exception(self, out, test, example, exc_info):
            self._report_item_name(out)
            return doctest.DocTestRunner.report_unexpected_exception(
                self, out, test, example, exc_info)

        def report_failure(self, out, test, example, got):
            self._report_item_name(out)
            return doctest.DocTestRunner.report_failure(self, out, test,
                                                        example, got)

    class Checker(doctest.OutputChecker):
        obj_pattern = re.compile('at 0x[0-9a-fA-F]+>')
        vanilla = doctest.OutputChecker()
        rndm_markers = {'# random', '# Random', '#random', '#Random', "# may vary"}
        stopwords = {'plt.', '.hist', '.show', '.ylim', '.subplot(',
                     'set_title', 'imshow', 'plt.show', 'ax.axis', 'plt.plot(',
                     '.bar(', '.title', '.ylabel', '.xlabel', 'set_ylim', 'set_xlim'}

        def __init__(self, parse_namedtuples=True, atol=1e-8, rtol=1e-2):
            self.parse_namedtuples = parse_namedtuples
            self.atol, self.rtol = atol, rtol

        def check_output(self, want, got, optionflags):
            # cut it short if they are equal
            if want == got:
                return True

            # skip stopwords in source
            if any(word in self._source for word in self.stopwords):
                return True

            # skip random stuff
            if any(word in want for word in self.rndm_markers):
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
    results = []

    all_success = True
    for name in get_all_dict(module)[0]:
        full_name = module.__name__ + '.' + name

        if full_name in DOCTEST_SKIPLIST:
            continue

        try:
            obj = getattr(module, name)
        except AttributeError:
            import traceback
            results.append((full_name, False,
                            "Missing item!\n" +
                            traceback.format_exc()))
            continue

        finder = doctest.DocTestFinder()
        try:
            tests = finder.find(obj, name, globs=dict(ns))
        except:
            import traceback
            results.append((full_name, False,
                            "Failed to get doctests!\n" +
                            traceback.format_exc()))
            continue

        flags = NORMALIZE_WHITESPACE | ELLIPSIS | IGNORE_EXCEPTION_DETAIL
        runner = DTRunner(full_name, checker=Checker(), optionflags=flags,
                          verbose=verbose)

        output = []
        success = True
        def out(msg):
            output.append(msg)

        class MyStderr(object):
            """Redirect stderr to the current stdout"""
            def write(self, msg):
                if doctest_warnings:
                    sys.stdout.write(msg)
                else:
                    out(msg)

        # Run tests, trying to restore global state afterward
        old_printoptions = np.get_printoptions()
        old_errstate = np.seterr()
        old_stderr = sys.stderr
        cwd = os.getcwd()
        tmpdir = tempfile.mkdtemp()
        sys.stderr = MyStderr()
        try:
            os.chdir(tmpdir)

            # try to ensure random seed is NOT reproducible
            np.random.seed(None)

            for t in tests:
                t.filename = short_path(t.filename, cwd)
                fails, successes = runner.run(t, out=out)
                if fails > 0:
                    success = False
                    all_success = False

            if have_matplotlib:
                plt.close('all')
        finally:
            sys.stderr = old_stderr
            os.chdir(cwd)
            shutil.rmtree(tmpdir)
            np.set_printoptions(**old_printoptions)
            np.seterr(**old_errstate)

        if dots:
            output_dot('.' if success else 'F')

        results.append((full_name, success, "".join(output)))

    return results


def main(argv):
    parser = ArgumentParser(usage=__doc__.lstrip())
    parser.add_argument("module_names", metavar="SUBMODULES", default=list(PUBLIC_SUBMODULES),
                        nargs='*', help="Submodules to check (default: all public)")
    parser.add_argument("--doctests", action="store_true", help="Run also doctests")
    parser.add_argument("-v", "--verbose", action="count", default=0)
    parser.add_argument("--doctest-warnings", action="store_true",
                        help="Enforce warning checking for doctests")
    args = parser.parse_args(argv)

    modules = []
    names_dict = {}

    os.environ['SCIPY_PIL_IMAGE_VIEWER'] = 'true'

    module_names = list(args.module_names)
    for name in list(module_names):
        if name in OTHER_MODULE_DOCS:
            name = OTHER_MODULE_DOCS[name]
            if name not in module_names:
                module_names.append(name)

    for submodule_name in module_names:
        module_name = BASE_MODULE + '.' + submodule_name
        __import__(module_name)
        module = sys.modules[module_name]

        if submodule_name not in OTHER_MODULE_DOCS:
            find_names(module, names_dict)

        if submodule_name in args.module_names:
            modules.append(module)

    dots = True
    success = True
    results = []

    print("Running checks for %d modules:" % (len(modules),))

    for module in modules:
        if dots:
            if module is not modules[0]:
                sys.stderr.write(' ')
            sys.stderr.write(module.__name__ + ' ')
            sys.stderr.flush()

        all_dict, deprecated, others = get_all_dict(module)
        names = names_dict.get(module.__name__, set())

        mod_results = []
        mod_results += check_items(all_dict, names, deprecated, others, module.__name__)
        mod_results += check_rest(module, set(names).difference(deprecated),
                                  dots=dots)
        if args.doctests:
            mod_results += check_doctests(module, (args.verbose >= 2), dots=dots,
                                          doctest_warnings=args.doctest_warnings)

        for v in mod_results:
            assert isinstance(v, tuple), v

        results.append((module, mod_results))

    if dots:
        sys.stderr.write("\n")
        sys.stderr.flush()

    # Report results
    all_success = True

    for module, mod_results in results:
        success = all(x[1] for x in mod_results)
        all_success = all_success and success

        if success and args.verbose == 0:
            continue

        print("")
        print("=" * len(module.__name__))
        print(module.__name__)
        print("=" * len(module.__name__))
        print("")

        for name, success, output in mod_results:
            if name is None:
                if not success or args.verbose >= 1:
                    print(output.strip())
                    print("")
            elif not success or (args.verbose >= 2 and output.strip()):
                print(name)
                print("-"*len(name))
                print("")
                print(output.strip())
                print("")

    if all_success:
        print("\nOK: refguide and doctests checks passed!")
        sys.exit(0)
    else:
        print("\nERROR: refguide or doctests have errors")
        sys.exit(1)


if __name__ == '__main__':
    main(argv=sys.argv[1:])
