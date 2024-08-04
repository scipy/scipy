#!/usr/bin/env python3
"""
refguide_check.py [OPTIONS] [-- ARGS]

Check for a Scipy submodule whether the objects in its __all__ dict
correspond to the objects included in the reference guide.

Example of usage::

    $ python3 refguide_check.py optimize

Note that this is a helper script to be able to check if things are missing;
the output of this script does need to be checked manually.  In some cases
objects are left out of the refguide for a good reason (it's an alias of
another function, or deprecated, or ...)

"""
import copy
import inspect
import io
import os
import re
import sys
import warnings
from argparse import ArgumentParser

import docutils.core
from docutils.parsers.rst import directives

from numpydoc.docscrape_sphinx import get_doc_object
from numpydoc.docscrape import NumpyDocString
from scipy.stats._distr_params import distcont, distdiscrete
from scipy import stats


# Enable specific Sphinx directives
from sphinx.directives.other import SeeAlso, Only
directives.register_directive('seealso', SeeAlso)
directives.register_directive('only', Only)

BASE_MODULE = "scipy"

PUBLIC_SUBMODULES = [
    'cluster',
    'cluster.hierarchy',
    'cluster.vq',
    'constants',
    'datasets',
    'differentiate',
    'fft',
    'fftpack',
    'fftpack.convolve',
    'integrate',
    'interpolate',
    'io',
    'io.arff',
    'io.matlab',
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
    'signal.windows',
    'sparse',
    'sparse.csgraph',
    'sparse.linalg',
    'spatial',
    'spatial.distance',
    'spatial.transform',
    'special',
    'stats',
    'stats.mstats',
    'stats.contingency',
    'stats.qmc',
    'stats.sampling'
]

# Docs for these modules are included in the parent module
OTHER_MODULE_DOCS = {
    'fftpack.convolve': 'fftpack',
    'io.wavfile': 'io',
    'io.arff': 'io',
}

# these names are not required to be present in ALL despite being in
# autosummary:: listing
REFGUIDE_ALL_SKIPLIST = [
    r'scipy\.sparse\.csgraph',
    r'scipy\.sparse\.linalg',
    r'scipy\.linalg\.blas\.[sdczi].*',
    r'scipy\.linalg\.lapack\.[sdczi].*',
]

# these names are not required to be in an autosummary:: listing
# despite being in ALL
REFGUIDE_AUTOSUMMARY_SKIPLIST = [
    r'scipy\.special\..*_roots',  # old aliases for scipy.special.*_roots
    r'scipy\.special\.jn',  # alias for jv
    r'scipy\.ndimage\.sum',   # alias for sum_labels
    r'scipy\.linalg\.solve_lyapunov',  # deprecated name
    r'scipy\.stats\.contingency\.chi2_contingency',
    r'scipy\.stats\.contingency\.expected_freq',
    r'scipy\.stats\.contingency\.margins',
    r'scipy\.stats\.reciprocal',  # alias for lognormal
    r'scipy\.stats\.trapz',   # alias for trapezoid
]
# deprecated windows in scipy.signal namespace
for name in ('barthann', 'bartlett', 'blackmanharris', 'blackman', 'bohman',
             'boxcar', 'chebwin', 'cosine', 'exponential', 'flattop',
             'gaussian', 'general_gaussian', 'hamming', 'hann', 'hanning',
             'kaiser', 'nuttall', 'parzen', 'triang', 'tukey'):
    REFGUIDE_AUTOSUMMARY_SKIPLIST.append(r'scipy\.signal\.' + name)


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
        res = re.search(
            r"^\s*\.\. (?:currentmodule|module):: ([a-z0-9A-Z_.]+)\s*$",
            line
        )
        if res:
            module_name = res.group(1)
            continue

        for pattern in patterns:
            res = re.match(pattern, line)
            if res is not None:
                name = res.group(1)
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

    others = (set(dir(module))
              .difference(set(deprecated))
              .difference(set(not_deprecated)))

    return not_deprecated, deprecated, others


def compare(all_dict, others, names, module_name):
    """Return sets of objects only in __all__, refguide, or completely missing."""
    only_all = set()
    for name in all_dict:
        if name not in names:
            for pat in REFGUIDE_AUTOSUMMARY_SKIPLIST:
                if re.match(pat, module_name + '.' + name):
                    break
            else:
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
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("error")
        try:
            f(**{"not a kwarg":None})
        except DeprecationWarning:
            return True
        except Exception:
            pass
        return False


def check_items(all_dict, names, deprecated, others, module_name, dots=True):
    num_all = len(all_dict)
    num_ref = len(names)

    output = ""

    output += "Non-deprecated objects in __all__: %i\n" % num_all
    output += "Objects in refguide: %i\n\n" % num_ref

    only_all, only_ref, missing = compare(all_dict, others, names, module_name)
    dep_in_ref = only_ref.intersection(deprecated)
    only_ref = only_ref.difference(deprecated)

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
            output += (
                f"ERROR: objects in {module_name}.__all__ but not in refguide::\n\n"
            )
            for name in sorted(only_all):
                output += "    " + name + "\n"

            output += "\nThis issue can be fixed by adding these objects to\n"
            output += "the function listing in __init__.py for this module\n"

        if len(only_ref) > 0:
            output += (
                f"ERROR: objects in refguide but not in {module_name}.__all__::\n\n"
            )
            for name in sorted(only_ref):
                output += "    " + name + "\n"

            output += "\nThis issue should likely be fixed by removing these objects\n"
            output += "from the function listing in __init__.py for this module\n"
            output += "or adding them to __all__.\n"

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
        return False, f"ERROR: {name}: no documentation"

    ok_unknown_items = set([
        'mod', 'currentmodule', 'autosummary', 'data', 'legacy',
        'obj', 'versionadded', 'versionchanged', 'module', 'class', 'meth',
        'ref', 'func', 'toctree', 'moduleauthor', 'deprecated',
        'sectionauthor', 'codeauthor', 'eq', 'doi', 'DOI', 'arXiv', 'arxiv',
        'versionremoved',
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

        m = re.match(
            r'.*Unknown (?:interpreted text role|directive type) "(.*)".*$',
            lines[0]
        )
        if m:
            if m.group(1) in ok_unknown_items:
                continue

        m = re.match(
            r'.*Error in "math" directive:.*unknown option: "label"', " ".join(lines),
            re.S
        )
        if m:
            continue

        output += (
            name + lines[0] + "::\n    " + "\n    ".join(lines[1:]).rstrip() + "\n"
        )
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
    skip_types = (dict, str, float, int)

    results = []

    if module.__name__[6:] not in OTHER_MODULE_DOCS:
        results += [(module.__name__,) +
                    validate_rst_syntax(inspect.getdoc(module),
                                        module.__name__, dots=dots)]

    for name in names:
        full_name = module.__name__ + '.' + name
        obj = getattr(module, name, None)

        if obj is None:
            results.append((full_name, False, f"{full_name} has no docstring"))
            continue
        elif isinstance(obj, skip_types):
            continue

        if inspect.ismodule(obj):
            text = inspect.getdoc(obj)
        else:
            try:
                text = str(get_doc_object(obj))
            except Exception:
                import traceback
                results.append((full_name, False,
                                "Error in docstring format!\n" +
                                traceback.format_exc()))
                continue

        m = re.search(".*?([\x00-\x09\x0b-\x1f]).*", text)
        if m:
            msg = ("Docstring contains a non-printable character "
                   f"{m.group(1)!r} in the line\n\n{m.group(0)!r}\n\n"
                   "Maybe forgot r\"\"\"?")
            results.append((full_name, False, msg))
            continue

        try:
            src_file = short_path(inspect.getsourcefile(obj))
        except TypeError:
            src_file = None

        if src_file:
            file_full_name = src_file + ':' + full_name
        else:
            file_full_name = full_name

        results.append(
            (full_name,) + validate_rst_syntax(text, file_full_name, dots=dots)
        )

    return results


def check_dist_keyword_names():
    # Look for collisions between names of distribution shape parameters and
    # keywords of distribution methods. See gh-5982.
    distnames = set(distdata[0] for distdata in distcont + distdiscrete)
    mod_results = []
    for distname in distnames:
        dist = getattr(stats, distname)

        method_members = inspect.getmembers(dist, predicate=inspect.ismethod)
        method_names = [method[0] for method in method_members
                        if not method[0].startswith('_')]
        for methodname in method_names:
            method = getattr(dist, methodname)
            try:
                params = NumpyDocString(method.__doc__)['Parameters']
            except TypeError:
                result = (f'stats.{distname}.{methodname}', False,
                          "Method parameters are not documented properly.")
                mod_results.append(result)
                continue

            if not dist.shapes:  # can't have collision if there are no shapes
                continue
            shape_names = dist.shapes.split(', ')

            param_names1 = set(param.name for param in params)
            param_names2 = set(inspect.signature(method).parameters)
            param_names = param_names1.union(param_names2)

            # # Disabling this check in this PR;
            # # these discrepancies are a separate issue.
            # no_doc_params = {'args', 'kwds', 'kwargs'}  # no need to document
            # undoc_params = param_names2 - param_names1 - no_doc_params
            # if un_doc_params:
            #     result = (f'stats.{distname}.{methodname}', False,
            #               f'Parameter(s) {undoc_params} are not documented.')
            #     mod_results.append(result)
            #     continue

            intersection = param_names.intersection(shape_names)

            if intersection:
                message = ("Distribution/method keyword collision: "
                           f"{intersection} ")
                result = (f'stats.{distname}.{methodname}', False, message)
            else:
                result = (f'stats.{distname}.{methodname}', True, '')
            mod_results.append(result)

    return mod_results


def main(argv):
    parser = ArgumentParser(usage=__doc__.lstrip())
    parser.add_argument("module_names", metavar="SUBMODULES", default=[],
                        nargs='*', help="Submodules to check (default: all public)")
    parser.add_argument("-v", "--verbose", action="count", default=0)
    args = parser.parse_args(argv)

    modules = []
    names_dict = {}

    if not args.module_names:
        args.module_names = list(PUBLIC_SUBMODULES)

    os.environ['SCIPY_PIL_IMAGE_VIEWER'] = 'true'

    module_names = list(args.module_names)
    for name in list(module_names):
        if name in OTHER_MODULE_DOCS:
            name = OTHER_MODULE_DOCS[name]
            if name not in module_names:
                module_names.append(name)

    for submodule_name in module_names:
        prefix = BASE_MODULE + '.'
        if not submodule_name.startswith(prefix):
            module_name = prefix + submodule_name
        else:
            module_name = submodule_name

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
        if module.__name__ == 'scipy.stats':
            mod_results += check_dist_keyword_names()

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
        print("\nOK: refguide checks passed!")
        sys.exit(0)
    else:
        print("\nERROR: refguide have errors")
        sys.exit(1)


if __name__ == '__main__':
    main(argv=sys.argv[1:])
