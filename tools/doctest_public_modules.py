import os
import glob
import sys
import importlib

from scpdt import testmod, testfile, DTConfig
from scpdt.util import get_all_list


BASE_MODULE = "scipy"

PUBLIC_SUBMODULES = [
    'cluster',
    'cluster.hierarchy',
    'cluster.vq',
    'constants',
    'fft',
    'fftpack',
    'fftpack.convolve',   # has zero doctests?
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
    'datasets',
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

################### A user ctx mgr to turn warnings to errors ###################
from scpdt.util import warnings_errors
from contextlib import contextmanager
import warnings


# FIXME: populate the dict once
@contextmanager
def warnings_errors_and_rng(test):
    """Temporarily turn (almost) all warnings to errors.

    Filter out known warnings which we allow.
    """
    known_warnings = dict()

    # these functions are known to emit "divide by zero" RuntimeWarnings
    divide_by_zero = [
        'scipy.linalg.norm', 'scipy.ndimage.center_of_mass',
    ]
    for name in divide_by_zero:
        known_warnings[name] = dict(category=RuntimeWarning,
                                    message='divide by zero')

    # Deprecated stuff in scipy.signal and elsewhere
    deprecated = [
        'scipy.signal.cwt', 'scipy.signal.morlet', 'scipy.signal.morlet2',
        'scipy.signal.ricker',
        'scipy.integrate.simpson',
        'scipy.interpolate.interp2d',
    ]
    for name in deprecated:
        known_warnings[name] = dict(category=DeprecationWarning)

    from scipy import integrate
    # the funcions are known to emit IntergrationWarnings
    integration_w = ['scipy.special.ellip_normal',
                     'scipy.special.ellip_harm_2',
    ]
    for name in integration_w:
        known_warnings[name] = dict(category=integrate.IntegrationWarning,
                                    message='The occurrence of roundoff')

    # scipy.stats deliberately emits UserWarnings sometimes
    user_w = ['scipy.stats.anderson_ksamp', 'scipy.stats.kurtosistest',
              'scipy.stats.normaltest']
    for name in user_w:
        known_warnings[name] = dict(category=UserWarning)

    # additional one-off warnings to filter
    dct = {
        'scipy.sparse.linalg.norm':
            dict(category=UserWarning, message="Exited at iteration"),
        # tutorials
        'linalg.rst':
            dict(message='the matrix subclass is not',
                 category=PendingDeprecationWarning),
        'stats.rst':
            dict(message='The maximum number of subdivisions',
                 category=integrate.IntegrationWarning),
    }
    known_warnings.update(dct)

    # these legitimately emit warnings in examples
    from scipy.signal._filter_design import BadCoefficients
    legit = set('scipy.signal.normalize')

    # Now, the meat of the matter: filter warnings,
    # also control the random seed for each doctest.

    # XXX: this matches the refguide-check behavior, but is a tad strange:
    # makes sure that the seed the old-fashioned np.random* methods is *NOT*
    # reproducible but the new-style `default_rng()` *IS* repoducible.
    # Should these two be either both repro or both not repro?

    from scipy._lib._util import _fixed_default_rng
    import numpy as np
    with _fixed_default_rng():
        np.random.seed(None)
        with warnings.catch_warnings():
            if test.name in known_warnings:
                warnings.filterwarnings('ignore', **known_warnings[test.name])
                yield
            elif test.name in legit:
                yield
            else:
                warnings.simplefilter('error', Warning)
                yield


config = DTConfig()
config.user_context_mgr = warnings_errors_and_rng
config.skiplist = set([
    'scipy.linalg.LinAlgError',     # comes from numpy
    'scipy.fftpack.fftshift',       # fftpack stuff is also from numpy
    'scipy.fftpack.ifftshift',
    'scipy.fftpack.fftfreq',
    'scipy.special.sinc',           # sinc is from numpy
    'scipy.optimize.show_options',  # does not have much to doctest
    'scipy.signal.normalize',       # manipulates warnings (XXX temp skip)
])
############################################################################

LOGFILE = open('doctest.log', 'a')



def doctest_submodules(module_names, verbose, fail_fast):
    all_success = True
    for submodule_name in module_names:
        prefix = BASE_MODULE + '.'
        if not submodule_name.startswith(prefix):
            module_name = prefix + submodule_name
        else:
            module_name = submodule_name

        module = importlib.import_module(module_name)

        full_name = module.__name__
        line = '='*len(full_name)
        sys.stderr.write(f"\n\n{line}\n")
        sys.stderr.write(full_name)
        sys.stderr.write(f"\n{line}\n")

        result, history = testmod(module, strategy='api',
                                  verbose=verbose,
                                  raise_on_error=fail_fast, config=config) 

        LOGFILE.write(module_name + '\n')
        LOGFILE.write("="*len(module_name)  + '\n')
        for entry in history:
            LOGFILE.write(entry[len(module_name)+1:] + '\n')

        sys.stderr.write(str(result))
        all_success = all_success and (result.failed == 0)
    return all_success


def doctest_single_file(fname, verbose, fail_fast):
    result, history = testfile(fname, config=config, module_relative=False,
                               verbose=verbose, raise_on_error=fail_fast)
    return result.failed == 0


def doctest_tutorial(verbose, fail_fast):
    base_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), '..')
    tut_path = os.path.join(base_dir, 'doc', 'source', 'tutorial', '*.rst')
    tut_path2 = os.path.join(base_dir, 'doc', 'source', 'tutorial', '*', '*.rst')
    sys.stderr.write('\nChecking tutorial files at %s:\n'
                     % os.path.relpath(tut_path, os.getcwd()))

    tutorials = [
        f for f in sorted(glob.glob(tut_path)) + sorted(glob.glob(tut_path2))
    ]

    # set up scipy-specific config
    config.pseudocode = set(['integrate.nquad(func,'])

    io_matfiles = glob.glob(os.path.join(tut_path.replace('.rst', '.mat')))
    config.local_resources = {'io.rst': io_matfiles}

    skip_these = [
        'sampling_srou', 'sampling_pinv',   # various failures
        'ND_regular_grid',    # needs rst parser
        'ND_unstructured',
        'extrapolation_examples'
    ]

    all_success = True
    for filename in tutorials:
        sys.stderr.write('\n' + filename + '\n')
        sys.stderr.write("="*len(filename) + '\n')

        # XXX: temporary skips
        if any(s in filename for s in skip_these):
            sys.stderr.write("SKIPPED\n")
            continue

        result, history = testfile(filename, module_relative=False,
                                    verbose=verbose, raise_on_error=fail_fast,
                                    report=True, config=config)
        all_success = all_success and (result.failed == 0)
    return all_success


def main(args):
    if args.submodule and args.filename:
        raise ValueError("Specify either a submodule or a single file, not both.")

    if args.filename:
        all_success = doctest_single_file(args.filename,
                                          verbose=args.verbose,
                                          fail_fast=args.fail_fast)
    else:
        all_success = True
        name = args.submodule   # XXX : dance w/ subsubmodules : cluster.vq etc
        submodule_names = [name]  if name else list(PUBLIC_SUBMODULES)
        all_success = doctest_submodules(submodule_names,
                                         verbose=args.verbose,
                                         fail_fast=args.fail_fast)

        # if full run: also check the tutorial
        if not args.submodule:
            tut_success = doctest_tutorial(verbose=args.verbose,
                                           fail_fast=args.fail_fast)
            all_success = all_success and tut_success

    LOGFILE.close()

    # final report
    if all_success:
        sys.stderr.write('\n\n>>>> OK: doctests PASSED\n')
        sys.exit(0)
    else:
        sys.stderr.write('\n\n>>>> ERROR: doctests FAILED\n')
        sys.exit(1)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="doctest runner")
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help='print verbose (`-v`) or very verbose (`-vv`) '
                              'output for all tests')
    parser.add_argument('-x', '--fail-fast', action='store_true',
                        help=('stop running tests after first failure'))
    parser.add_argument( "-s", "--submodule", default=None,
                        help="Submodule whose tests to run (cluster,"
                             " constants, ...)")
    parser.add_argument( "-t", "--filename", default=None,
                        help="Specify a .py file to check")
    args = parser.parse_args()

    main(args)

