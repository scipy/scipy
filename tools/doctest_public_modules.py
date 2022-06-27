import os
import sys
import importlib

from scpdt import testmod, testfile, DTConfig
from scpdt.util import get_all_list


BASE_MODULE = "scipy"

PUBLIC_SUBMODULES = [
    'cluster',
    'cluster.hierarchy',
    'cluster.vq',
    'constants',     ## obj.__name__
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

################### A user ctx mgr to turn warnings to errors ###################
from scpdt.util import warnings_errors

config = DTConfig()
config.user_context_mgr = warnings_errors
############################################################################

#LOGFILE = open('doctest.log', 'a')



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

#        LOGFILE.write(module_name + '\n')
#        LOGFILE.write("="*len(module_name)  + '\n')
#        for entry in history:
#            LOGFILE.write(str(entry) + '\n')

        sys.stderr.write(str(result))
        all_success = all_success and (result.failed == 0)
    return all_success


def doctest_single_file(fname):
    result, history = testfile(fname, config=config, module_relative=False)
    return result.failed == 0


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
    parser.add_argument( "-t", "--tests", default=None,
                        help="Specify a .py file to check")
    args = parser.parse_args()

    # TODO: verbosity

    if args.submodule and args.tests:
        raise ValueError("Specify either a submodule or a single file, not both.")

    if args.tests:
        all_success = doctest_single_file(args.tests)
    else:
        name = args.submodule   # XXX : dance w/ subsubmodules : cluster.vq etc
        submodule_names = [name]  if name else list(PUBLIC_SUBMODULES)
        all_success = doctest_submodules(submodule_names,
                                         verbose=args.verbose,
                                         fail_fast=args.fail_fast)

    # final report
    if all_success:
        sys.stderr.write('\n\n>>>> OK: doctests PASSED\n')
        sys.exit(0)
    else:
        sys.stderr.write('\n\n>>>> ERROR: doctests FAILED\n')
        sys.exit(1)

