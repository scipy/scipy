import os
import sys
import importlib

from scpdt import testmod, DTConfig
from scpdt._util import get_all_list


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
from scpdt._util import warnings_errors

config = DTConfig()
config.user_context_mgr = warnings_errors
############################################################################

module_names = PUBLIC_SUBMODULES

LOGFILE = open('doctest.log', 'a')

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

    result, history = testmod(module, strategy='api', config=config) 

    LOGFILE.write(module_name + '\n')
    LOGFILE.write("="*len(module_name)  + '\n')
    for entry in history:
        LOGFILE.write(str(entry) + '\n')

    sys.stderr.write(str(result))

    all_success = all_success and (result.failed == 0)

# final report
if all_success:
    sys.stderr.write('\n\n>>>> OK: doctests PASSED\n')
    sys.exit(0)
else:
    sys.stderr.write('\n\n>>>> ERROR: doctests FAILED\n')
    sys.exit(-1)
