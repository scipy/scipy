'''Create shared libraries for use within scipy.'''

# Define some things for the module
MODULE_NAME = 'pyHiGHS'
VERSION = '0.4.0'

# Dependencies
CYTHON_VERSION = '0.29.16'
NUMPY_VERSION = '1.18.2'

from distutils.core import setup
from distutils.extension import Extension
from distutils.ccompiler import new_compiler
from distutils.util import get_platform
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext as _build_ext
from Cython.Build import cythonize

import os
import sys
from datetime import datetime
import pathlib
import sysconfig

# see https://stackoverflow.com/questions/14320220/testing-python-c-libraries-get-build-path
def get_distutils_lib_path():
    PLAT_SPEC = "%s-%d.%d" % (get_platform(), *sys.version_info[:2])
    return os.path.join("build", "lib.%s" % PLAT_SPEC)

class build_ext(_build_ext):
    '''Subclass build_ext to bootstrap numpy.'''

    def get_export_symbols(self, ext):
        return ext.export_symbols

    def finalize_options(self):

        # Make sure directory of shared libraries is in rpath
        # for inplace builds
        self.user = True

        _build_ext.finalize_options(self)

        # Modify the rpath to insert the directory where the SOs are installed;
        # it will be a relative path, so module has it as it's working
        # directory
        if sys.platform != 'win32':
            self.rpath.append('./')

        # Prevent numpy from thinking it's still in its setup process
        import numpy as np
        self.include_dirs.append(np.get_include())

# Create HConfig.h: this is usually created by cmake,
# but we just need an empty file and we'll do the
# pound defines here in setup.py
HConfig_h = pathlib.Path('src/HConfig.h')
if not HConfig_h.exists():
    HConfig_h.touch()

def get_sources(CMakeLists, start_token, end_token):
    # Read in sources from CMakeLists.txt
    with open(CMakeLists, 'r') as f:
        s = f.read()

        # Find block where sources are listed
        start_idx = s.find(start_token) + len(start_token)
        end_idx = s[start_idx:].find(end_token) + len(s[:start_idx])
        sources = s[start_idx:end_idx].split('\n')
        sources = [s.strip() for s in sources if s[0] != '#']

    # Make relative to setup.py
    sources = [str(pathlib.Path('src/' + s)) for s in sources]
    return sources

# Get sources for each extension
sources = get_sources('src/CMakeLists.txt', 'set(sources\n', ')')
basiclu_sources = get_sources('src/CMakeLists.txt', 'set(basiclu_sources\n', ')')
ipx_sources = get_sources('src/CMakeLists.txt', 'set(ipx_sources\n', ')')

# Grab some more info about HiGHS from root CMakeLists
def get_version(CMakeLists, start_token, end_token=')'):
    with open('CMakeLists.txt', 'r') as f:
        s = f.read()
        start_idx = s.find(start_token) + len(start_token) + 1
        end_idx = s[start_idx:].find(end_token) + len(s[:start_idx])
    return s[start_idx:end_idx].strip()
HIGHS_VERSION_MAJOR = get_version('CMakeLists.txt', 'HIGHS_VERSION_MAJOR')
HIGHS_VERSION_MINOR = get_version('CMakeLists.txt', 'HIGHS_VERSION_MINOR')
HIGHS_VERSION_PATCH = get_version('CMakeLists.txt', 'HIGHS_VERSION_PATCH')

# Get path to shared libraries (for local build only)
CYTHON_DIR = pathlib.Path(__file__).parent / MODULE_NAME
HIGHS_DIR = str(CYTHON_DIR.resolve().parent)

LIBRARY_DIRS = []
LIBRARY_DIRS.append(str(CYTHON_DIR.parent / get_distutils_lib_path() / MODULE_NAME))

# Windows can't have runtime_library_dirs
RUNTIME_LIBRARY_DIRS = []
if sys.platform != 'win32':
    RUNTIME_LIBRARY_DIRS = LIBRARY_DIRS

# Read in current GITHASH
#with open('GITHASH', 'r') as f:
#    GITHASH = f.read().strip()
GITHASH = "n/a"

# Here are the pound defines that HConfig.h would usually provide;
# We provide an empty HConfig.h file and do the defs and undefs
# here:
TODAY_DATE = datetime.today().strftime('%Y-%m-%d')
DEFINE_MACROS = [
    #('OPENMP', None),
    ('CMAKE_BUILD_TYPE', '"Release"'),
    ('HiGHSRELEASE', None),
    ('IPX_ON', None),
    ('HIGHS_GITHASH', '"%s"' % GITHASH),
    ('HIGHS_COMPILATION_DATE', '"' + TODAY_DATE +'"'),
    ('HIGHS_VERSION_MAJOR', HIGHS_VERSION_MAJOR),
    ('HIGHS_VERSION_MINOR', HIGHS_VERSION_MINOR),
    ('HIGHS_VERSION_PATCH', HIGHS_VERSION_PATCH),
    ('HIGHS_DIR', '"' + HIGHS_DIR + '"'),
]
UNDEF_MACROS = [
    'OPENMP', # unconditionally disable openmp
    'EXT_PRESOLVE',
    'SCIP_DEV',
    'HiGHSDEV',
    'OSI_FOUND',
]

# Naming conventions of shared libraries differ platform to platform (used for unix build):
SO_PREFIX = str(pathlib.Path(new_compiler().library_filename('', lib_type='shared')).with_suffix(''))
SO_SUFFIX = str(pathlib.Path(sysconfig.get_config_var('EXT_SUFFIX')).with_suffix(''))
if SO_SUFFIX is None:
    # https://bugs.python.org/issue19555
    SO_SUFFIX = str(pathlib.Path(sysconfig.get_config_var('SO')).with_suffix(''))

# No prefix or suffix for Windows, for some reason it's ".dll" and ".pyd"
if sys.platform == 'win32':
    SO_PREFIX = ''
    SO_SUFFIX = ''

# We use some modern C++, as you should. HiGHS uses C++11, no penalty for going to C++14
EXTRA_COMPILE_ARGS = []
if sys.platform == 'win32':
    EXTRA_COMPILE_ARGS.append('/std:c++14')
else:
    EXTRA_COMPILE_ARGS.append('-std=c++14')

# Shared include dirs
HIGHS_INCLUDE_DIRS = [
    str(pathlib.Path('src/ipm/basiclu/include/')),
    str(pathlib.Path('external/')),
    str(pathlib.Path('src/')),
    str(pathlib.Path('src/ipm/ipx/include/')),
    str(pathlib.Path('src/lp_data/')),
    str(pathlib.Path('src/io/')),
    str(pathlib.Path('src/mip/')),
    str(pathlib.Path('src/interfaces/')),
]

# For Linux/Mac we can split up into shared libraries that link to each other.
# For Mac we need to do this so C sources can be compiled without -std=c++14 option
if sys.platform != 'win32':
    # Create shared library for each component
    extensions = [

        # BASICLU
        Extension(
            MODULE_NAME + '.' + SO_PREFIX + 'basiclu',
            basiclu_sources,
            include_dirs=[
                str(pathlib.Path('src/')),
                str(pathlib.Path('src/ipm/basiclu/include/')),
            ],
            language="c",
            define_macros=DEFINE_MACROS,
            undef_macros=UNDEF_MACROS,
        ),

        # IPX
        Extension(
            MODULE_NAME + '.' + SO_PREFIX + 'ipx',
            ipx_sources,
            include_dirs=[
                str(pathlib.Path('src/')),
                str(pathlib.Path('src/ipm/ipx/include/')),
                str(pathlib.Path('src/ipm/basiclu/include/')),
            ],
            language="c++",
            libraries=['basiclu' + SO_SUFFIX],
            library_dirs=LIBRARY_DIRS,
            runtime_library_dirs=LIBRARY_DIRS,
            define_macros=DEFINE_MACROS,
            undef_macros=UNDEF_MACROS,
            extra_compile_args=EXTRA_COMPILE_ARGS,
        ),

        # HiGHS
        Extension(
            MODULE_NAME + '.libhighs',
            sources,
            include_dirs=[
                str(pathlib.Path(MODULE_NAME + '/src/')),
                str(pathlib.Path('src/')),
                str(pathlib.Path('src/ipm/ipx/include/')),
                str(pathlib.Path('src/lp_data/')),
            ],
            language="c++",
            library_dirs=LIBRARY_DIRS,
            runtime_library_dirs=LIBRARY_DIRS,
            libraries=['ipx' + SO_SUFFIX],
            define_macros=DEFINE_MACROS,
            undef_macros=UNDEF_MACROS,
            extra_compile_args=EXTRA_COMPILE_ARGS
        ),

        # Wrapper over HiGHS C++ API
        Extension(
            MODULE_NAME + '.highs_wrapper',
            [
                str(pathlib.Path(MODULE_NAME + '/src/highs_wrapper.pyx'))
            ],
            include_dirs=[
                str(pathlib.Path(MODULE_NAME+ '/src/')),
            ] + HIGHS_INCLUDE_DIRS,
            language="c++",
            define_macros=DEFINE_MACROS,
            undef_macros=UNDEF_MACROS,
            extra_compile_args=EXTRA_COMPILE_ARGS,
            libraries=['highs' + SO_SUFFIX],
            library_dirs=LIBRARY_DIRS,
            runtime_library_dirs=LIBRARY_DIRS,
        ),
    ]
else:
    # For Windows, it's a pain in the butt to export all symbols to DLLs;
    # we'll take the easy way out and compile everything from source to
    # create a single extension
    extensions = [
        Extension(
            MODULE_NAME + '.highs_wrapper',
            [
                str(pathlib.Path(MODULE_NAME + '/src/highs_wrapper.pyx'))
            ] + basiclu_sources + ipx_sources + sources,
            include_dirs=[
                str(pathlib.Path(MODULE_NAME+ '/src/')),
            ] + HIGHS_INCLUDE_DIRS,
            language="c++",
            define_macros=DEFINE_MACROS,
            undef_macros=UNDEF_MACROS,
            extra_compile_args=EXTRA_COMPILE_ARGS,
        ),
    ]

setup(
    name='scikit-highs',
    version=VERSION,
    author='Nicholas McKibben',
    author_email='nicholas.bgp@gmail.com',
    packages=find_packages(),
    scripts=[],
    url='https://github.com/mckib2/HiGHS',
    license='MIT',
    description='Cython interface to HiGHS.',
    long_description=open('PYREADME.rst', encoding='utf-8').read(),
    install_requires=[
        "numpy>=" + NUMPY_VERSION,
        "Cython>=" + CYTHON_VERSION,
    ],
    cmdclass={
        'build_ext': build_ext,
    },
    setup_requires=['numpy', 'Cython'],
    python_requires='>=3',

    ext_modules=cythonize(extensions),
)
