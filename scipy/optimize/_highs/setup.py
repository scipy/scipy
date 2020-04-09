
import sys
import pathlib
from datetime import datetime

def _get_sources(CMakeLists, start_token, end_token):
    # Read in sources from CMakeLists.txt
    CMakeLists = pathlib.Path(__file__).parent / CMakeLists
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

# Grab some more info about HiGHS from root CMakeLists
def _get_version(CMakeLists, start_token, end_token=')'):
    CMakeLists = pathlib.Path(__file__).parent / CMakeLists
    with open(CMakeLists, 'r') as f:
        s = f.read()
        start_idx = s.find(start_token) + len(start_token) + 1
        end_idx = s[start_idx:].find(end_token) + len(s[:start_idx])
    return s[start_idx:end_idx].strip()

def configuration(parent_package='', top_path=None):

    from numpy.distutils.misc_util import Configuration
    config = Configuration('_highs', parent_package, top_path)

    # HiGHS info
    HIGHS_VERSION_MAJOR = _get_version('CMakeLists.txt', 'HIGHS_VERSION_MAJOR')
    HIGHS_VERSION_MINOR = _get_version('CMakeLists.txt', 'HIGHS_VERSION_MINOR')
    HIGHS_VERSION_PATCH = _get_version('CMakeLists.txt', 'HIGHS_VERSION_PATCH')
    GITHASH = 'n/a'
    HIGHS_DIR = str(pathlib.Path(__file__).parent.resolve())

    # Here are the pound defines that HConfig.h would usually provide;
    # We provide an empty HConfig.h file and do the defs and undefs
    # here:
    TODAY_DATE = datetime.today().strftime('%Y-%m-%d')
    DEFINE_MACROS = [
        #('OPENMP', None),
        ('CMAKE_BUILD_TYPE', '"Release"'),
        ('HiGHSRELEASE', None),
        ('IPX_ON', 'ON'),
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

    # HiGHS written using C++11 standard, we use C++14 for std::make_unique
    EXTRA_COMPILE_ARGS = []
    if sys.platform == 'win32':
        EXTRA_COMPILE_ARGS.append('/std:c++14')
    else:
        EXTRA_COMPILE_ARGS.append('-std=c++14')

    # -O3 takes longer to compile, but passes unit tests -- could use?
    #EXTRA_COMPILE_ARGS.append('-O3')

    # Compile BASICLU as a static library to appease clang:
    # (won't allow -std=c++14 option for C sources)
    basiclu_sources = _get_sources('src/CMakeLists.txt', 'set(basiclu_sources\n', ')')
    config.add_library(
        'basiclu',
        sources=basiclu_sources,
        include_dirs=[
            str(pathlib.Path('src/')),
            str(pathlib.Path('src/ipm/basiclu/include/')),
        ],
        language='c',
        macros=DEFINE_MACROS,
    )

    # Compile the rest of the sources all together, linking the BASICLU static library
    ipx_sources = _get_sources('src/CMakeLists.txt', 'set(ipx_sources\n', ')')
    highs_sources = _get_sources('src/CMakeLists.txt', 'set(sources\n', ')')
    WRAPPER_INCLUDE_DIRS = [
        str(pathlib.Path('src/ipm/basiclu/include/')),
        str(pathlib.Path('external/')),
        str(pathlib.Path('src/')),
        str(pathlib.Path('src/ipm/ipx/include/')),
        str(pathlib.Path('src/lp_data/')),
        str(pathlib.Path('src/io/')),
        str(pathlib.Path('src/mip/')),
        str(pathlib.Path('src/interfaces/')),
        str(pathlib.Path('pyHiGHS/src/')),
    ]
    config.add_extension(
        'highs_wrapper',
        sources=[
            str(pathlib.Path('pyHiGHS/src/highs_wrapper.cxx'))
        ] + ipx_sources + highs_sources,
        include_dirs=[
            str(pathlib.Path('pyHiGHS/src/')),
        ] + WRAPPER_INCLUDE_DIRS,
        language='c++',
        libraries=['basiclu'],
        define_macros=DEFINE_MACROS,
        undef_macros=UNDEF_MACROS,
        extra_compile_args=EXTRA_COMPILE_ARGS,
    )

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
