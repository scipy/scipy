import inspect
import pathlib

import scipy.stats
import numpy as np

def pre_build_hook(build_ext, ext):
    from scipy._build_utils.compiler_helper import get_cxx_std_flag
    std_flag = get_cxx_std_flag(build_ext._cxx_compiler)
    if std_flag is not None:
        ext.extra_compile_args.append(std_flag)

def configuration(parent_package='', top_path=None):
    from scipy._build_utils.boostinator import get_include_dir
    from scipy.stats.boost._generate_pyx import klass_mapper
    from numpy.distutils.misc_util import Configuration
    config = Configuration('boost', parent_package, top_path)

    DEFINES = [
        ('BOOST_MATH_DOMAIN_ERROR_POLICY', 'ignore_error'),  # return nan instead of throwing
        ('BOOST_MATH_PROMOTE_DOUBLE_POLICY', 'false'),
    ]
    INCLUDES = [
        'include/',
        'src/',
        str(get_include_dir()),
        np.get_include(),
    ]

    # generate the PXD and PYX wrappers
    src_dir = pathlib.Path(__file__).parent / 'src'
    for b, s in klass_mapper.items():
        if s is None:
            print(f'{b} has no scipy equivalent! Skipping!')
            continue
        scipy_name = s.__class__.__name__.split('_gen')[0]
        ext = config.add_extension(
            f'{scipy_name}_ufunc',
            sources=[f'{src_dir}/{scipy_name}_ufunc.cxx'],
            include_dirs=INCLUDES,
            define_macros=DEFINES,
            language='c++',
        )
        # Add c++11/14 support:
        ext._pre_build_hook = pre_build_hook

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
