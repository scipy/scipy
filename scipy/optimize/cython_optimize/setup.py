from __future__ import division, print_function, absolute_import
import os


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs
    config = Configuration('cython_optimize', parent_package, top_path)
    config.add_data_dir('tests')
    newton_dir = 'Newton'
    newton_src = [os.path.join(newton_dir, '*.c')]
    newton_hdr = [os.path.join(newton_dir, 'newton.h')]
    config.add_library('newton',
                       sources=newton_src,
                       headers=newton_hdr)
    config.add_extension('zeros', sources=['zeros.c'],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('zeros_struct', sources=['zeros_struct.c'],
                         libraries=['newton'], depends=(newton_src + newton_hdr),
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('zeros_array', sources=['zeros_array.c'],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('examples.zeros_examples',
                         sources=['examples/zeros_examples.c'],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('examples.zeros_struct_examples',
                         sources=['examples/zeros_struct_examples.c'],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('examples.zeros_struct_alt_examples',
                         sources=['examples/zeros_struct_alt_examples.c'],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('examples.zeros_array_examples',
                         sources=['examples/zeros_array_examples.c'],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_subpackage('examples')
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
