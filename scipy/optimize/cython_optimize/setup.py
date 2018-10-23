from __future__ import division, print_function, absolute_import
import os


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs
    config = Configuration('cython_optimize', parent_package, top_path)
    config.add_data_dir('tests')
    config.add_extension('zeros_tuple', sources=['zeros_tuple.c'],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('zeros_struct', sources=['zeros_struct.c'],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('zeros_array', sources=['zeros_array.c'],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('tests.examples.zeros_tuple_examples',
                         sources=[os.path.join('tests', 'examples',
                             'zeros_tuple_examples.c')],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('tests.examples.zeros_struct_examples',
                         sources=[os.path.join('tests', 'examples',
                             'zeros_struct_examples.c')],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('tests.examples.zeros_struct_alt_examples',
                         sources=[os.path.join('tests', 'examples',
                             'zeros_struct_alt_examples.c')],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('tests.examples.zeros_array_examples',
                         sources=[os.path.join('tests', 'examples',
                             'zeros_array_examples.c')],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_subpackage('tests.examples')
    config.add_data_files('*.pxd')
    config.add_data_files(os.path.join('examples', '*.pxd'))
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
