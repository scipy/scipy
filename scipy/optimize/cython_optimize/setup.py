from __future__ import division, print_function, absolute_import


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs
    config = Configuration('cython_optimize', parent_package, top_path)
    config.add_data_dir('tests')
    config.add_extension('zeros', sources=['zeros.c'],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('zeros_struct', sources=['zeros_struct.c'],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('tests.zeros_examples',
                         sources=['tests/zeros_examples.c'],
                         include_dirs=[get_numpy_include_dirs()])
    config.add_extension('tests.zeros_struct_examples',
                         sources=['tests/zeros_struct_examples.c'],
                         include_dirs=[get_numpy_include_dirs()])
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
