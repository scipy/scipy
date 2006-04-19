import numpy
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration
from os.path import join, dirname

def configuration(parent_package='', top_path=None):

    config = Configuration('montecarlo', parent_package, top_path)

    # This code requires 'randomkit' to have been built using 'add_extension' in
    # numpy/random/setup.py.

    random_lib_dir = dirname(numpy.random.__file__)

    config.add_extension('_intsampler',
              include_dirs = [numpy.get_numpy_include(), random_lib_dir],
              libraries=['randomkit'],
              library_dirs=[random_lib_dir],
              runtime_library_dirs=[random_lib_dir],
              sources = [join('src', f) for f in
                        ['_intsamplermodule.c', 'compact5table.c']])
    
    config.add_data_dir('tests')
    config.add_data_dir('examples')
    config.add_data_dir('doc')

    return config

if __name__ == '__main__':
    setup(**configuration(top_path='').todict())
