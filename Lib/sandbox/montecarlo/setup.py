import numpy
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration
from os.path import join

def configuration(parent_package='', top_path=None):

    config = Configuration('montecarlo', parent_package, top_path)

    config.add_extension('_intsampler',
                         include_dirs = [numpy.get_numpy_include()],
                         sources = [join('src',f) for f in
                         ['_intsamplermodule.c', 'sampler5tbl.c']] )

    config.add_data_dir('tests')
    config.add_data_dir('examples')
    config.add_data_dir('doc')

    return config

if __name__ == '__main__':
    setup(**configuration(top_path='').todict())
