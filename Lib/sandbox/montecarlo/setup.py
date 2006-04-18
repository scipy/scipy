import numpy
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration
from os.path import join

def configuration(parent_package='', top_path=None):

    config = Configuration('montecarlo', parent_package, top_path)

    config.add_library('randomkit', 
            # The following path needs to be extracted in a portable way.
            # This points to the randomkit.c file in the numpy source tree.
            sources=['/home/schofield/Install/numpy/numpy/random/mtrand/randomkit.c'])

    config.add_extension('_intsampler',
              include_dirs = [numpy.get_numpy_include(),
                    # The following path needs to be extracted in a portable way.
                    # This points to the default installation location used by
                    # config.add_header().
                   '/usr/include/python2.4/numpy/random/'],
              libraries=['randomkit'],
              sources = [join('src', f) for f in
                        ['_intsamplermodule.c', 'compact5table.c']]
              ) 
    config.add_data_dir('tests')
    config.add_data_dir('examples')
    config.add_data_dir('doc')

    return config

if __name__ == '__main__':
    setup(**configuration(top_path='').todict())
