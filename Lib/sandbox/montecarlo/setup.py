import numpy
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration
from os.path import join, dirname

def configuration(parent_package='', top_path=None):

    config = Configuration('montecarlo', parent_package, top_path)

    # This code requires 'randomkit.c' and 'randomkit.h' to have been copied
    # to (or symlinked to) montecarlo/src/.

    config.add_extension('_intsampler',
              sources = [join('src', f) for f in
                        ['_intsamplermodule.c', 'compact5table.c', 'randomkit.c']])
    
    config.add_data_dir('tests')
    config.add_data_dir('examples')
    config.add_data_dir('doc')

    return config

if __name__ == '__main__':
    setup(**configuration(top_path='').todict())
