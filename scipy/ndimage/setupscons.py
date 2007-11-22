from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration
from numpy import get_include

def configuration(parent_package='', top_path=None):

    config = Configuration('ndimage', parent_package, top_path, 
                           setup_name = 'setupscons.py')

    config.add_sconscript('SConstruct')
    config.add_subpackage('segment')
    config.add_data_dir('tests')

    return config

if __name__ == '__main__':
    setup(**configuration(top_path='').todict())
