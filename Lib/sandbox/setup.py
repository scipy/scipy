def configuration(parent_package='',top_path=None):
    from scipy.distutils.misc_util import Configuration
    config = Configuration('sandbox',parent_package,top_path)

    # An example package:
    #config.add_subpackage('exmplpackage')

    # Maximum entropy package
    #config.add_subpackage('maxent')
    
    # Robert Kern's corner:
    #config.add_subpackage('rkern')
    
    # ODRPACK
    #config.add_subpackage('odr')

    config.add_subpackage('xplt')
    
    return config

if __name__ == '__main__':
    from scipy.distutils.core import setup
    setup(**configuration(top_path='').todict())
