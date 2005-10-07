def configuration(parent_package='',top_path=None):
    from scipy.distutils.misc_util import Configuration
    config = Configuration('sandbox',parent_package,top_path)

    # An example package:
    #config.add_subpackage('exmplpackage')

    # Robert Kern's corner:
    #config.add_subpackage('rkern')
    
    # ODRPACK
    #config.add_subpackage('odr')
    
    return config

if __name__ == '__main__':
    from scipy.distutils.core import setup
    setup(**configuration(top_path='').todict())
