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

    # Delaunay triangulation and Natural Neighbor interpolation
    #config.add_subpackage('delaunay')

    # Gist-based plotting library for X11
    #config.add_subpackage('xplt')
    
    #config.add_subpackage('nd_image')
    
    return config

if __name__ == '__main__':
    from scipy.distutils.core import setup
    setup(**configuration(top_path='').todict())
