def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('sandbox',parent_package,top_path)

    # All subpackages should be commented out in the version
    # committed to the repository. This prevents build problems
    # for people who are not actively working with these 
    # potentially unstable packages.

    # An example package:
    #config.add_subpackage('exmplpackage')

    # Monte Carlo package
    #config.add_subpackage('montecarlo')

    # PySparse fork with NumPy compatibility
    #config.add_subpackage('pysparse')
    
    # Robert Kern's corner:
    #config.add_subpackage('rkern')

    # ODRPACK
    #config.add_subpackage('odr')

    # Delaunay triangulation and Natural Neighbor interpolation
    #config.add_subpackage('delaunay')

    # Gist-based plotting library for X11
    #config.add_subpackage('xplt')

    #config.add_subpackage('nd_image')

    # elementwise numerical expressions
    #config.add_subpackage('numexpr')

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
