#!/usr/bin/env python
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('buildgrid',parent_package,top_path,
                           version='0.5',
                           description='Build grid from scattered '\
                           'data in file/lists',
                           author='Eugene Druker',
                           author_email='eugene.druker@gmail.com')
                          
    netcdf = netcdf_info().get_info()
    config.add_extension('build_grid',
                         sources = ['build_grid.c'])
    config.add_data_files('README')
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    config = configuration(top_path='').todict() 
    setup(**config)
