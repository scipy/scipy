#!/usr/bin/env python
import os
from numpy.distutils.system_info import system_info, dict_append, \
    NotFoundError

class NetCDFNotFound(NotFoundError):
    """
    NetCDF (http://www.unidata.ucar.edu/software/netcdf/) not found.
    Directories to search for the libraries can be specified in the
    numpy/distutils/site.cfg file (section [netcdf]) or by setting
    the NETCDF variable.
    """

class netcdf_info(system_info):
    section = 'netcdf'
    dir_env_var = 'NETCDF'
    notfounderror = NetCDFNotFound

    def __init__(self, *args, **kw):
        system_info.__init__(self, *args, **kw)
        self.cp.defaults()['libraries'] = 'netcdf'

    def calc_info(self):
        info = self.calc_libraries_info()
        include_dirs = self.get_include_dirs()
        inc_dir = None
        for d in include_dirs:
            if self.combine_paths(d, 'netcdf.h'):
                inc_dir = d
                break
        if inc_dir is not None:
            dict_append(info, include_dirs=[inc_dir])
        self.set_info(**info)

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('netcdf',parent_package,top_path)
    netcdf = netcdf_info().get_info()
    config.add_extension('_netcdf',
                         sources = ['_netcdf.c'],
                         depends = ['_netcdf.h'],
                         extra_info = netcdf)
    config.add_data_files('demomodule.c','Scientific_LICENSE',
                          'README')
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    config = configuration(top_path='').todict() 
    setup(**config)
