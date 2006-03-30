#!/usr/bin/env python
import os

netcdf_prefix = None

if netcdf_prefix is None:
    try:
        netcdf_prefix = os.environ['NETCDF_PREFIX']
    except KeyError:
        for netcdf_prefix in ['/usr/local', '/usr', '/sw']:
            netcdf_include = os.path.join(netcdf_prefix, 'include')
            netcdf_lib = os.path.join(netcdf_prefix, 'lib')
            if os.path.exists(os.path.join(netcdf_include, 'netcdf.h')):
                break
        else:
            netcdf_prefix = None

if netcdf_prefix is None:
    print "netCDF not found, the netCDF module will not be built!"
    print "If netCDF is installed somewhere on this computer,"
    print "please set NETCDF_PREFIX to the path where"
    print "include/netcdf.h and lib/netcdf.a are located"
    print "and re-run the build procedure."
        
def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('netcdf',parent_package,top_path)
    config.add_extension('_netcdf',
                         sources = ['_netcdf.c'],
                         include_dirs=[netcdf_include],
                         library_dirs=[netcdf_lib],
                         libraries = ['netcdf'])
    config.add_data_files(['netcdf',('demomodule.c','Scientific_LICENSE',
	                             'README')])
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    config = configuration(top_path='').todict() 
    setup(**config)
