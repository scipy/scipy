#!/usr/bin/env python
# 05.12.2005, c
# last change: 27.03.2006
def configuration(parent_package='',top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    import os.path as op

    config = Configuration( 'umfpack', parent_package, top_path )
    config.add_data_dir('tests')

    umf_info = get_info( 'umfpack', notfound_action = 1 )
    if umf_info:
        print 'Umfpack present, ok.'
    else:
        return None

    scipyInclude = numpy.get_numpy_include()
    umfpackInclude = umf_info['include_dirs'][0]

    config.add_extension( '__umfpack',
                          sources = ['umfpack.i'],
                          swig_opts = ['-I' + umfpackInclude],
                          include_dirs = [umfpackInclude, scipyInclude],
                          libraries = ['cblas'],
                          extra_objects = umf_info['extra_objects'] )
    # config.add_scripts( 'test_umfpack.py' )

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
