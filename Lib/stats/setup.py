#!/usr/bin/env python

import os

def configuration(parent_package='',top_path=None):
    from scipy.distutils.misc_util import Configuration
    config = Configuration('stats', parent_package, top_path)

    config.add_data_dir('tests')

    config.add_library('statlib',
                       sources=[os.path.join('statlib', '*.f')])
    
    # add statlib module
    config.add_extension('statlib',
        sources=['statlib.pyf'],
        f2py_options=['--no-wrap-functions'],
        libraries=['statlib'],
    )
    
    # add futil module
    config.add_extension('futil',
        sources=['futil.f'],
    )

    # add mvn module
    config.add_extension('mvn',
        sources=['mvn.pyf','mvndst.f'],
    )

    return config

if __name__ == '__main__':    
    from scipy.distutils.core import setup
    setup(**configuration(top_path='').todict())
