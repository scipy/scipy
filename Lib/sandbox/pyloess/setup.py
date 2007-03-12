#!/usr/bin/env python
__version__ = '1.0'
__revision__ = "$Revision: 2811 $"
__date__     = '$Date: 2007-03-02 06:30:02 -0500 (Fri, 02 Mar 2007) $'

import os
from os.path import join

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    confgr = Configuration('pyloess',parent_package,top_path)
    confgr.add_extension('_lowess',
                         sources=[join('src', 'f_lowess.pyf'),
                                  join('src', 'lowess.f'),]
                         )
    confgr.add_extension('_stl',
                         sources=[join('src', 'f_stl.pyf'),
                                  join('src', 'stl.f')],
                         )
    confgr.add_data_dir('tests')
    return confgr

if __name__ == "__main__":
    from numpy.distutils.core import setup
    config = configuration(top_path='').todict() 
    setup(**config)