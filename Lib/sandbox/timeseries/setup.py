#!/usr/bin/env python
__version__ = '1.0'
__revision__ = "$Revision: 37 $"
__date__     = '$Date: 2006-12-08 14:30:29 -0500 (Fri, 08 Dec 2006) $'

import os
from os.path import join

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs
    nxheader = join(get_numpy_include_dirs()[0],'numpy',)
    confgr = Configuration('timeseries',parent_package,top_path)
    sources = join('src', 'cseries.c')
    confgr.add_extension('cseries',
                         sources=[sources,],
                         include_dirs=[nxheader],
                         )
    confgr.add_data_dir('doc')
    confgr.add_data_dir('examples')
    return confgr

if __name__ == "__main__":
    from numpy.distutils.core import setup
    #setup.update(nmasetup)
    config = configuration(top_path='').todict() 
    setup(**config)