#!/usr/bin/env python

import os
from scipy_distutils.misc_util import get_path, default_config_dict, dot_join
from scipy_distutils.core import Extension
from scipy_distutils.x11_info import get_x11_info

def configuration(parent_package=''):
    """
       gist only works with an X-windows server
       This will install *.gs and *.gp files to
       '%spython%s/site-packages/scipy/xplt' % (sys.prefix,sys.version[:3])
    """
    x11_info = get_x11_info()
    if not x11_info:
        return

    config = default_config_dict('xplt',parent_package)
    local_path = get_path(__name__)
    
    sources = ['gistCmodule.c']
    sources = [os.path.join(local_path,x) for x in sources]
                                               
    ext = Extension(dot_join(parent_package,'xplt.gistC'),
                    sources,
                    include_dirs = x11_info.get('include_dirs',[]),
                    library_dirs = x11_info.get('library_dirs',[]),
                    libraries = x11_info.get('libraries',[]) + ['m'])
    config['ext_modules'].append(ext)
        
    from glob import glob
    gist = glob(os.path.join(local_path,'gist','*.c'))
    # libraries are C static libraries
    config['libraries'].append(('gist',{'sources':gist,
                                        'macros':[('STDC_HEADERS',1)]}))
                                   
    file_ext = ['*.gs','*.gp', '*.ps', '*.help']
    xplt_files = [glob(os.path.join(local_path,x)) for x in file_ext]
    xplt_files = reduce(lambda x,y:x+y,xplt_files,[])
    xplt_path = os.path.join(local_path,'xplt')
    config['data_files'].extend( [(xplt_path,xplt_files)])
    
    return config

if __name__ == '__main__':    
    from scipy_distutils.core import setup
    setup(**configuration())
