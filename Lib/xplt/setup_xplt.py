#!/usr/bin/env python

import os
from scipy_distutils.misc_util import get_path, default_config_dict, dot_join
from scipy_distutils.misc_util import dict_append
from scipy_distutils.core import Extension
from scipy_distutils.system_info import x11_info


def configuration(parent_package=''):
    """
       gist only works with an X-windows server
       This will install *.gs and *.gp files to
       '%spython%s/site-packages/scipy/xplt' % (sys.prefix,sys.version[:3])
    """
    x11 = x11_info().get_info()
    if not x11:
        return {}

    config = default_config_dict('xplt',parent_package)
    local_path = get_path(__name__)
    
    sources = ['gistCmodule.c']
    sources = [os.path.join(local_path,x) for x in sources]

    ext_arg = {'name':dot_join(parent_package,'xplt.gistC'),
               'sources':sources}
    dict_append(ext_arg,**x11)
    dict_append(ext_arg,libraries=['m'])
    ext = Extension (**ext_arg)
    config['ext_modules'].append(ext)
        
    from glob import glob
    gist = glob(os.path.join(local_path,'gist','*.c'))
    # libraries are C static libraries
    config['libraries'].append(('gist',{'sources':gist,
                                        'macros':[('STDC_HEADERS',1)]}))
                                   
    file_ext = ['*.gs','*.gp', '*.ps', '*.help']
    xplt_files = [glob(os.path.join(local_path,x)) for x in file_ext]
    xplt_files = reduce(lambda x,y:x+y,xplt_files,[])
    xplt_path = os.path.join(parent_package,'xplt')
    config['data_files'].extend( [(xplt_path,xplt_files)])
    
    return config

if __name__ == '__main__':    
    from scipy_distutils.core import setup
    setup(**configuration())
