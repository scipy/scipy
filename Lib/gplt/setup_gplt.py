import os, sys
from scipy_distutils.misc_util import get_path, default_config_dict

def configuration(parent_package=''):
    parent_path = parent_package
    if parent_package:
        parent_package += '.'
    local_path = get_path(__name__)

    config = default_config_dict()
    config['packages'].append(parent_package+'gplt')
    #config['packages'].append(parent_package+'gplt.tests') 

    if sys.platform == 'win32':
        data = ['wgnuplot.exe', 'gnuplot_helper.exe']
        data = [os.path.join(local_path,x) for x in data]               
        config['data_files'].extend( [(os.path.join(parent_path,'gplt'),data)])    
    return config

