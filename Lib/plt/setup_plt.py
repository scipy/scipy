import os
from scipy_distutils.misc_util import get_path, default_config_dict

def configuration(parent_package=''):
    parent_path = parent_package
    if parent_package:
        parent_package += '.'
    local_path = get_path(__name__)

    config = default_config_dict()
    config['packages'].append(parent_package+'plt')
    #config['packages'].append(parent_package+'plt.tests') 

    data = ['lena.dat', 'colormaps.dir','colormaps.dat' ]
    data = [os.path.join(local_path,x) for x in data]               
    config['data_files'].extend( [(os.path.join(parent_path,'plt'),data)])
  
    return config
