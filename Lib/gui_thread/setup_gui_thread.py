from scipy_distutils.misc_util import get_path, default_config_dict

def configuration(parent_package=''):
    if parent_package:
        parent_package += '.'
    local_path = get_path(__name__)

    config = default_config_dict()
    config['packages'].append(parent_package+'gui_thread')
    #config['packages'].append(parent_package+'gui_thread.tests') 

    #I think this is needed
    config['package_dir']['gui_thread'] = local_path
    
    return config
