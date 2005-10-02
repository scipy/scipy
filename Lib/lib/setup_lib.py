#!/usr/bin/env python

def configuration(parent_package='',parent_path=None):
    from scipy.distutils.misc_util import dot_join, get_path,\
         default_config_dict, merge_config_dicts, get_subpackages

    package = 'lib'
    config = default_config_dict(package,parent_package)
    local_path = get_path(__name__,parent_path)

    config_list = [config]
    config_list += get_subpackages(local_path,
                                   parent=config['name'],
                                   parent_path=parent_path)

    config_dict = merge_config_dicts(config_list)

    config = merge_config_dicts(config_list)
    return config

if __name__ == '__main__':
    from scipy.distutils.core import setup

    setup(**configuration(parent_path=''))
