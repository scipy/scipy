
def configuration(parent_package=''):

    # The following three lines constitute minimal contents
    # of configuration(..) that is suitable for pure Python
    # packages.
    package = 'xxx'
    from scipy_distutils.misc_util import default_config_dict,dot_join
    config = default_config_dict(package,parent_package)

    # xxx contains pure Python sub-package yyy
    config['packages'].append(dot_join(parent_package,package,'yyy'))
    # And don't forget adding sub-packages of yyy:
    config['packages'].append(dot_join(parent_package,package,'yyy.tests'))

    # xxx has extension module ..
    # xxx has f2py generated extension module ..
    # xxx has swig generated extension module ..

    return config

if __name__ == '__main__':
    from scipy_distutils.core import setup
    setup(**configuration())
