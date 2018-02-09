

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    config = Configuration('_ivp', parent_package, top_path)

    config.add_extension('base',
                         sources=['base.c'])

    config.add_extension('bdf',
                         sources=['bdf.c'])

    config.add_extension('common',
                         sources=['common.c'])
    
    config.add_extension('ivp',
                         sources=['ivp.c'])

    config.add_extension('lsoda',
                         sources=['lsoda.c'])

    config.add_extension('radau',
                         sources=['radau.c'])

    config.add_extension('rk',
                         sources=['rk.c'])

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
