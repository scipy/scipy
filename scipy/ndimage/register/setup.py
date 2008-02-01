
#!/usr/bin/env python

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('register', parent_package, top_path)

    config.add_extension('_register',
                         sources=['Register_EXT.c',
                                  'Register_IMPL.c']
    )

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())


