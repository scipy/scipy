
def configuration(parent_package='io',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('matlab', parent_package, top_path)
    config.add_extension('_streams', sources=['_streams.c'])
    config.add_extension('_mio_utils', sources=['_mio_utils.c'])
    config.add_extension('_mio5_utils', sources=['_mio5_utils.c'])
    config.add_data_dir('tests')
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
