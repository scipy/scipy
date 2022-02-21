
def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('transform', parent_package, top_path)

    config.add_data_dir('tests')

    config.add_data_files('_rotation.pyi')
    config.add_extension('_rotation',
                         sources=['_rotation.c'])

    return config
