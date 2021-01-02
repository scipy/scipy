'''Build a test extension using Boost headers.'''


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    from scipy._build_utils.boostinator import get_include_dir
    config = Configuration('boost_test_build', parent_package, top_path)
    config.add_library(
        'boost_test_build',
        sources=['boost_test_build.cpp'],
        include_dirs=[str(get_include_dir())],
        language='c++')
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
