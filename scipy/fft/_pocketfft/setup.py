from setuptools import setup, Extension
import sys


class _deferred_pybind11_include(object):
    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)


# TODO: Compiler arguments should be determined by the compiler used

include_dirs = [_deferred_pybind11_include(True), _deferred_pybind11_include()]
#extra_compile_args = ['--std=c++11', '-march=native', '-O3']
extra_compile_args = []
python_module_link_args = []


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('_pocketfft', parent_package, top_path)
    config.add_extension('pypocketfft',
                         sources=['pypocketfft.cxx'],
                         depends=['pocketfft.h'],
                         include_dirs=include_dirs,
                         language='c++',
                         extra_compile_args=extra_compile_args,
                         extra_link_args=python_module_link_args)

    config.add_data_files('LICENSE.md')
    config.add_data_dir('tests')
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
