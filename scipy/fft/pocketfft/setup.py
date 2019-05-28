from setuptools import setup, Extension
import sys


class _deferred_pybind11_include(object):
    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)


include_dirs = ['./', _deferred_pybind11_include(True),
                _deferred_pybind11_include()]
extra_compile_args = ['--std=c++11', '-march=native', '-O3']
python_module_link_args = []

if sys.platform == 'darwin':
    import distutils.sysconfig
    extra_compile_args += ['--stdlib=libc++', '-mmacosx-version-min=10.9', '-openmp']
    vars = distutils.sysconfig.get_config_vars()
    vars['LDSHARED'] = vars['LDSHARED'].replace('-bundle', '')
    python_module_link_args+=['-bundle', '-openmp']
else:
    extra_compile_args += ['-DPOCKETFFT_OPENMP', '-fopenmp', '-Wfatal-errors', '-Wfloat-conversion' ,'-Wsign-conversion', '-Wconversion' ,'-W', '-Wall', '-Wstrict-aliasing=2', '-Wwrite-strings', '-Wredundant-decls', '-Woverloaded-virtual', '-Wcast-qual', '-Wcast-align', '-Wpointer-arith']
    python_module_link_args += ['-march=native', '-Wl,-rpath,$ORIGIN', '-fopenmp']

# if you don't want debugging info, add "-s" to python_module_link_args


def get_extension_modules():
    return [Extension('pypocketfft',
                      language='c++',
                      sources=['pypocketfft.cc'],
                      depends=['pocketfft_hdronly.h'],
                      include_dirs=include_dirs,
                      extra_compile_args=extra_compile_args,
                      extra_link_args=python_module_link_args)]


setup(name='pypocketfft',
      version='0.0.1',
      description='Python interface for pocketfft',
      include_package_data=True,
      author='Martin Reinecke',
      author_email='martin@mpa-garching.mpg.de',
      packages=[],
      setup_requires=['numpy>=1.15.0', 'pybind11>=2.2.4'],
      ext_modules=get_extension_modules(),
      install_requires=['numpy>=1.15.0', 'pybind11>=2.2.4']
      )
