from distutils.core import setup
from Cython.Build import cythonize

setup(
  name='cython optimize api',
  ext_modules=cythonize(["newton_example.pyx", "zeros.pyx"]),
)