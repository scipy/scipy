from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy as np

EXTENSIONS = []

# newton_zero
EXTENSIONS.append(
    Extension("zeros", ["zeros.pyx"],
              include_dirs=[np.get_include()])
)

# ivexample
EXTENSIONS.append(
    Extension("tests.zeros_examples", ["tests/zeros_examples.pyx"],
              include_dirs=[np.get_include()])
)

# zeros structs
EXTENSIONS.append(
    Extension("zeros_struct", ["zeros_struct.pyx"],
              include_dirs=[np.get_include()])
)
EXTENSIONS.append(
    Extension("tests.zeros_struct_examples", ["tests/zeros_struct_examples.pyx"],
              include_dirs=[np.get_include()])
)

setup(
    name='cython optimize api',
    ext_modules=cythonize(EXTENSIONS),
    include_dirs=[np.get_include()]
)
