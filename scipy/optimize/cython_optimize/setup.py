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

# structs
EXTENSIONS.append(
    Extension("test.test_zeros_struct", ["tests/test_zeros_struct.pyx"],
              include_dirs=[np.get_include()])
)

setup(
    name='cython optimize api',
    ext_modules=cythonize(EXTENSIONS),
    include_dirs=[np.get_include()]
)
