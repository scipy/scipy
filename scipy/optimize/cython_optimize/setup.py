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
    Extension("newton_example", ["newton_example.pyx"],
              include_dirs=[np.get_include()])
)

setup(
    name='cython optimize api',
    ext_modules=cythonize(EXTENSIONS),
    include_dirs=[np.get_include()]
)
