#!/usr/bin/env python
"""Install file for example on how to use Pyrex with Numpy.

For more details, see:
http://www.scipy.org/Cookbook/Pyrex_and_NumPy
http://www.scipy.org/Cookbook/ArrayStruct_and_Pyrex
"""

from distutils.core import setup
from distutils.extension import Extension

# Make this usable by people who don't have pyrex installed (I've committed
# the generated C sources to SVN).
try:
    from Pyrex.Distutils import build_ext
    has_pyrex = True
except ImportError:
    has_pyrex = False
import numpy

# Define a pyrex-based extension module, using the generated sources if pyrex
# is not available.
if has_pyrex:
    pyx_sources = ['cloess.pyx']
    cmdclass    = {'build_ext': build_ext}
else:
    pyx_sources = ['cloess.c']
    cmdclass    = {}

f_sources = ['loessf.f', 'linpack_lite.f']
c_sources = ['loess.c', 'loessc.c']


pyx_ext = Extension('cloess',
                    pyx_sources + c_sources + f_sources,
                    include_dirs = [numpy.get_include()])

# Call the routine which does the real work
setup(name        = 'cloess',
      description = 'Small example on using Pyrex to write a Numpy extension',
      url         = 'http://www.scipy.org/Cookbook/Pyrex_and_NumPy',
      ext_modules = [pyx_ext],
      cmdclass    = cmdclass,
      )
