# Last Change: Wed Mar 05 03:00 PM 2008 J
# vim:syntax=python
from os.path import join

from numpy.distutils.misc_util import get_numpy_include_dirs
from numscons import GetNumpyEnvironment

env = GetNumpyEnvironment(ARGUMENTS)

env.AppendUnique(CPPPATH = get_numpy_include_dirs())
env.NumpyPythonExtension('numpyio', source = 'numpyiomodule.c')
