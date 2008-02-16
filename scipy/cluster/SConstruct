# Last Change: Thu Oct 18 09:00 PM 2007 J
# vim:syntax=python
from os.path import join

from numpy.distutils.misc_util import get_numpy_include_dirs
from numscons import GetNumpyEnvironment

env = GetNumpyEnvironment(ARGUMENTS)

env.AppendUnique(CPPPATH = get_numpy_include_dirs())
env.NumpyPythonExtension('_vq', source = [join('src', 'vq_module.c'),
                                          join('src', 'vq.c')])
