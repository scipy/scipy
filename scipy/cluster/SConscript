# Last Change: Mon Nov 03 07:00 PM 2008 J
# vim:syntax=python
from os.path import join

from numscons import GetNumpyEnvironment

env = GetNumpyEnvironment(ARGUMENTS)

env.NumpyPythonExtension('_hierarchy_wrap',
                         source = [join('src', 'hierarchy_wrap.c'),
                                   join('src', 'hierarchy.c')])

env.NumpyPythonExtension('_vq',
                         source = [join('src', 'vq_module.c'),
                                   join('src', 'vq.c')])
