"""
SciPy --- A scientific computing package for Python
===================================================

Available subpackages
---------------------

"""

"""
SciPy has levels
Level 0 -- Numeric and core routines in scipy_base

Level 1 -- Level 0 + fft, special, linalg, stats
Level 1a -- Core routines which depend on Level 1.

Level 2 -- plotting interface.

Level 3
Packages which define own functions plus depend on Levels 0-2.

Level 0, 1, 2 should be imported in order and then other levels imported
  as available.
"""

# Level 0
# modules to import under the scipy namespace
from scipy_version import scipy_version as __version__
from scipy_base import *

from helpmod import *

try:
    help
except NameError:
    help = info

#-------- doc string hooks --------#

def _level_docs(module=None,_cache=[]):
    if module is None:
        return _cache
    try:
        title = module.__dict__['__doc_title__']
    except KeyError:
        try:
            title = module.__dict__['__doc__']
            title = title.lstrip().split('\n',1)[0]
        except KeyError:
            title = ' N/A'
    _cache.append((module.__name__,title))

def _pkg_titles():
    level_docs = _level_docs()
    lengths = [len(name) for (name,title) in level_docs]
    max_length = max(lengths)
    lines = []
    for (name,title) in level_docs:
        w = max_length - len(name)
        lines.append('%s%s --- %s' % (name, w*' ', title))
    return '\n'.join(lines)

#----------------------------------#

# Level 1
# these modules will just be imported (not subsumed)
special = ppimport('special');  _level_docs(special)
io = ppimport('io');            _level_docs(io)
linalg = ppimport('linalg');    _level_docs(linalg)
stats = ppimport('stats');      _level_docs(stats)
fftpack = ppimport('fftpack');  _level_docs(fftpack)

for n in ['fft','fftn','fft2','ifft','ifft2','ifftn','fftshift',
          'ifftshift','fftfreq']:
    exec '%s = ppimport_attr(fftpack,%s)' % (n,`n`)

for n in ['mean','median','std','cov','corrcoef']:
    exec '%s = ppimport_attr(stats,%s)' % (n,`n`)

for n in ['isinf','isfinite','isnan']:
    exec '%s = ppimport_attr(special,%s)' % (n,`n`)

# Functions to be subsumed that need Level 0 and Level 1
from common import *
from pilutil import *

# Level 2
optimize = ppimport('optimize');       _level_docs(optimize)
integrate = ppimport('integrate');     _level_docs(integrate)
signal = ppimport('signal');           _level_docs(signal)

interpolate = ppimport('interpolate'); _level_docs(interpolate)
cow = ppimport('cow');                 _level_docs(cow)
ga = ppimport('ga');                   _level_docs(ga)
cluster = ppimport('cluster');         _level_docs(cluster)
weave = ppimport('weave');             _level_docs(weave)

# Level 3
xplt = ppimport('xplt');               _level_docs(xplt)
gplt = ppimport('gplt');               _level_docs(gplt)
plt = ppimport('plt');                 _level_docs(plt)

#----- update doc string -------#

__doc__ += _pkg_titles()
            
#---- testing ----#

def test(level=1):
    """ From this top level, there are possibly many many tests.
        Test only the quick tests by default.
    """
    import unittest
    runner = unittest.TextTestRunner()
    runner.run(test_suite(level))
    return runner

def test_all(level=10):
    test(level)
    
def test_suite(level = 1):
    import scipy_test.testing
    import scipy
    ignore = ['xplt','plt','gplt','gui_thread','sparse','scipy_version']
    suites = [scipy_test.testing.harvest_test_suites(scipy,ignore,level=level)]
    import unittest
    return unittest.TestSuite(suites)
