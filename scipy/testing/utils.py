"""Simple testing utilities
"""

__all__ = ['set_local_path', 'restore_path', 'measure', 'info', 'warn', 'error']

import os
import sys

from numpy.distutils.misc_util import yellow_text, red_text
from numpy.testing.utils import jiffies

def set_local_path():
    """ Prepend local directory to sys.path.

    The caller is responsible for removing this path by using

      restore_path()
    """
    f = sys._getframe(1)
    if f.f_locals['__name__']=='__main__':
        testfile = sys.argv[0]
    else:
        testfile = f.f_locals['__file__']
    local_path = os.path.normpath(os.path.dirname(os.path.abspath(testfile)))
    sys.path.insert(0,local_path)
    return

def restore_path():
    del sys.path[0]
    return

def measure(code_str,times=1,label=None):
    """ Return elapsed time for executing code_str in the
    namespace of the caller for given times.
    """
    frame = sys._getframe(1)
    locs,globs = frame.f_locals,frame.f_globals

    code = compile(code_str,
                   'Test name: %s ' % label,
                   'exec')
    i = 0
    elapsed = jiffies()
    while i<times:
        i += 1
        exec code in globs,locs
    elapsed = jiffies() - elapsed
    return 0.01*elapsed

def info(message):
    print >> sys.stdout, message
    sys.stdout.flush()

def warn(message):
    print >> sys.stderr,yellow_text('Warning: %s' % (message))
    sys.stderr.flush()

def error(message):
    print >> sys.stderr,red_text('Error: %s' % (message))
    sys.stderr.flush()
