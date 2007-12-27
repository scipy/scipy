"""Simple testing utilities, to eventually be put into numpy.testing
"""

import sys

from numpy.distutils.misc_util import yellow_text, red_text
from numpy.testing.utils import jiffies

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
