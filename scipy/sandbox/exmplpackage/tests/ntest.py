"""Simple testing utilities, to eventually be put into numpy.testing
"""

import sys

from numpy.distutils.misc_util import yellow_text
from numpy.testing.utils import jiffies

def measure(code_str,times=1,test_name=None):
    """ Return elapsed time for executing code_str in the
    namespace of the caller for given times.
    """
    frame = sys.get_frame(1)
    locs,globs = frame.f_locals,frame.f_globals

    code = compile(code_str,
                   'Test name: %s '+test_name,
                   'exec')
    i = 0
    elapsed = jiffies()
    while i<times:
        i += 1
        exec code in globs,locs
    elapsed = jiffies() - elapsed
    return 0.01*elapsed

def warn(message):
    print >> sys.stderr,yellow_text('Warning: %s' % (message))
    sys.stderr.flush()

def info(message):
    print >> sys.stdout, message
    sys.stdout.flush()
