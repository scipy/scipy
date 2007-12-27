"""
Cleaned-up NumpyTestCase class, to be moved into numpy.testing later...
"""

import unittest

class NumpyTestCase2(unittest.TestCase):
    """New implementation of NumpyTestCase that is nose-friendly.

    It is slated for replacing the base one eventually, but it will be first
    used in scipy's testing.  Once everything is OK there, it will be used
    into all of numpy."""

    def measure(self,code_str,times=1):
        """ Return elapsed time for executing code_str in the
        namespace of the caller for given times.
        """
        frame = get_frame(1)
        locs,globs = frame.f_locals,frame.f_globals
        code = compile(code_str,
                       'NumpyTestCase runner for '+self.__class__.__name__,
                       'exec')
        i = 0
        elapsed = jiffies()
        while i<times:
            i += 1
            exec code in globs,locs
        elapsed = jiffies() - elapsed
        return 0.01*elapsed

    def warn(self, message):
        from numpy.distutils.misc_util import yellow_text
        print>>sys.stderr,yellow_text('Warning: %s' % (message))
        sys.stderr.flush()

    def info(self, message):
        print>>sys.stdout, message
        sys.stdout.flush()
