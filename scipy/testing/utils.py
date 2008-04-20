"""Simple testing utilities
"""

__all__ = ['set_local_path', 'restore_path', 'measure', 'info', 'warn',\
           'error', 'decorate_methods']

import os
import sys
import re
from inspect import isfunction

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

def decorate_methods(cls, decorator, testmatch=None):
    ''' Apply decorator to all methods in class matching testmatch

    Parameters
    ----------
    cls : class
        Class to decorate methods for
    decorator : function
        Decorator to apply to methods
    testmatch : compiled regexp or string to compile to regexp
        Decorators are applied if testmatch.search(methodname)
        is not None.  Default value is
        re.compile(r'(?:^|[\\b_\\.%s-])[Tt]est' % os.sep)
        (the default for nose)
    '''
    if testmatch is None:
        testmatch = re.compile(r'(?:^|[\\b_\\.%s-])[Tt]est' % os.sep)
    else:
        testmatch = re.compile(testmatch)
    cls_attr = cls.__dict__
    methods = filter(isfunction, cls_attr.values())
    for function in methods:
        try:
            if hasattr(function, 'compat_func_name'):
                funcname = function.compat_func_name
            else:
                funcname = function.__name__
        except AttributeError:
            # not a function
            continue
        if testmatch.search(funcname) and not funcname.startswith('_'):
            setattr(cls, funcname, decorator(function))
    return
