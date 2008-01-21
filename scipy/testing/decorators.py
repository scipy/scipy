"""Decorators for labeling test objects."""

try:
    import nose
except ImportError:
    pass

def slow(t):
    """Labels a test as 'slow'.

    The exact definition of a slow test is obviously both subjective and
    hardware-dependent, but in general any individual test that requires more
    than a second or two should be labeled as slow (the whole suite consits of
    thousands of tests, so even a second is significant)."""
    
    t.slow = True
    return t

def willfail(t):
    ''' Labels test as known failure

    This label allows the tester to deselect the test in standard cases '''
    t.willfail = True
    return t

def setastest(tf=True):
    ''' Signals to nose that this function is or is not a test
    
    e.g
    >>> @setastest(False)
    >>> def func_with_test_in_name(arg1, arg2): pass
    ...
    >>>
    
    Note that this decorator cannot use the nose namespace, because it
    can be called from a non-test.
    '''
    def set_test(t):
        t.__test__ = tf
        return t
    return set_test

def skipif(skip_condition, msg=None):
    ''' Make function raise SkipTest exception if skip_condition is true

    Parameters
    ---------
    skip_condition : bool
        Flag to determine whether to skip test (True) or not (False)
    msg : string
        Message to give on raising a SkipTest exception
    '''
    if msg is None:
        msg = 'Test skipped due to test condition (see code)'
    def skip_decorator(f):
        def skipper(*args, **kwargs):
            if skip_condition:
                raise nose.SkipTest, msg
            else:
                return f(*args, **kwargs)
        return nose.tools.make_decorator(f)(skipper)
    return skip_decorator
