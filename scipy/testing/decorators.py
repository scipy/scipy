"""Decorators for labeling test objects."""

def slow(t):
    """Labels a test as 'slow'.

    The exact definition of a slow test is obviously both subjective and
    hardware-dependent, but in general any individual test that requires more
    than a second or two should be labeled as slow (the whole suite consits of
    thousands of tests, so even a second is significant)."""
    
    t.slow = True
    return t

def bench(t):
    ''' Labels a test as a benchmark.

    Benchmark tests are often slow, and intended to test timings
    between different algorithms rather than validity of the result. '''
    
    t.bench = True
    return t

def willfail(t):
    ''' Labels test as known failure

    This label allows the tester to deselect the test in standard cases '''
    t.willfail = True
    return t

def is_nosetest(tf):
    ''' Signals to nose that this function is or is not a test

    e.g
    >>> @is_nosetest(False)
    >>> def func_with_test_in_name(arg1, arg2): pass
    ...
    >>> 
    '''
    def set_test(t):
        t.__test__ = tf
        return t
    return set_test
