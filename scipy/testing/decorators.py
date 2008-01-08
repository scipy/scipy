"""Decorators for labeling test objects."""

''' String expression which is true for all undecorated tests '''
undecorated_def = 'not slow'

def slow(t):
    """Labels a test as 'slow'.

    The exact definition of a slow test is obviously both subjective and
    hardware-dependent, but in general any individual test that requires more
    than a second or two should be labeled as slow (the whole suite consits of
    thousands of tests, so even a second is significant)."""
    
    t.slow = True
    return t

