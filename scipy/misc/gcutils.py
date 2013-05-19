""" Module for testing automatic garbage collection of objects
"""
import weakref
import gc

from contextlib import contextmanager

__all__ = ['set_gc_state', 'gc_state', 'check_refs_for']

class ReferenceError(Exception):
    pass


def set_gc_state(state):
    """ Set status of garbage collector """
    if gc.isenabled() == state:
        return
    if state:
        gc.enable()
    else:
        gc.disable()


@contextmanager
def gc_state(state):
    """ Context manager to set state of garbage collector to `state`

    Parameters
    ----------
    state : bool
        True for gc enabled, False for disabled

    Examples
    --------
    >>> with gc_state(False):
    ...     assert not gc.isenabled()
    >>> with gc_state(True):
    ...     assert gc.isenabled()
    """
    orig_state = gc.isenabled()
    set_gc_state(state)
    yield
    set_gc_state(orig_state)


@contextmanager
def check_refs_for(func, *args, **kwargs):
    """ Context manager to check for remaining references after deletion

    Parameters
    ----------
    func : callable
        Callable to create object to check
    \*args : sequence
        positional arguments to `func` in order to create object to check
    \*\*kwargs : dict
        keyword arguments to `func` in order to create object to check

    Examples
    --------
    >>> class C(object): pass
    >>> with check_refs_for(C) as c:
    ...     # do something
    ...     del c

    >>> class C(object):
    ...     def __init__(self):
    ...         self._circular = self # Make circular reference
    >>> with check_refs_for(C) as c: #doctest: +IGNORE_EXCEPTION_DETAIL
    ...     # do something
    ...     del c
    Traceback (most recent call last):
        ...
    ReferenceError: Remaining reference(s) to object
    """
    with gc_state(False):
        obj = func(*args, **kwargs)
        ref = weakref.ref(obj)
        yield obj
        del obj
        if ref() != None:
            raise ReferenceError("Remaining reference(s) to object")
