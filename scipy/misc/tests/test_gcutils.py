""" Test for check_refs_for context manager and gc utilities
"""
import gc

from scipy.misc.gcutils import set_gc_state, gc_state, check_refs_for, ReferenceError

from nose.tools import assert_equal, raises

def test_set_gc_state():
    gc_status = gc.isenabled()
    try:
        for state in (True, False):
            gc.enable()
            set_gc_state(state)
            assert_equal(gc.isenabled(), state)
            gc.disable()
            set_gc_state(state)
            assert_equal(gc.isenabled(), state)
    finally:
        if gc_status:
            gc.enable()


def test_gc_state():
    # Test gc_state context manager
    gc_status = gc.isenabled()
    try:
        for pre_state in (True, False):
            set_gc_state(pre_state)
            for with_state in (True, False):
                # Check the gc state is with_state in with block
                with gc_state(with_state):
                    assert_equal(gc.isenabled(), with_state)
                # And returns to previous state outside block
                assert_equal(gc.isenabled(), pre_state)
                # Even if the gc state is set explicitly within the block
                with gc_state(with_state):
                    assert_equal(gc.isenabled(), with_state)
                    set_gc_state(not with_state)
                assert_equal(gc.isenabled(), pre_state)
    finally:
        if gc_status:
            gc.enable()


def test_check_refs_for():
    # Ordinary use
    class C(object):
        def __init__(self, arg0, arg1, name='myname'):
            self.name = name
    for gc_current in (True, False):
        with gc_state(gc_current):
            # We are deleting from with-block context, so that's OK
            with check_refs_for(C, 0, 2, 'another name') as c:
                assert_equal(c.name, 'another name')
                del c
            # Or not using the thing in with-block context, also OK
            with check_refs_for(C, 0, 2, name='third name'):
                pass
            assert_equal(gc.isenabled(), gc_current)


@raises(ReferenceError)
def test_check_refs_for_nodel():
    class C(object): pass
    # Need to delete after using if in with-block context
    with check_refs_for(C) as c:
        pass


@raises(ReferenceError)
def test_check_refs_for_circular():
    class C(object):
        def __init__(self):
            self._circular = self
    # Circular reference, no automatic garbage collection
    with check_refs_for(C) as c:
        del c


@raises(ReferenceError)
def test_check_refs_for_circular2():
    class C(object):
        def __init__(self):
            self._circular = self
    # Still circular reference, no automatic garbage collection
    with check_refs_for(C):
        pass
