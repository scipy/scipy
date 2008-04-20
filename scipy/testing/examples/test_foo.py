"""Example SciPy test module, using nose features.

For information on nose, see:

  http://somethingaboutorange.com/mrl/projects/nose

To run it in its simplest form::
  nosetests test_foo.py


Labels will be implemented as one decorator per test.  All pre-defined labels
are implemented in numpy.testing.decorators

This would run all tests, since nose by default discovers everything::
  nosetests -s  test_foo.py

this::
  nosetests -s -A slow test_foo.py

will only run tests labeled 'slow' (for which there's a decorator) and this::
  nosetests -s -A "not slow" test_foo.py

will *exclude* all slow tests from a run.
"""

# This single import statement should provide all the common functionality for
# scipy tests in a single location.
from scipy.testing import *

def setup():
    """Module-level setup"""
    print 'doing setup'

def teardown():
    """Module-level teardown"""
    print 'doing teardown'


class ClassicTest(TestCase):
    """A regular unittest, with the extra Numpy features."""
    def test_1(self):
        print 'First test'

    def test_2(self):
        print 'Second test'

    # This is how to tag a test as slow
    @dec.slow
    def test_big(self,level=5):
        print 'Big, slow test'


def test_simplefunction():
    """A simple test function."""
    assert True

def check_even(n, nn):
    """A check function to be used in a test generator."""
    assert n % 2 == 0 or nn % 2 == 0

# Test generators are best left without docstring, so nose (in verbose mode)
# shows the actual call they produce automatically, including arguments.
def test_evens():
    for i in range(0,4,2):
        yield check_even, i, i*3

def setup_func():
    """A trivial setup function."""
    print "In setup_func"

def teardown_func():
    """A trivial teardown function."""
    print "In teardown_func"

@nose.with_setup(setup_func, teardown_func)
def test_with_extras():
    """This test uses the setup/teardown functions."""
    print "  In test_with_extras"

# Setup and teardown functions may be used with test generators. The setup and
# teardown attributes must be attached to the generator function:
@nose.with_setup(setup_func, teardown_func)
def test_generator():
    for i in range(6,10,2):
        yield check_even, i, i*3

@dec.slow
def test_nasty():
    "A difficult test that takes a long time..."
    print '*** nasty slow test ***'


def test_time():
    "A simple test that times things"
    x = 1
    time = measure("x+1",times=100,label='test_time')
    info('Time taken: %s' % time)

def test_warn():
    "A simple test that prints a warning."
    warn('Bad things are happening...')

def test_error():
    "A simple test that prints an error message."
    error('Really bad things are happening...')
