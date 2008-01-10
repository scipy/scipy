============
 Nose tests
============

SciPy will use, and require nose for all testing.  It will be possible for
users to use all of scipy without nose on their systems, but *not* to run the
test suite.

Labels
======

We are going to use labels for tests, instead of the old level system.  For
now, we'll use 'slow' as a label for slow-running tests, and will make more if
needed in the future for speed considerations.

Labels will be implemented as one decorator per test.  All pre-defined labels
are implemented in numpy.testing.decorators

This would run all tests, since nose by default discovers everything::
  nosetests -s  test_foo.py

this::
  nosetests -s -A slow test_foo.py

will only run tests labeled 'slow' (for which there's a decorator) and this::
  nosetests -s -A "not slow" test_foo.py

will *exclude* all slow tests from a run.

The scipy.test() function call will expose the above with convenient syntax.

Initially we only have the @slow decorator, later we may provide new ones as
the need for them arises in actual use.


Benchmark tests
===============

Routines that don't actually test correctness but instead do performance
benchmarking will live in a benchmarks/ directory next to the tests/ directory
of each module.  There will be a scipy.benchmark() call that does benchmarking,
similar to scipy.test() but separate from it.

Scipy test

For each package, there will be a function that takes level arguments,
and performs tests per level

import scipy.mypackage
scipy.mypackage.test() # fast, unlabeled tests 
scipy.mypackage.test('full') # unlabeled, not slow, not villfail
scipy.mypackage.test('slow') # just slow tests
scipy.mypackage.test('bench') # just benchmarks
scipy.mypackage.test(None) # all possible tests, including benchmarks
scipy.mypackage.test(doctests=True) # fast tests, with doctests

At the base level, scipy.test(*args) collects the test suite from each
package, and runs it, with *args as above.

scipy.mypackage.test()

Runs all plausible tests in this package, and package test directory,
and runs any subpackage tests such as
scipy.mypackage.mysubpackage.test()


