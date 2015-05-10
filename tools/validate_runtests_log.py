#!/usr/bin/env python
"""
Take the test runner log output from the stdin, looking for
the magic line nose runner prints when the test run was successful.

In an ideal world, this should be done directly in runtests.py  using the
nose API, some failure modes are fooling nose to terminate the python process
with zero exit code, see, eg, https://github.com/scipy/scipy/issues/4736
In short, lapack's xerbla can terminate the process with a fortran level STOP
command, which (i) aborts the py process so that runtests.py does not finish,
and (ii) the exit code is implementation-defined.

Also check that the number of tests run is larger than some baseline number
(taken from the state of the master branch at some random point in time.)
This probably could/should be made less brittle.

"""
from __future__ import print_function

import sys
import re


if __name__ == "__main__":
    # full or fast test suite?
    try:
        testmode = sys.argv[1]
        if testmode not in ('fast', 'full'):
            raise IndexError
    except IndexError:
        raise ValueError("Usage: validate.py {full|fast} < logfile.")

    # fetch the expected number of tests
    # these numbers are for 6abad09
    # XXX: this should probably track the commit hash or commit date
    expected_size = {'full': 19055,
                     'fast': 17738}

    # read in the log, parse for the success line.
    r = re.compile("Ran (?P<num_tests>\d+) tests in (?P<time>\d+\S+)")

    found_it = False
    for line in sys.stdin:
        m = r.search(line)
        if m:
            found_it = True
            break

    if found_it:
        # now check that the number of tests run is reasonable
        expected = expected_size[testmode]
        actual = int(m.group('num_tests'))
        if actual < expected:
            print("Too few tests: expected %s, run %s" % (expected, actual))
            sys.exit(1)
        else:
            sys.exit(0)
    else:
        print('Test runner validation errored: did the run really finish?')
        sys.exit(-1)
