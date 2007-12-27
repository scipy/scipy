"""Common test support for all scipy test scripts.

This single module should provide all the common functionality for scipy tests
in a single location, so that test script can just import it and work right
away.
"""

from unittest import TestCase

import nose

# These two modules will need to be later put in the right places...
import decorators as dec

from ntest import *
