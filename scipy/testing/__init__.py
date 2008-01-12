"""Common test support for all scipy test scripts.

This single module should provide all the common functionality for scipy tests
in a single location, so that test script can just import it and work right
away.
"""

import unittest
from unittest import TestCase

try:
    import nose
except ImportError:
    pass

import decorators as dec
from numpy.testing.utils import *
from utils import *
