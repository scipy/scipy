""" Test refcounting and behavior of SCXX.
"""
import unittest
import time
import os,sys

from numpy.testing import *
from test_scxx_object import *
from test_scxx_sequence import *
from test_scxx_dict import *


if __name__ == "__main__":
    nose.run(argv=['', __file__])
