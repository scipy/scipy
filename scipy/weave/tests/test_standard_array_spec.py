from __future__ import absolute_import, print_function

from numpy import arange
from numpy.testing import TestCase, assert_, run_module_suite

from scipy.weave import standard_array_spec


class TestArrayConverter(TestCase):

    def test_type_match_string(self):
        s = standard_array_spec.array_converter()
        assert_(not s.type_match('string'))

    def test_type_match_int(self):
        s = standard_array_spec.array_converter()
        assert_(not s.type_match(5))

    def test_type_match_array(self):
        s = standard_array_spec.array_converter()
        assert_(s.type_match(arange(4)))


if __name__ == "__main__":
    run_module_suite()
