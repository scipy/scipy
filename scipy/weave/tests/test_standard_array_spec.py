from numpy import *
from numpy.testing import *

from scipy.weave import standard_array_spec

def remove_whitespace(in_str):
    out = in_str.replace(" ","")
    out = out.replace("\t","")
    out = out.replace("\n","")
    return out

class TestArrayConverter(TestCase):
    def test_type_match_string(self):
        s = standard_array_spec.array_converter()
        assert( not s.type_match('string') )
    def test_type_match_int(self):
        s = standard_array_spec.array_converter()
        assert(not s.type_match(5))
    def test_type_match_array(self):
        s = standard_array_spec.array_converter()
        assert(s.type_match(arange(4)))

if __name__ == "__main__":
    nose.run(argv=['', __file__])
