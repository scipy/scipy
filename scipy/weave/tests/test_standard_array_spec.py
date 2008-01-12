
from numpy import *
from scipy.testing import *

from scipy.weave import standard_array_spec

def remove_whitespace(in_str):
    out = in_str.replace(" ","")
    out = out.replace("\t","")
    out = out.replace("\n","")
    return out

def print_assert_equal(test_string,actual,desired):
    """this should probably be in scipy_test.testing
    """
    import pprint
    try:
        assert(actual == desired)
    except AssertionError:
        import cStringIO
        msg = cStringIO.StringIO()
        msg.write(test_string)
        msg.write(' failed\nACTUAL: \n')
        pprint.pprint(actual,msg)
        msg.write('DESIRED: \n')
        pprint.pprint(desired,msg)
        raise AssertionError, msg.getvalue()

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
    unittest.main()
