from numpy.testing import *

from scipy.weave import slice_handler
from scipy.weave.slice_handler import indexed_array_pattern
from scipy.weave.ast_tools import *

class TestBuildSliceAtom(TestCase):
    def generic_check(self,slice_vars,desired):
        pos = slice_vars['pos']
        ast_list = slice_handler.build_slice_atom(slice_vars,pos)
        actual = ast_to_string(ast_list)
        assert_equal(actual,desired)
    def test_exclusive_end(self):
        slice_vars = {'begin':'1', 'end':'2', 'step':'_stp',
                      'single_index':'_index','pos':0}
        desired = 'slice(1,2-1)'
        self.generic_check(slice_vars,desired)

class TestSlice(TestCase):

    def generic_check(self,suite_string,desired):
        import parser
        ast_tuple = parser.suite(suite_string).totuple()
        found, data = find_first_pattern(ast_tuple,indexed_array_pattern)
        subscript = data['subscript_list'][1] #[0] is symbol, [1] is the supscript
        actual = slice_handler.slice_ast_to_dict(subscript)
        assert_equal(actual,desired,suite_string)

    def test_empty_2_slice(self):
        """match slice from a[:]"""
        test ="a[:]"
        desired = {'begin':'_beg', 'end':'_end', 'step':'_stp',
                   'single_index':'_index'}
        self.generic_check(test,desired)
    def test_begin_2_slice(self):
        """match slice from a[1:]"""
        test ="a[1:]"
        desired = {'begin':'1', 'end':'_end', 'step':'_stp',
                   'single_index':'_index'}
        self.generic_check(test,desired)
    def test_end_2_slice(self):
        """match slice from a[:2]"""
        test ="a[:2]"
        desired = {'begin':'_beg', 'end':'2', 'step':'_stp',
                   'single_index':'_index'}
        self.generic_check(test,desired)
    def test_begin_end_2_slice(self):
        """match slice from a[1:2]"""
        test ="a[1:2]"
        desired = {'begin':'1', 'end':'2', 'step':'_stp',
                   'single_index':'_index'}
        self.generic_check(test,desired)
    def test_empty_3_slice(self):
        """match slice from a[::]"""
        test ="a[::]"
        desired = {'begin':'_beg', 'end':'_end', 'step':'_stp',
                   'single_index':'_index'}
        self.generic_check(test,desired)
    def test_begin_3_slice(self):
        """match slice from a[1::]"""
        test ="a[1::]"
        desired = {'begin':'1', 'end':'_end', 'step':'_stp',
                   'single_index':'_index'}
        self.generic_check(test,desired)
    def test_end_3_slice(self):
        """match slice from a[:2:]"""
        test ="a[:2:]"
        desired = {'begin':'_beg', 'end':'2', 'step':'_stp',
                   'single_index':'_index'}
        self.generic_check(test,desired)
    def test_stp3_slice(self):
        """match slice from a[::3]"""
        test ="a[::3]"
        desired = {'begin':'_beg', 'end':'_end', 'step':'3',
                   'single_index':'_index'}
        self.generic_check(test,desired)
    def test_begin_end_3_slice(self):
        """match slice from a[1:2:]"""
        test ="a[1:2:]"
        desired = {'begin':'1', 'end':'2','step':'_stp',
                   'single_index':'_index'}
        self.generic_check(test,desired)
    def test_begin_step_3_slice(self):
        """match slice from a[1::3]"""
        test ="a[1::3]"
        desired = {'begin':'1', 'end':'_end','step':'3',
                   'single_index':'_index'}
        self.generic_check(test,desired)
    def test_end_step_3_slice(self):
        """match slice from a[:2:3]"""
        test ="a[:2:3]"
        desired = {'begin':'_beg', 'end':'2', 'step':'3',
                   'single_index':'_index'}
        self.generic_check(test,desired)
    def test_begin_end_stp3_slice(self):
        """match slice from a[1:2:3]"""
        test ="a[1:2:3]"
        desired = {'begin':'1', 'end':'2', 'step':'3','single_index':'_index'}
        self.generic_check(test,desired)
    def test_expr_3_slice(self):
        """match slice from a[:1+i+2:]"""
        test ="a[:1+i+2:]"
        desired = {'begin':'_beg', 'end':"1+i+2",'step':'_stp',
                   'single_index':'_index'}
        self.generic_check(test,desired)
    def test_single_index(self):
        """match slice from a[0]"""
        test ="a[0]"
        desired = {'begin':'_beg', 'end':"_end",'step':'_stp',
                   'single_index':'0'}
        self.generic_check(test,desired)

def replace_whitespace(in_str):
    out = in_str.replace(" ","")
    out = out.replace("\t","")
    out = out.replace("\n","")
    return out

class TestTransformSlices(TestCase):
    def generic_check(self,suite_string,desired):
        import parser
        ast_list = parser.suite(suite_string).tolist()
        slice_handler.transform_slices(ast_list)
        actual = ast_to_string(ast_list)
        # Remove white space from expressions so that equivelant
        # but differently formatted string will compare equally
        actual = replace_whitespace(actual)
        desired = replace_whitespace(desired)
        assert_equal(actual,desired,suite_string)

    def test_simple_expr(self):
        """transform a[:] to slice notation"""
        test ="a[:]"
        desired = 'a[slice(_beg,_end,_stp)]'
        self.generic_check(test,desired)
    def test_simple_expr(self):
        """transform a[:,:] = b[:,1:1+2:3] *(c[1-2+i:,:] - c[:,:])"""
        test ="a[:,:] = b[:,1:1+2:3] *(c[1-2+i:,:] - c[:,:])"
        desired = " a[slice(_beg,_end),slice(_beg,_end)] = "\
                                    " b[slice(_beg,_end), slice(1,1+2-1,3)] *"\
                                    " (c[slice(1-2+i,_end), slice(_beg,_end)] -"\
                                    "  c[slice(_beg,_end), slice(_beg,_end)])"
        self.generic_check(test,desired)


if __name__ == "__main__":
    import nose
    nose.run(argv=['', __file__])
