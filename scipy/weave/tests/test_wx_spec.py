"""
check_var_in -- tests whether a variable is passed in correctly
                and also if the passed in variable can be reassigned
check_var_local -- tests wheter a variable is passed in , modified,
                   and returned correctly in the local_dict dictionary
                   argument
check_return -- test whether a variable is passed in, modified, and
                then returned as a function return value correctly
"""

from numpy.testing import *

e = None
DONOTRUN = False
try:
    from scipy.weave import ext_tools, wx_spec
    import wx
except ImportError, e:
    wx = None
    DONOTRUN = True
except RuntimeError, e:
    wx = None
    DONOTRUN = True

skip = dec.skipif(DONOTRUN, "(error was %s)" % str(e))

class TestWxConverter(TestCase):
    def setUp(self):
        if not DONOTRUN:
            self.app = wx.App()
            self.s = wx_spec.wx_converter()

    @dec.slow
    def test_type_match_string(self):
        assert(not self.s.type_match('string') )

    @dec.slow
    def test_type_match_int(self):
        assert(not self.s.type_match(5))

    @dec.slow
    def test_type_match_float(self):
        assert(not self.s.type_match(5.))

    @dec.slow
    def test_type_match_complex(self):
        assert(not self.s.type_match(5.+1j))

    @dec.slow
    def test_type_match_complex(self):
        assert(not self.s.type_match(5.+1j))

    @dec.slow
    def test_type_match_wxframe(self):
        f=wx.Frame(None,-1,'bob')
        assert(self.s.type_match(f))

    @dec.slow
    def test_var_in(self):
        mod = ext_tools.ext_module('wx_var_in',compiler='')
        mod.customize.add_header('<wx/string.h>')
        mod.customize.add_extra_compile_arg(' '.join(self.s.extra_compile_args))
        mod.customize.add_extra_link_arg(' '.join(self.s.extra_link_args))

        a = wx.Frame(None,-1,'bob')
        code = """
               py::tuple args(1);
               args[0] = py::object("jim");
               a.mcall("SetTitle",args);
               """
        test = ext_tools.ext_function('test',code,['a'],locals(),globals())
        mod.add_function(test)
        mod.compile()
        import wx_var_in
        b=a
        wx_var_in.test(b)
        assert(b.GetTitle() == "jim")
        try:
            b = 1.
            wx_var_in.test(b)
        except AttributeError:
            pass
        try:
            b = 1
            wx_var_in.test(b)
        except AttributeError:
            pass

    @dec.slow
    def no_check_var_local(self):
        mod = ext_tools.ext_module('wx_var_local')
        a = 'string'
        code = 'a="hello";'
        var_specs = ext_tools.assign_variable_types(['a'],locals())
        test = ext_tools.ext_function_from_specs('test',code,var_specs)
        mod.add_function(test)
        mod.compile()
        import wx_var_local
        b='bub'
        q={}
        wx_var_local.test(b,q)
        assert('a' == 'string')

    @dec.slow
    def no_test_no_check_return(self):
        mod = ext_tools.ext_module('wx_return')
        a = 'string'
        code = """
               a= Py::wx("hello");
               return_val = Py::new_reference_to(a);
               """
        test = ext_tools.ext_function('test',code,['a'],locals())
        mod.add_function(test)
        mod.compile()
        import wx_return
        b='bub'
        c = wx_return.test(b)
        assert(c == 'hello')

decorate_methods(TestWxConverter, skip)

if __name__ == "__main__":
    nose.run(argv=['', __file__])
