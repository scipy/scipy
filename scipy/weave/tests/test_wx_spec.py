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
set_package_path()
from weave import ext_tools, wx_spec
restore_path()

import wx

class TestWxConverter(NumpyTestCase):
    def setUp(self):
        self.app = wx.App()
        self.s = wx_spec.wx_converter()

    def check_type_match_string(self,level=5):
        assert(not self.s.type_match('string') )

    def check_type_match_int(self,level=5):
        assert(not self.s.type_match(5))

    def check_type_match_float(self,level=5):
        assert(not self.s.type_match(5.))

    def check_type_match_complex(self,level=5):
        assert(not self.s.type_match(5.+1j))

    def check_type_match_complex(self,level=5):
        assert(not self.s.type_match(5.+1j))

    def check_type_match_wxframe(self,level=5):
        f=wx.Frame(None,-1,'bob')
        assert(self.s.type_match(f))

    def check_var_in(self,level=5):
        mod = ext_tools.ext_module('wx_var_in',compiler='msvc')
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
        except TypeError:
            pass
        try:
            b = 1
            wx_var_in.test(b)
        except TypeError:
            pass

    def no_check_var_local(self,level=5):
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

    def no_test_no_check_return(self,level=5):
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

if __name__ == "__main__":
    import sys
    if len(sys.argv) == 1:
        sys.argv.extend(["--level=5"])
    NumpyTest().run()
