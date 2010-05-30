import nose
from numpy.testing import *

from scipy.weave import ext_tools, c_spec
try:
    from scipy.weave.standard_array_spec import array_converter
except ImportError:
    pass # requires numpy.numerix


from weave_test_utils import *

build_dir = empty_temp_dir()

class TestExtModule(TestCase):
    #should really do some testing of where modules end up
    @dec.slow
    def test_simple(self):
        """ Simplest possible module """
        mod = ext_tools.ext_module('simple_ext_module')
        mod.compile(location = build_dir)
        import simple_ext_module
    @dec.slow
    def test_multi_functions(self):
        mod = ext_tools.ext_module('module_multi_function')
        var_specs = []
        code = ""
        test = ext_tools.ext_function_from_specs('test',code,var_specs)
        mod.add_function(test)
        test2 = ext_tools.ext_function_from_specs('test2',code,var_specs)
        mod.add_function(test2)
        mod.compile(location = build_dir)
        import module_multi_function
        module_multi_function.test()
        module_multi_function.test2()
    @dec.slow
    def test_with_include(self):
        # decalaring variables
        a = 2.;

        # declare module
        mod = ext_tools.ext_module('ext_module_with_include')
        mod.customize.add_header('<iostream>')

        # function 2 --> a little more complex expression
        var_specs = ext_tools.assign_variable_types(['a'],locals(),globals())
        code = """
               std::cout.clear(std::ios_base::badbit);
               std::cout << std::endl;
               std::cout << "test printing a value:" << a << std::endl;
               std::cout.clear(std::ios_base::goodbit);
               """
        test = ext_tools.ext_function_from_specs('test',code,var_specs)
        mod.add_function(test)
        # build module
        mod.compile(location = build_dir)
        import ext_module_with_include
        ext_module_with_include.test(a)

    @dec.slow
    def test_string_and_int(self):
        # decalaring variables
        a = 2;b = 'string'
        # declare module
        mod = ext_tools.ext_module('ext_string_and_int')
        code = """
               a=b.length();
               return_val = PyInt_FromLong(a);
               """
        test = ext_tools.ext_function('test',code,['a','b'])
        mod.add_function(test)
        mod.compile(location = build_dir)
        import ext_string_and_int
        c = ext_string_and_int.test(a,b)
        assert(c == len(b))

    @dec.slow
    def test_return_tuple(self):
        # decalaring variables
        a = 2
        # declare module
        mod = ext_tools.ext_module('ext_return_tuple')
        var_specs = ext_tools.assign_variable_types(['a'],locals())
        code = """
               int b;
               b = a + 1;
               py::tuple returned(2);
               returned[0] = a;
               returned[1] = b;
               return_val = returned;
               """
        test = ext_tools.ext_function('test',code,['a'])
        mod.add_function(test)
        mod.compile(location = build_dir)
        import ext_return_tuple
        c,d = ext_return_tuple.test(a)
        assert(c==a and d == a+1)

class TestExtFunction(TestCase):
    #should really do some testing of where modules end up
    @dec.slow
    def test_simple(self):
        """ Simplest possible function """
        mod = ext_tools.ext_module('simple_ext_function')
        var_specs = []
        code = ""
        test = ext_tools.ext_function_from_specs('test',code,var_specs)
        mod.add_function(test)
        mod.compile(location = build_dir)
        import simple_ext_function
        simple_ext_function.test()

class TestAssignVariableTypes(TestCase):
    def test_assign_variable_types(self):
        try:
            from numpy.numerix import arange, Float32, Float64
        except:
            # skip this test if numpy.numerix not installed
            return

        import types
        a = arange(10,typecode = Float32)
        b = arange(5,typecode = Float64)
        c = 5
        arg_list = ['a','b','c']
        actual = ext_tools.assign_variable_types(arg_list,locals())
        #desired = {'a':(Float32,1),'b':(Float32,1),'i':(Int32,0)}

        ad = array_converter()
        ad.name, ad.var_type, ad.dims = 'a', Float32, 1
        bd = array_converter()
        bd.name, bd.var_type, bd.dims = 'b', Float64, 1

        cd = c_spec.int_converter()
        cd.name, cd.var_type = 'c', types.IntType
        desired = [ad,bd,cd]
        assert_equal(actual,desired)

if __name__ == "__main__":
    import nose
    nose.run(argv=['', __file__])
