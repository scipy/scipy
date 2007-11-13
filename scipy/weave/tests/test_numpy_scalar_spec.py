import time
import os,sys

# Note: test_dir is global to this file.
#       It is made by setup_test_location()

#globals
global test_dir
test_dir = ''

import numpy
from numpy.testing import *
set_package_path()
from weave import inline_tools,ext_tools
from weave.build_tools import msvc_exists, gcc_exists
from weave.catalog import unique_file
from weave.numpy_scalar_spec import numpy_complex_scalar_converter

restore_path()


def unique_mod(d,file_name):
    f = os.path.basename(unique_file(d,file_name))
    m = os.path.splitext(f)[0]
    return m

def remove_whitespace(in_str):
    import string
    out = string.replace(in_str," ","")
    out = string.replace(out,"\t","")
    out = string.replace(out,"\n","")
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

#----------------------------------------------------------------------------
# Scalar conversion test classes
#   int, float, complex
#----------------------------------------------------------------------------

class TestNumpyComplexScalarConverter(NumpyTestCase):
    compiler = ''

    def setUp(self):
        self.converter = numpy_complex_scalar_converter()

    def check_type_match_string(self,level=5):
        assert( not self.converter.type_match('string') )
    def check_type_match_int(self,level=5):
        assert( not self.converter.type_match(5))
    def check_type_match_float(self,level=5):
        assert( not self.converter.type_match(5.))
    def check_type_match_complex128(self,level=5):
        assert(self.converter.type_match(numpy.complex128(5.+1j)))

    def check_complex_var_in(self,level=5):
        mod_name = sys._getframe().f_code.co_name + self.compiler
        mod_name = unique_mod(test_dir,mod_name)
        mod = ext_tools.ext_module(mod_name)
        a = numpy.complex(1.+1j)
        code = "a=std::complex<double>(2.,2.);"
        test = ext_tools.ext_function('test',code,['a'])
        mod.add_function(test)
        mod.compile(location = test_dir, compiler = self.compiler)
        exec 'from ' + mod_name + ' import test'
        b=numpy.complex128(1.+1j)
        test(b)
        try:
            b = 1.
            test(b)
        except TypeError:
            pass
        try:
            b = 'abc'
            test(b)
        except TypeError:
            pass

    def check_complex_return(self,level=5):
        mod_name = sys._getframe().f_code.co_name + self.compiler
        mod_name = unique_mod(test_dir,mod_name)
        mod = ext_tools.ext_module(mod_name)
        a = 1.+1j
        code = """
               a= a + std::complex<double>(2.,2.);
               return_val = PyComplex_FromDoubles(a.real(),a.imag());
               """
        test = ext_tools.ext_function('test',code,['a'])
        mod.add_function(test)
        mod.compile(location = test_dir, compiler = self.compiler)
        exec 'from ' + mod_name + ' import test'
        b=1.+1j
        c = test(b)
        assert( c == 3.+3j)

    def check_inline(self, level=5):
        a = numpy.complex128(1+1j)
        result = inline_tools.inline("return_val=1.0/a;",['a'])
        assert( result==.5-.5j)

class TestMsvcNumpyComplexScalarConverter(
                  TestNumpyComplexScalarConverter):
    compiler = 'msvc'
class TestUnixNumpyComplexScalarConverter(
                  TestNumpyComplexScalarConverter):
    compiler = ''
class TestGccNumpyComplexScalarConverter(
                  TestNumpyComplexScalarConverter):
    compiler = 'gcc'


def setup_test_location():
    import tempfile
    #test_dir = os.path.join(tempfile.gettempdir(),'test_files')
    test_dir = tempfile.mktemp()
    if not os.path.exists(test_dir):
        os.mkdir(test_dir)
    sys.path.insert(0,test_dir)
    return test_dir

test_dir = setup_test_location()

def teardown_test_location():
    import tempfile
    test_dir = os.path.join(tempfile.gettempdir(),'test_files')
    if sys.path[0] == test_dir:
        sys.path = sys.path[1:]
    return test_dir

def remove_file(name):
    test_dir = os.path.abspath(name)

if not msvc_exists():
    for _n in dir():
        if _n[:10]=='test_msvc_': exec 'del '+_n
else:
    for _n in dir():
        if _n[:10]=='test_unix_': exec 'del '+_n

if not (gcc_exists() and msvc_exists() and sys.platform == 'win32'):
    for _n in dir():
        if _n[:9]=='test_gcc_': exec 'del '+_n

if __name__ == "__main__":
    NumpyTest('weave.numpy_scalar_spec').run()
