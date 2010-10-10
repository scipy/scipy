import os
import sys

# Note: test_dir is global to this file.
#       It is made by setup_test_location()

#globals
global test_dir
test_dir = ''

import numpy
from numpy.testing import TestCase, dec, assert_

from scipy.weave import inline_tools,ext_tools
from scipy.weave.build_tools import msvc_exists, gcc_exists
from scipy.weave.catalog import unique_file
from scipy.weave.numpy_scalar_spec import numpy_complex_scalar_converter

def unique_mod(d,file_name):
    f = os.path.basename(unique_file(d,file_name))
    m = os.path.splitext(f)[0]
    return m

def remove_whitespace(in_str):
    out = in_str.replace(" ","")
    out = out.replace("\t","")
    out = out.replace("\n","")
    return out

#----------------------------------------------------------------------------
# Scalar conversion test classes
#   int, float, complex
#----------------------------------------------------------------------------

class NumpyComplexScalarConverter(TestCase):

    compiler = ''

    def setUp(self):
        self.converter = numpy_complex_scalar_converter()

    @dec.slow
    def test_type_match_string(self):
        assert_( not self.converter.type_match('string') )

    @dec.slow
    def test_type_match_int(self):
        assert_( not self.converter.type_match(5))

    @dec.slow
    def test_type_match_float(self):
        assert_( not self.converter.type_match(5.))

    @dec.slow
    def test_type_match_complex128(self):
        assert_(self.converter.type_match(numpy.complex128(5.+1j)))

    @dec.slow
    def test_complex_var_in(self):
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

    @dec.slow
    def test_complex_return(self):
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
        assert_( c == 3.+3j)

    @dec.slow
    def test_inline(self):
        a = numpy.complex128(1+1j)
        result = inline_tools.inline("return_val=1.0/a;",['a'])
        assert_( result==.5-.5j)

# class TestMsvcNumpyComplexScalarConverter(
#                   TestNumpyComplexScalarConverter):
#     compiler = 'msvc'
# class TestUnixNumpyComplexScalarConverter(
#                   TestNumpyComplexScalarConverter):
#     compiler = ''
# class TestGccNumpyComplexScalarConverter(
#                   TestNumpyComplexScalarConverter):
#     compiler = 'gcc'
for _n in dir():
    if _n[-9:]=='Converter':
        if msvc_exists():
            exec "class Test%sMsvc(%s):\n    compiler = 'msvc'"%(_n,_n)
        else:
            exec "class Test%sUnix(%s):\n    compiler = ''"%(_n,_n)
        if gcc_exists():
            exec "class Test%sGcc(%s):\n    compiler = 'gcc'"%(_n,_n)


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
        if _n[:8]=='TestMsvc': exec 'del '+_n
else:
    for _n in dir():
        if _n[:8]=='TestUnix': exec 'del '+_n

if not (gcc_exists() and msvc_exists() and sys.platform == 'win32'):
    for _n in dir():
        if _n[:7]=='TestGcc': exec 'del '+_n


if __name__ == "__main__":
    import nose
    nose.run(argv=['', __file__])
