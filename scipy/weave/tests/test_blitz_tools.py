from __future__ import absolute_import, print_function

import time
import parser
import warnings

from numpy import (float32, float64, complex64, complex128,
                   zeros, random, array)

from numpy.testing import (TestCase, assert_equal,
                           assert_allclose, run_module_suite)

from scipy.weave import blitz_tools, blitz, BlitzWarning
from scipy.weave.ast_tools import harvest_variables
from weave_test_utils import remove_whitespace, debug_print, TempdirBlitz, dec


class TestAstToBlitzExpr(TestCase):
    def generic_check(self,expr,desired):
        ast = parser.suite(expr)
        ast_list = ast.tolist()
        actual = blitz_tools.ast_to_blitz_expr(ast_list)
        actual = remove_whitespace(actual)
        desired = remove_whitespace(desired)
        assert_equal(actual,desired,expr)

    def test_simple_expr(self):
        # convert simple expr to blitz
        expr = "a[:1:2] = b[:1+i+2:]"
        desired = "a(blitz::Range(_beg,1-1,2))="\
                  "b(blitz::Range(_beg,1+i+2-1));"
        self.generic_check(expr,desired)

    def test_fdtd_expr(self):
        # Convert fdtd equation to blitz.
        # Note:  This really should have "\" at the end of each line to
        # indicate continuation.
        expr = "ex[:,1:,1:] =   ca_x[:,1:,1:] * ex[:,1:,1:]" \
                             "+ cb_y_x[:,1:,1:] * (hz[:,1:,1:] - hz[:,:-1,:])"\
                             "- cb_z_x[:,1:,1:] * (hy[:,1:,1:] - hy[:,1:,:-1])"
        desired = 'ex(_all,blitz::Range(1,_end),blitz::Range(1,_end))='\
                  '  ca_x(_all,blitz::Range(1,_end),blitz::Range(1,_end))'\
                  ' *ex(_all,blitz::Range(1,_end),blitz::Range(1,_end))'\
                  '+cb_y_x(_all,blitz::Range(1,_end),blitz::Range(1,_end))'\
                  '*(hz(_all,blitz::Range(1,_end),blitz::Range(1,_end))'\
                  '  -hz(_all,blitz::Range(_beg,Nhz(1)-1-1),_all))'\
                  ' -cb_z_x(_all,blitz::Range(1,_end),blitz::Range(1,_end))'\
                  '*(hy(_all,blitz::Range(1,_end),blitz::Range(1,_end))'\
                  '-hy(_all,blitz::Range(1,_end),blitz::Range(_beg,Nhy(2)-1-1)));'
        self.generic_check(expr,desired)


class TestBlitz(TestCase):
    """These are long running tests...

    Would be useful to benchmark these things somehow.
    """
    def generic_check(self, expr, arg_dict, type, size):
        clean_result = array(arg_dict['result'],copy=1)
        t1 = time.time()
        exec(expr, globals(),arg_dict)
        t2 = time.time()
        standard = t2 - t1
        desired = arg_dict['result']
        arg_dict['result'] = clean_result
        t1 = time.time()
        blitz_tools.blitz(expr,arg_dict,{},verbose=0)
        t2 = time.time()
        compiled = t2 - t1
        actual = arg_dict['result']
        # TODO: this isn't very stringent.  Need to tighten this up and
        # learn where failures are occurring.
        assert_allclose(abs(actual.ravel()), abs(desired.ravel()),
                        rtol=1e-4, atol=1e-6)
        return standard, compiled

    def generic_2d(self,expr,typ):
        # The complex testing is pretty lame...
        ast = parser.suite(expr)
        arg_list = harvest_variables(ast.tolist())
        all_sizes = [(10,10), (50,50), (100,100), (500,500), (1000,1000)]
        debug_print('\nExpression:', expr)
        with TempdirBlitz():
            for size in all_sizes:
                arg_dict = {}
                for arg in arg_list:
                    arg_dict[arg] = random.normal(0,1,size).astype(typ)
                    # set imag part of complex values to non-zero value
                    try:
                        arg_dict[arg].imag = arg_dict[arg].real
                    except:
                        pass
                debug_print('Run:', size,typ)
                standard,compiled = self.generic_check(expr,arg_dict,type,size)
                try:
                    speed_up = standard/compiled
                except:
                    speed_up = -1.
                debug_print("1st run(numpy,compiled,speed up):  %3.4f, %3.4f, "
                            "%3.4f" % (standard,compiled,speed_up))
                standard,compiled = self.generic_check(expr,arg_dict,type,size)
                try:
                    speed_up = standard/compiled
                except:
                    speed_up = -1.
                debug_print("2nd run(numpy,compiled,speed up):  %3.4f, %3.4f, "
                            "%3.4f" % (standard,compiled,speed_up))

    @dec.slow
    def test_5point_avg_2d_float(self):
        expr = "result[1:-1,1:-1] = (b[1:-1,1:-1] + b[2:,1:-1] + b[:-2,1:-1]" \
                                  "+ b[1:-1,2:] + b[1:-1,:-2]) / 5."
        self.generic_2d(expr,float32)

    @dec.slow
    def test_5point_avg_2d_double(self):
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=BlitzWarning)
            expr = "result[1:-1,1:-1] = (b[1:-1,1:-1] + b[2:,1:-1] + b[:-2,1:-1]" \
                                      "+ b[1:-1,2:] + b[1:-1,:-2]) / 5."
            self.generic_2d(expr,float64)

    @dec.slow
    def _check_5point_avg_2d_complex_float(self):
        """ Note: THIS TEST is KNOWN TO FAIL ON GCC 3.x.
            It will not adversely affect 99.99 percent of weave

        result[1:-1,1:-1] = (b[1:-1,1:-1] + b[2:,1:-1] + b[:-2,1:-1]
                             + b[1:-1,2:] + b[1:-1,:-2]) / 5.

        Note: THIS TEST is KNOWN TO FAIL ON GCC 3.x.  The reason is that
        5. is a double and b is a complex32.  blitz doesn't know
        how to handle complex32/double.  See:
        http://www.oonumerics.org/MailArchives/blitz-support/msg00541.php
        Unfortunately, the fix isn't trivial.  Instead of fixing it, I
        prefer to wait until we replace blitz++ with Pat Miller's code
        that doesn't rely on blitz..
        """
        expr = "result[1:-1,1:-1] = (b[1:-1,1:-1] + b[2:,1:-1] + b[:-2,1:-1]" \
                                  "+ b[1:-1,2:] + b[1:-1,:-2]) / 5."
        self.generic_2d(expr,complex64)

    @dec.slow
    def test_5point_avg_2d_complex_double(self):
        expr = "result[1:-1,1:-1] = (b[1:-1,1:-1] + b[2:,1:-1] + b[:-2,1:-1]" \
                                  "+ b[1:-1,2:] + b[1:-1,:-2]) / 5."
        self.generic_2d(expr,complex128)


@dec.slow
def test_blitz_bug():
    # Assignment to arr[i:] used to fail inside blitz expressions.
    with TempdirBlitz():
        N = 4
        expr_buggy = 'arr_blitz_buggy[{0}:] = arr[{0}:]'
        expr_not_buggy = 'arr_blitz_not_buggy[{0}:{1}] = arr[{0}:]'
        random.seed(7)
        arr = random.randn(N)
        sh = arr.shape[0]
        for lim in [0, 1, 2]:
            arr_blitz_buggy = zeros(N)
            arr_blitz_not_buggy = zeros(N)
            arr_np = zeros(N)
            blitz(expr_buggy.format(lim))
            blitz(expr_not_buggy.format(lim, 'sh'))
            arr_np[lim:] = arr[lim:]
            assert_allclose(arr_blitz_buggy, arr_np)
            assert_allclose(arr_blitz_not_buggy, arr_np)


if __name__ == "__main__":
    run_module_suite()
