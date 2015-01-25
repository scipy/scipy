from __future__ import absolute_import, print_function

import random
import parser

import numpy as np
from numpy.testing import TestCase, assert_array_equal, run_module_suite

from scipy.weave import size_check
from scipy.weave.ast_tools import harvest_variables


class TestMakeSameLength(TestCase):
    def generic_check(self,x,y,desired):
        actual = size_check.make_same_length(x,y)
        desired = desired
        assert_array_equal(actual,desired)

    def test_scalar(self):
        x,y = (),()
        desired = np.array(()), np.array(())
        self.generic_check(x,y,desired)

    def test_x_scalar(self):
        x,y = (),(1,2)
        desired = np.array((1,1)), np.array((1,2))
        self.generic_check(x,y,desired)

    def test_y_scalar(self):
        x,y = (1,2),()
        desired = np.array((1,2)), np.array((1,1))
        self.generic_check(x,y,desired)

    def test_x_short(self):
        x,y = (1,2),(1,2,3)
        desired = np.array((1,1,2)), np.array((1,2,3))
        self.generic_check(x,y,desired)

    def test_y_short(self):
        x,y = (1,2,3),(1,2)
        desired = np.array((1,2,3)), np.array((1,1,2))
        self.generic_check(x,y,desired)


class TestBinaryOpSize(TestCase):
    def generic_check(self,x,y,desired):
        actual = size_check.binary_op_size(x,y)
        desired = desired
        assert_array_equal(actual,desired)

    def generic_error_check(self,x,y):
        self.assertRaises(ValueError, size_check.binary_op_size, x, y)

    def desired_type(self,val):
        return np.array(val)

    def test_scalar(self):
        x,y = (),()
        desired = self.desired_type(())
        self.generic_check(x,y,desired)

    def test_x1(self):
        x,y = (1,),()
        desired = self.desired_type((1,))
        self.generic_check(x,y,desired)

    def test_y1(self):
        x,y = (),(1,)
        desired = self.desired_type((1,))
        self.generic_check(x,y,desired)

    def test_x_y(self):
        x,y = (5,),(5,)
        desired = self.desired_type((5,))
        self.generic_check(x,y,desired)

    def test_x_y2(self):
        x,y = (5,10),(5,10)
        desired = self.desired_type((5,10))
        self.generic_check(x,y,desired)

    def test_x_y3(self):
        x,y = (5,10),(1,10)
        desired = self.desired_type((5,10))
        self.generic_check(x,y,desired)

    def test_x_y4(self):
        x,y = (1,10),(5,10)
        desired = self.desired_type((5,10))
        self.generic_check(x,y,desired)

    def test_x_y5(self):
        x,y = (5,1),(1,10)
        desired = self.desired_type((5,10))
        self.generic_check(x,y,desired)

    def test_x_y6(self):
        x,y = (1,10),(5,1)
        desired = self.desired_type((5,10))
        self.generic_check(x,y,desired)

    def test_x_y7(self):
        x,y = (5,4,3,2,1),(3,2,1)
        desired = self.desired_type((5,4,3,2,1))
        self.generic_check(x,y,desired)

    def test_error1(self):
        x,y = (5,),(4,)
        self.generic_error_check(x,y)

    def test_error2(self):
        x,y = (5,5),(4,5)
        self.generic_error_check(x,y)


class TestDummyArray(TestBinaryOpSize):
    def generic_check(self,x,y,desired):
        if isinstance(x, tuple):
            x = np.ones(x)
        if isinstance(y, tuple):
            y = np.ones(y)
        xx = size_check.dummy_array(x)
        yy = size_check.dummy_array(y)
        ops = ['+', '-', '/', '*', '<<', '>>']
        for op in ops:
            actual = eval('xx' + op + 'yy')
            desired = desired
            assert_array_equal(actual,desired)

    def desired_type(self,val):
        return size_check.dummy_array(np.array(val),1)


class TestDummyArrayIndexing(TestCase):
    def generic_check(self,ary,expr,desired):
        a = size_check.dummy_array(ary)
        actual = eval(expr).shape
        assert_array_equal(actual,desired, expr)

    def generic_wrap(self,a,expr):
        desired = np.array(eval(expr).shape)
        try:
            self.generic_check(a,expr,desired)
        except IndexError:
            if 0 not in desired:
                msg = '%s raised IndexError in dummy_array, but forms\n' \
                      'valid array shape -> %s' % (expr, str(desired))
                raise AttributeError(msg)

    def generic_1d(self,expr):
        a = np.arange(10)
        self.generic_wrap(a,expr)

    def generic_2d(self,expr):
        a = np.ones((10,20))
        self.generic_wrap(a,expr)

    def generic_3d(self,expr):
        a = np.ones((10,20,1))
        self.generic_wrap(a,expr)

    def generic_1d_index(self,expr):
        a = np.arange(10)
        desired = np.array(())
        self.generic_check(a,expr,desired)

    def test_1d_index_0(self):
        self.generic_1d_index('a[0]')

    def test_1d_index_1(self):
        self.generic_1d_index('a[4]')

    def test_1d_index_2(self):
        self.generic_1d_index('a[-4]')

    def test_1d_index_3(self):
        try:
            self.generic_1d('a[12]')
        except IndexError:
            pass

    def test_1d_index_calculated(self):
        self.generic_1d_index('a[0+1]')

    def test_1d_0(self):
        self.generic_1d('a[:]')

    def test_1d_1(self):
        self.generic_1d('a[1:]')

    def test_1d_2(self):
        self.generic_1d('a[-1:]')

    def test_1d_3(self):
        self.generic_1d('a[-11:]')

    def test_1d_4(self):
        self.generic_1d('a[:1]')

    def test_1d_5(self):
        self.generic_1d('a[:-1]')

    def test_1d_6(self):
        self.generic_1d('a[:-11]')

    def test_1d_7(self):
        self.generic_1d('a[1:5]')

    def test_1d_8(self):
        self.generic_1d('a[1:-5]')

    def test_1d_9(self):
        # don't support zero length slicing at the moment.
        try:
            self.generic_1d('a[-1:-5]')
        except IndexError:
            pass

    def test_1d_10(self):
        self.generic_1d('a[-5:-1]')

    def test_1d_stride_0(self):
        self.generic_1d('a[::1]')

    def test_1d_stride_1(self):
        self.generic_1d('a[::-1]')

    def test_1d_stride_2(self):
        self.generic_1d('a[1::1]')

    def test_1d_stride_3(self):
        self.generic_1d('a[1::-1]')

    def test_1d_stride_4(self):
        # don't support zero length slicing at the moment.
        try:
            self.generic_1d('a[1:5:-1]')
        except IndexError:
            pass

    def test_1d_stride_5(self):
        self.generic_1d('a[5:1:-1]')

    def test_1d_stride_6(self):
        self.generic_1d('a[:4:1]')

    def test_1d_stride_7(self):
        self.generic_1d('a[:4:-1]')

    def test_1d_stride_8(self):
        self.generic_1d('a[:-4:1]')

    def test_1d_stride_9(self):
        self.generic_1d('a[:-4:-1]')

    def test_1d_stride_10(self):
        self.generic_1d('a[:-3:2]')

    def test_1d_stride_11(self):
        self.generic_1d('a[:-3:-2]')

    def test_1d_stride_12(self):
        self.generic_1d('a[:-3:-7]')

    def test_1d_random(self):
        # throw a bunch of different indexes at it for good measure.
        choices = map(lambda x: repr(x),range(50)) + range(50) + ['']*50
        for i in range(100):
            try:
                beg = random.choice(choices)
                end = random.choice(choices)
                step = random.choice(choices)
                if step in ['0',0]:
                    step = 'None'
                self.generic_1d('a[%s:%s:%s]' % (beg,end,step))
            except IndexError:
                pass

    def test_2d_0(self):
        self.generic_2d('a[:]')

    def test_2d_1(self):
        self.generic_2d('a[:2]')

    def test_2d_2(self):
        self.generic_2d('a[:,:]')

    def test_2d_random(self):
        # throw a bunch of different indexes at it for good measure.
        choices = map(lambda x: repr(x),range(50)) + range(50) + ['']*50
        for i in range(100):
            try:
                beg = random.choice(choices)
                end = random.choice(choices)
                step = random.choice(choices)
                beg2 = random.choice(choices)
                end2 = random.choice(choices)
                step2 = random.choice(choices)
                if step in ['0',0]:
                    step = 'None'
                if step2 in ['0',0]:
                    step2 = 'None'
                expr = 'a[%s:%s:%s,%s:%s:%s]' % (beg,end,step,beg2,end2,step2)
                self.generic_2d(expr)
            except IndexError:
                pass

    def test_3d_random(self):
        # throw a bunch of different indexes at it for good measure.
        choices = map(lambda x: repr(x),range(50)) + range(50) + ['']*50
        for i in range(100):
            try:
                idx = []
                for i in range(9):
                    val = random.choice(choices)
                    if (i+1) % 3 == 0 and val in ['0',0]:
                        val = 'None'
                    idx.append(val)
                expr = 'a[%s:%s:%s,%s:%s:%s,%s:%s:%s]' % tuple(idx)
                self.generic_3d(expr)
            except IndexError:
                pass


class TestReduction(TestCase):
    def test_1d_0(self):
        a = np.ones((5,))
        actual = size_check.reduction(a,0)
        desired = size_check.dummy_array((),1)
        assert_array_equal(actual.shape,desired.shape)

    def test_2d_0(self):
        a = np.ones((5,10))
        actual = size_check.reduction(a,0)
        desired = size_check.dummy_array((10,),1)
        assert_array_equal(actual.shape,desired.shape)

    def test_2d_1(self):
        a = np.ones((5,10))
        actual = size_check.reduction(a,1)
        desired = size_check.dummy_array((5,),1)
        assert_array_equal(actual.shape,desired.shape)

    def test_3d_0(self):
        a = np.ones((5,6,7))
        actual = size_check.reduction(a,1)
        desired = size_check.dummy_array((5,7),1)
        assert_array_equal(actual.shape,desired.shape)

    def test_error0(self):
        a = np.ones((5,))
        try:
            size_check.reduction(a,-2)
        except ValueError:
            pass

    def test_error1(self):
        a = np.ones((5,))
        try:
            size_check.reduction(a,1)
        except ValueError:
            pass


class TestExpressions(TestCase):
    def generic_check(self,expr,desired,**kw):
        ast_list = parser.expr(expr).tolist()
        args = harvest_variables(ast_list)
        loc = locals().update(kw)
        for var in args:
            s = '%s = size_check.dummy_array(%s)' % (var,var)
            exec(s,loc)
        try:
            actual = eval(expr,locals()).shape
        except:
            actual = 'failed'

        if actual is 'failed' and desired is 'failed':
            return

        assert_array_equal(actual,desired, expr)

    def generic_wrap(self,expr,**kw):
        try:
            x = np.array(eval(expr,kw))
            try:
                desired = x.shape
            except:
                desired = np.zeros(())
        except:
            desired = 'failed'
        self.generic_check(expr,desired,**kw)

    def test_generic_1d(self):
        a = np.arange(10)
        expr = 'a[:]'
        self.generic_wrap(expr,a=a)
        expr = 'a[:] + a'
        self.generic_wrap(expr,a=a)
        bad_expr = 'a[4:] + a'
        self.generic_wrap(bad_expr,a=a)
        a = np.arange(10)
        b = np.ones((1,10))
        expr = 'a + b'
        self.generic_wrap(expr,a=a,b=b)
        bad_expr = 'a[:5] + b'
        self.generic_wrap(bad_expr,a=a,b=b)

    def test_single_index(self):
        a = np.arange(10)
        expr = 'a[5] + a[3]'
        self.generic_wrap(expr,a=a)

    def test_calculated_index(self):
        a = np.arange(10)
        nx = 0
        expr = 'a[5] + a[nx+3]'
        size_check.check_expr(expr,locals())

    def test_calculated_index2(self):
        a = np.arange(10)
        nx = 0
        expr = 'a[1:5] + a[nx+1:5+nx]'
        size_check.check_expr(expr,locals())

    def generic_2d(self,expr):
        a = np.ones((10,20))
        self.generic_wrap(a,expr)

    def generic_3d(self,expr):
        a = np.ones((10,20,1))
        self.generic_wrap(a,expr)


if __name__ == "__main__":
    run_module_suite()
