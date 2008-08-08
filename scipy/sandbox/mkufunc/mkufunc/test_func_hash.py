import unittest

from mkufunc.api import func_hash


class Tests(unittest.TestCase):

    # These tests are very (Python) version specific.

    def test_simple(self):

        def f(x):
            return 2.5 * x * x + 4.7 * x

        self.assertEqual(func_hash(f),
                         '5f12e97debf1d2cb9e0a2f92e045b1fb')


    def test_extra(self):

        def f(x):
            return 2.5 * x * x + 4.7 * x

        self.assertEqual(func_hash(f, salt=[(int, int), (float, float)]),
                         'e637d9825ef20cb56d364041118ca72e')

    def test_const(self):

        def add_a(b):
            return a + b   # a in globals

        self.assertEqual(func_hash(add_a),
                         '9ff237f372bf233470ce940edd58f60d')

    def test_inner(self):

        def foo(x):
            inner1 = lambda t: t/3.0
            def inner2(n):
                return n + 3
            return inner1(x) + inner2(int(x))

        self.assertEqual(func_hash(foo),
                         '814c113dfc77e7ebb52915dd3ce9c37a')


if __name__ == '__main__':
    unittest.main()
