import unittest

from funcutil import func_hash


class Tests(unittest.TestCase):
    
    def test_simple(self):
        
        def f(x):
            return 2.5 * x * x + 4.7 * x
        
        self.assertEqual(func_hash(f),
                         'f8c94c2e2dd69226706f90c2f4293497')
        
        
    def test_extra(self):
        
        def f(x):
            return 2.5 * x * x + 4.7 * x
        
        self.assertEqual(func_hash(f, salt=[(int, int), (float, float)]),
                         'd81db2e37ade51a430e47b72c55e197e')
        
    def test_const(self):
        
        def add_a(b):
            return a + b   # a in globals
        
        self.assertEqual(func_hash(add_a),
                         '55a68633f905a1373f61659b41402f02')
        
    def test_inner(self):

        def foo(x):
            inner1 = lambda t: t/3.0
            def inner2(n):
                return n + 3
            return inner1(x) + inner2(int(x))

        #func_hash(foo, verbose=1)
        self.assertEqual(func_hash(foo),
                         'a836c2dbe1b202bd68e1fe3affe1ce7a')


if __name__ == '__main__':
    unittest.main()
