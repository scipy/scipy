''' Define test function for scipy package '''
try:
    import nose
except ImportError:
    from scipy.testing.nulltester import NullTester as Tester
else:
    from scipy.testing.nosetester import NoseTester as Tester

