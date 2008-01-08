''' Define test function for scipy package '''
try:
    import nose
    from scipy.testing.nosetester import NoseTester as Tester
except ImportError:
    from scipy.testing.nulltester import NullTester as Tester

