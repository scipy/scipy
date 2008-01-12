''' Null tester (when nose not importable '''

class NullTester(object):
    def __init__(self, *args, **kwargs):
        pass
    def test(self, labels=None, *args, **kwargs):
        raise ImportError, 'Need nose testing on path for tests'
    
