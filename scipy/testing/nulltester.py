''' Null tester (when nose not importable)

Merely returns error reporting lack of nose package

See pkgtester, nosetester modules

'''

nose_url = 'http://somethingaboutorange.com/mrl/projects/nose'

class NullTester(object):
    def __init__(self, *args, **kwargs):
        pass
    def test(self, labels=None, *args, **kwargs):
        raise ImportError, 'Need nose for tests - see %s' % nose_url
    def bench(self, labels=None, *args, **kwargs):
        raise ImportError, 'Need nose for benchmarks - see %s' % nose_url
    
