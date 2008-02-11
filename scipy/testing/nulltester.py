''' Null tester to signal nose tests disabled

Merely returns error reporting lack of nose package or version number
below requirements.

See pkgtester, nosetester modules

'''

class NullTester(object):
    _msg = 'Need nose >=0.10 for tests - see %s' % \
           'http://somethingaboutorange.com/mrl/projects/nose'
    def __init__(self, *args, **kwargs):
        pass
    def test(self, labels=None, *args, **kwargs):
        raise ImportError, self._msg
    def bench(self, labels=None, *args, **kwargs):
        raise ImportError, self._msg
    
