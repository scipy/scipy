import sys
from _image import *
from combine import *

__version__ = '2.0'
if sys.version_info < (2,4):
    def test():
        import doctest, _image, combine

        t = doctest.Tester(globs = globals())

        t.rundict(_image.__dict__, "_image")
        t.rundict(combine.__dict__, "combine")

        return t.summarize()

else:
    def test():
        import doctest, _image, combine

        finder=doctest.DocTestFinder()
        tests=finder.find(_image)
        tests.extend(finder.find(combine))

        runner=doctest.DocTestRunner(verbose=False)

        for test in tests:
            runner.run(test)
        return runner.summarize()
