''' Nose tester object '''
import os
import sys

import nose

class NoseTester(object):
    """ Scipy nose test runner.

    Usage: NoseTester(<package>).test()

    <package> is package path or module - None finds calling module path
    """
    def __init__(self, package=None):
        if package is None:
            f = sys._getframe(1)
            package = f.f_locals.get('__file__', None)
            assert package is not None
            package = os.path.dirname(package)
        elif isinstance(package, type(os)):
            package = os.path.dirname(package.__file__)
        self.package_path = package
        
    def test(self, labels='fast', verbose=1, doctests=False, extra_argv=None):
        ''' Module testing function

        labels - identifies tests to run.  This can be a string to
          pass to the nosetests executable with the '-A'
          option, or one of several special values.
          Special values are:
          'fast' - the default - which corresponds to
             nosetests -A option of
             'not slow and not bench and not willfail'.
          'full' - fast (as above) and slow tests as in
             nosetests -A option of 'not bench and not willfail'.             
          None or '' - run all tests and benchmarks

        verbose - verbosity value 1-10
        doctests - if True, run doctests in module
        extra_argv - list with any extra args to pass to nosetest
        '''
        argv = ['scipy module test', self.package_path, '-s']
        if labels:
            if labels == 'fast':
                labels = 'not slow and not bench and not willfail'
            elif labels == 'full':
                labels = 'not bench and not willfail'
            argv += ['-A', labels]
        argv += ['--verbosity', str(verbose)]
        if doctests:
            argv+=['--with-doctest']
        if extra_argv:
            argv+= extra_argv
        nose.run(argv=argv)

        
