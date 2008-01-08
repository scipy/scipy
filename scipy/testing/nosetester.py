''' Nose tester object '''
import os
import sys

import nose

class NoseTester(object):
    """ Scipy nose tests suite manager.

    Usage: NoseTester(<package>).test()

    <package> is package path, package name or its module object.
    """
    def __init__(self, package=None):
        if package is None:
            f = sys._getframe(1)
            package = f.f_locals.get('__file__', None)
            assert package is not None
            self.package_path = os.path.dirname(package)
        else:
            self.package_path = self._process_package(package)

    def _process_package(self, package):
        ''' Package can be module, path, or package name '''
        try:
            pfile = package.__file__
        except AttributeError:
            pass
        else:
            return os.path.dirname(pfile)
        if os.path.isabs(package):
            return package
        if package.find(os.path.sep) == -1:
            # Try scipy package name
            import scipy
            scipy.pkgload(package)
            try:
                module = getattr(scipy, package)
            except AttributeError:
                pass
            else:
                self.package_path = os.path.dirname(module.__file__)
                return
        # Default to relative path
        self.package_path = os.path.abspath(package)
        return

    def test(self, labels='fast', verbose=1, doctests=False, extra_argv=None):
        ''' Module testing function

        labels - identifies tests to run.  This can be a string to
          pass to the nostests executable with the '-a'
          option, or one of several special values.
          Special values are:
          'fast' - the default - which corresponds to
             nosetests -a option of 'not slow and not bench'.
          None or '' - run all tests and benchmarks

        verbose - verbosity value 1-10
        doctests - if True, run doctests in module
        extra_argv - list with any extra args to pass to nosetest
        '''
        argv = ['scipy module test', self.package_path, '-s']
        if labels:
            if labels == 'fast':
                labels = 'not slow and not bench'
            argv += ['-A', labels]
        argv += ['--verbosity', str(verbose)]
        if doctests:
            argv+=['--with-doctest']
        if extra_argv:
            argv+= extra_argv
        nose.run(argv=argv)
        
