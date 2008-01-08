''' Nose tester object '''
import os
import sys

import nose

class NoseTester(object):
    """ Scipy nose tests site manager.

    Usage: NoseTester(<package>).test()

    <package> is package path, package name or its module object.
    """
    def __init__(self, package=None):
        if package is None:
            f = sys._getframe(1)
            package = f.f_locals.get('__file__',f.f_globals.get('__file__',None))
            assert package is not None
            self.package_path = os.path.dirname(package)
        else:
            self._process_package(package)

    def _process_package(self, package):
        ''' Package can be module, path, or package name '''
        try:
            pfile = package.__file__
        except AttributeError:
            pass
        else:
            self.package_path = os.path.dirname(pfile)
            return
        if os.path.isabs(package):
            self.package_path = package
            return
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

    def test(self, labels=None, *args, **kwargs):
        if labels is None:
            labels = []
        argv = ['', self.package_path]
        nose.run(argv=argv)
        
