''' Nose test running

Implements .test functions for modules.

'''
import os
import sys

import nose

class NoseTester(object):
    """ Nose test runner.

    Usage: NoseTester(<package>).test()
    
    <package> is package path or module Default for package is None. A
    value of None finds calling module path.

    Typical call is from module __init__, and corresponds to this:
    
    >>> test = NoseTester().test
    
    In practice, because nose may not be importable, the __init__
    files actually have:
    
    >>> from scipy.testing.pkgtester import Tester
    >>> test = Tester().test
    
    The pkgtester module checks for the presence of nose on the path,
    returning this class if nose is present, and a null class
    otherwise.
    """
    def __init__(self, package=None):
        ''' Test class init

        Parameters
        ----------
        package : string or module
            If string, gives full path to package
            If None, extract calling module path
            Default is None
            
        '''
        if package is None:
            f = sys._getframe(1)
            package = f.f_locals.get('__file__', None)
            assert package is not None
            package = os.path.dirname(package)
        elif isinstance(package, type(os)):
            package = os.path.dirname(package.__file__)
        self.package_path = package
        
    def test(self, label='fast', verbose=1, doctests=False, extra_argv=None):
        ''' Module testing function

        Parameters
        ----------
        label : {'fast', 'full', '', attribute identifer}
            Identifies tests to run.  This can be a string to pass to
            the nosetests executable with the'-A' option, or one of
            several special values.
            Special values are:
            'fast' - the default - which corresponds to
                nosetests -A option of
                'not slow and not bench and not willfail'.
            'full' - fast (as above) and slow tests as in
                nosetests -A option of 'not bench and not willfail'. 
            None or '' - run all tests and benchmarks
            attribute_identifier - string passed directly to
                nosetests as '-A' 
        verbose : integer
            verbosity value for test outputs, 1-10
        doctests : boolean
            If True, run doctests in module, default False
        extra_argv : list
            List with any extra args to pass to nosetests
        '''
        argv = ['scipy module test', self.package_path, '-s']
        if label:
            if not isinstance(label, basestring):
                raise TypeError, 'Test selection label should be a string'
            if label == 'fast':
                label = 'not slow and not bench and not willfail'
            elif label == 'full':
                label = 'not bench and not willfail'
            argv += ['-A', label]
        argv += ['--verbosity', str(verbose)]
        if doctests:
            argv+=['--with-doctest']
        if extra_argv:
            argv+= extra_argv
        nose.run(argv=argv)

        
