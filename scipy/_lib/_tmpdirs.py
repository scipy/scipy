''' Contexts for *with* statement providing temporary directories
'''
from contextlib import chdir, contextmanager
from tempfile import TemporaryDirectory

@contextmanager
def in_tempdir():
    ''' Create, return, and change directory to a temporary directory

    Examples
    --------
    >>> import os
    >>> my_cwd = os.getcwd()
    >>> with in_tempdir() as tmpdir:
    ...     _ = open('test.txt', 'wt').write('some text')
    ...     assert os.path.isfile('test.txt')
    ...     assert os.path.isfile(os.path.join(tmpdir, 'test.txt'))
    >>> os.path.exists(tmpdir)
    False
    >>> os.getcwd() == my_cwd
    True
    '''
    with TemporaryDirectory() as td:
        with chdir(td):
            yield td
