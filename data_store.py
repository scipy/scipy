""" Load or save values to a file.

    Shelves work well for storing data, but they are slow to access
    repeatedly - especially for large data sets.  This module allows
    you to store data to a file and then load it back into the workspace.
    When the data is stored, a python module is also created as the
    "namespace for the data"
    >>> import scipy.io
    >>> import os
    >>> a = 1
    >>> scipy.io.save_as_module('c:/temp/junker',{'a':a})
    >>> os.chdir('c:/temp')
    >>> import junker
    >>> junker.a
    1
"""

__all__ = ['save_as_module',
           # The rest of these are all deprecated
           'save', 'create_module',
           'create_shelf', 'load']

import dumb_shelve
import os

from numpy import deprecate_with_doc, deprecate

def _load(module):
    """ Load data into module from a shelf with
        the same name as the module.
    """
    dir,filename = os.path.split(module.__file__)
    filebase = filename.split('.')[0]
    fn = os.path.join(dir, filebase)
    f = dumb_shelve.open(fn, "r")
    #exec( 'import ' + module.__name__)
    for i in f.keys():
        exec( 'import ' + module.__name__+ ';' +
              module.__name__+'.'+i + '=' + 'f["' + i + '"]')
#       print i, 'loaded...'
#   print 'done'

load = deprecate_with_doc("""
This is an internal function used with scipy.io.save_as_module

If you are saving arrays into a module, you should think about using
HDF5 or .npz files instead.
""")(_load)


def _create_module(file_name):
    """ Create the module file.
    """
    if not os.path.exists(file_name+'.py'): # don't clobber existing files
        module_name = os.path.split(file_name)[-1]
        f = open(file_name+'.py','w')
        f.write('import scipy.io.data_store as data_store\n')
        f.write('import %s\n' % module_name)
        f.write('data_store._load(%s)' % module_name)
        f.close()

create_module = deprecate_with_doc("""
This is an internal function used with scipy.io.save_as_module

If you are saving arrays into a module, you should think about
using HDF5 or .npz files instead.
""")(_create_module)

def _create_shelf(file_name,data):
    """Use this to write the data to a new file
    """
    shelf_name = file_name.split('.')[0]
    f = dumb_shelve.open(shelf_name,'w')
    for i in data.keys():
#       print 'saving...',i
        f[i] = data[i]
#   print 'done'
    f.close()

create_shelf = deprecate_with_doc("""
This is an internal function used with scipy.io.save_as_module

If you are saving arrays into a module, you should think about using
HDF5 or .npz files instead.
""")(_create_shelf)


def save_as_module(file_name=None,data=None):
    """ Save the dictionary "data" into
        a module and shelf named save
    """
    _create_module(file_name)
    _create_shelf(file_name,data)

save = deprecate(save_as_module, 'save', 'save_as_module')
