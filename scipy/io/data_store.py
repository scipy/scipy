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

__all__ = ['save_as_module']

import dumb_shelve
import os


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


def save_as_module(file_name=None,data=None):
    """
    Save the dictionary "data" into a module and shelf named save.

    Parameters
    ----------
    file_name : str, optional
        File name of the module to save.
    data : dict, optional
        The dictionary to store in the module.

    """
    _create_module(file_name)
    _create_shelf(file_name,data)

