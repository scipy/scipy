""" Load or save values to a file.  

    Shelves work well for storing data, but they are slow to access
    repeatedly - especially for large data sets.  This module allows
    you to store data to a file and then load it back into the workspace.
    When the data is stored, a python module is also created as the
    "namespace for the data"
    >>> import data_store
    >>> import os
    >>> a = 1
    >>> data_store.save('c:/temp/junker',{'a':a})
    >>> os.chdir('c:/temp')
    >>> import junker
    >>> junker.a
    1
"""

__all__ = ['load', 'save', 'create_module', 'create_shelf']
import dumb_shelve
import string
import os

def load(module):
    """ Load data into module from a shelf with
        the same name as the module.
    """
    dir,filename = os.path.split(module.__file__)
    filebase = string.split(filename,'.')[0]
    fn = os.path.join(dir, filebase)
    f = dumb_shelve.open(fn, "r")
    #exec( 'import ' + module.__name__)
    for i in f.keys():
        exec( 'import ' + module.__name__+ ';' +
              module.__name__+'.'+i + '=' + 'f["' + i + '"]')
#       print i, 'loaded...'
#   print 'done'    

def save(file_name=None,data=None):
    """ Save the dictionary "data" into
        a module and shelf named save
    """
    import dumb_shelve
    create_module(file_name)
    create_shelf(file_name,data)

def create_module(file_name):
    """ Create the module file.
    """
    if not os.path.exists(file_name+'.py'): # don't clobber existing files
        module_name = os.path.split(file_name)[-1]
        f = open(file_name+'.py','w')   
        f.write('import scipy.io.data_store as data_store\n')
        f.write('import %s\n' % module_name)
        f.write('data_store.load(%s)' % module_name)
        f.close()

def create_shelf(file_name,data):
    """Use this to write the data to a new file
    """
    shelf_name = string.split(file_name,'.')[0]
    f = dumb_shelve.open(shelf_name,'w')
    for i in data.keys():
#       print 'saving...',i
        f[i] = data[i]
#   print 'done'
    f.close()



    
