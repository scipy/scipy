"""
SciPy --- A scientific computing package for Python
===================================================

Available subpackages
---------------------

"""
# Modules to import under the scipy namespace
from scipy_version import scipy_version as __version__
from scipy_base import *
from helpmod import *

from scipy_test.testing import ScipyTest

#-------- doc string hooks --------#

def _level_docs(module=None,_cache=[]):
    if module is None:
        return _cache
    try:
        title = module.__dict__['__doc_title__']
    except KeyError:
        try:
            title = module.__dict__['__doc__']
            title = title.lstrip().split('\n',1)[0]
        except KeyError:
            title = ' N/A'
    _cache.append((module.__name__,title))

def _pkg_titles():
    level_docs = _level_docs()
    lengths = [len(name) for (name,title) in level_docs]
    max_length = max(lengths)
    lines = []
    for (name,title) in level_docs:
        w = max_length - len(name)
        lines.append('%s%s --- %s' % (name, w*' ', title))
    return '\n'.join(lines)

#----- Import packages ----#

def _import_packages():
    """ Import packages in scipy directory that implement
    info_<packagename>.py. See DEVELOPERS.txt for more info.
    """
    from glob import glob
    import os
    frame = sys._getframe(1)

    for info_file in glob(os.path.join(__path__[0],'*','info_*.py')):

        package_name = os.path.basename(os.path.dirname(info_file))
        if package_name != os.path.splitext(os.path.basename(info_file))[0][5:]:
            print '  !! Mismatch of package name %r and %s' \
                  % (package_name, info_file)
            continue

        sys.path.insert(0,os.path.dirname(info_file))
        # TODO: catch exceptions here:
        exec 'import info_%s as info_module' % (package_name)
        del sys.path[0]
        
        if getattr(info_module,'ignore',0):
            continue

        global_symbols = getattr(info_module,'global_symbols',[])
        
        if getattr(info_module,'postpone_import',1):
            code = '%s = ppimport(%r)' % (package_name,package_name)

            for name in global_symbols:
                code += '\n%s = ppimport_attr(%s,%r)' % (name,package_name,name)
        else:
            code = 'import %s' % (package_name)
            # XXX: Should we check the existence of package.test? Warn?
            code += '\n%s.test = ScipyTest(%s).test' % (package_name,package_name)
            if global_symbols:
                code += '\nfrom '+package_name+' import '+','.join(global_symbols)
        # XXX: Should we catch exceptions here??
        exec (code, frame.f_globals,frame.f_locals)

        _level_docs(info_module)
        # XXX: Ugly hack to fix package name:
        code = '_level_docs()[-1] = (%s.__name__,_level_docs()[-1][1])' \
               % (package_name)
        exec (code, frame.f_globals,frame.f_locals)

_import_packages()
del _import_packages

#Packages installed standalone:
weave = ppimport('weave');             _level_docs(weave)
wxPython_thread = ppimport_attr(ppimport('gui_thread'),'wxPython_thread')

#----------------------------------#

from common import *
from pilutil import *

#---- testing ----#

test = ScipyTest('scipy').test

#----- update doc string -------#

__doc__ += _pkg_titles()
