#
# Title: Import wxPython to its own thread.
#
# Author: Pearu Peterson <pearu@cens.ioc.ee>
# Created: October, 2003
#

__all__ = ['wxPython_thread']

import re
import os
import os.path
import sys
import new
import types
import thread
import inspect
import atexit
import imp
from scipy_base import ParallelExec

def get_extmodules(module,pexec):
    # .. that need wrappers
    if module != 'wx': # The wx module is a dummy that we create.
        assert not sys.modules.has_key(module),\
               module+' is already imported, cannot proceed'
    state0 = sys.modules.keys()
    pexec('import '+module,wait=1)
    state1 = sys.modules.keys()
    return [k for k in state1 if k not in state0 \
            and (not hasattr(sys.modules[k],'__file__') \
                 or os.path.splitext(sys.modules[k].__file__)[1][:3]!='.py')]

def wrap_extmodule(module, call_holder):
    """ Return wrapper module.
    Current design assumes wx extension modules.
    """
    new_module = new.module(module.__name__)
    for m in inspect.getmembers(module):
        n,v = m
        t = type(v)
        if t in [types.NoneType,types.StringType,types.IntType]:
            # No wrappers are needed.
            setattr(new_module,n,v)
        elif t is types.BuiltinFunctionType:
            setattr(new_module,n,wrap_builtinfunction(v,call_holder))
        elif t is types.ClassType:
            if str(v)=='wxPython.wxc.wxPyAssertionError':
                # Do we need a wrapper for it?
                setattr(new_module,n,v)
            else:
                raise NotImplementedError,`t`
        elif n in ('wxPyCoreAPI', '_wxPyCoreAPI'):
            # wxPyCoreAPI is not used by the python part of wx, so
            # no wrapper is necessary
            setattr(new_module,n,v)
        elif n=='cvar':
            # In wx all cvar attributes are of string type:
            #print [type(getattr(v,a)) for a in get_cvar_attrs(v)]
            # So no wrappers are necessary.
            setattr(new_module,n,v)
        else:
            raise NotImplementedError,`t`
    # This gives the way how to recognize if module is a wrapper:
    #new_module._old_module = module 
    return new_module

def wrap_builtinfunction(func, call_holder):
    func_name = func.__name__
    if func_name in ('wxApp_CleanUp', 'App_CleanUp'):
        return func
    main_thread_id = call_holder.main_thread_id
    func_tmpl = """\
def %(func_name)s(*args,**kws):
    if get_ident()==%(main_thread_id)s:
        call_holder.call('call_holder.result=old_func(*args,**kws)')
        result = call_holder.result
        call_holder.result = None
        return result
    return old_func(*args,**kws)
func_code=%(func_name)s.func_code"""
    exec func_tmpl % (locals())
    globs = {'old_func':func,
             'get_ident':thread.get_ident,
             'call_holder':call_holder}
    return new.function(func_code,globs,func_name)

## if sys.version[:3]>='2.3':
##     from tempfile import mkstemp
##     mktemp = lambda :mkstemp()[1]
## else:
##     from tempfile import mktemp
## from cStringIO import StringIO
## re_cvar_attr = re.compile(r'Global variables {(?P<attrs>[\s,\w\W]*)}')
## cvar_attr_tmp_fn = mktemp()
## def get_cvar_attrs(cvar):
##     f = open(cvar_attr_tmp_fn,'w')
##     print >>f,cvar
##     f.close()
##     s = open(cvar_attr_tmp_fn,'r').read()
##     m = re_cvar_attr.match(s)
##     assert m is not None,`(`s`,cvar)`
##     return m.group('attrs').strip().split(', ')

def _import_wx_core(wx_pth, pexec):
    """Imports the core modules for wx.  This is necessary for
    wxPython-2.5.x. 
    """
    # Find the suffix.
    suffix = '.so'
    for x in [x[0] for x in imp.get_suffixes() if x[-1] is imp.C_EXTENSION]:
        if os.path.exists(os.path.join(wx_pth, '_core_' + x)):
            suffix = x
            break
    # Now import the modules manually.
    pexec('import imp, os.path')
    code="""\
for i in ["_core_", "_controls_", "_misc_", "_windows_", "_gdi_"]:
    p = os.path.join('%s', i + '%s')
    imp.load_dynamic('wx.' + i, p)
"""%(wx_pth, suffix)
    pexec(code)
    
    # Now create a dummy module in sys.modules to inhibit importing the
    # actual one.  We will reload(wx) to get it right later.
    m = new.module('wx')
    m.__file__ = os.path.join('wx_pth', '__init__.py')
    sys.modules['wx'] = m

    
def wxPython_thread():
    class AttrHolder: pass

    call_holder = AttrHolder()
    call_holder.call = None
    call_holder.main_thread_id = thread.get_ident()
    pexec = ParallelExec()

    # Check if 'wx' namespace based modules are used.
    mod_name = 'wxPython'
    for path in sys.path:
        wx_pth = os.path.join(path, 'wx')
        if os.path.exists(os.path.join(wx_pth, '__init__.py')):
            assert not sys.modules.has_key('wx'), \
                   'wx is already imported, cannot proceed'
            mod_name = 'wx'
            _import_wx_core(wx_pth, pexec)
            break
    
    # Create wrappers to wxPython extension modules:
    for name in get_extmodules(mod_name,pexec):
        module = sys.modules[name]
        if module is None: # happens with gui_thread.wxPython
            continue
        sys.modules[name] = wrap_extmodule(module,call_holder)
    
    if mod_name == 'wx':
        import wx
        pexec('reload(wx)')
    else:
        import wxPython
        pexec('reload(wxPython')

    pexec('from wxBackgroundApp import wxBackgroundApp')
    pexec('call_holder.call=wxBackgroundApp()',wait=1)
    # Start wxPython application in background:
    pexec('call_holder.call.run()')

    atexit.register(pexec.shutdown)
    return
