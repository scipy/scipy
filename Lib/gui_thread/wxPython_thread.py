#
# Title: Import wxPython to its own thread.
#
# Author: Pearu Peterson <pearu@cens.ioc.ee>
# Created: October, 2003
#

__all__ = ['wxPython_thread']

import os
import sys
import new
import types
import thread
import inspect
import atexit
import imp
import glob
from scipy_base import ParallelExec

def get_extmodules(pexec):
    # .. that need wrappers

    state0 = sys.modules.keys()
    _do_import(pexec)
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


def _do_import(pexec):
    """This imports the wxPython modules.  This code handles both
    wxPython 2.4.x and 2.5.x versions.

    """
    sfx = [x[0] for x in imp.get_suffixes() if x[-1] is imp.C_EXTENSION]
    suffix = '.so'

    # Handle the imports for wxPython-2.5.x
    flag = 0
    for path in sys.path:
        wx_pth = os.path.join(path, 'wx')
        if os.path.exists(os.path.join(wx_pth, '__init__.py')):
            for x in sfx:
                if os.path.exists(os.path.join(wx_pth, '_core_' + x)):
                    suffix = x
                    flag = 1
                    assert not sys.modules.has_key('wx'), \
                           'wx is already imported, cannot proceed'
                    # Now do necessary stuff for wx.
                    _import_wx_core(pexec, wx_pth, suffix)
                    break

    # wxPython-2.4.x
    if not flag: 
        assert not sys.modules.has_key('wxPython'),\
               'wxPython is already imported, cannot proceed'
        pexec('import wxPython', wait=1)
        path = os.path.dirname(sys.modules['wxPython'].__file__)
        suffix = os.path.splitext(sys.modules['wxPython.wxc'].__file__)[1]
        exts = glob.glob(os.path.join(path, '*' + suffix))
        for x in [os.path.splitext(os.path.basename(i))[0] for i in exts]:
            pexec('import wxPython.%s'%x, wait=1)
        

def _import_wx_core(pexec, wx_pth, suffix):
    """Imports the core modules for wx.  This is necessary for
    wxPython-2.5.x.

    """
    # Create a dummy module in sys.modules to inhibit importing the
    # actual one.  We will reload(wx) to get it right later.
    m = new.module('wx')
    m.__file__ = os.path.join('wx_pth', '__init__.py')
    sys.modules['wx'] = m

    # Now import all the extension modules manually.
    pexec('import imp')
    libs = glob.glob(os.path.join(wx_pth, '*' + suffix))
    done = []
    for i in ['_core_', '_controls_', '_misc_', '_windows_', '_gdi_']:
        p = os.path.join(wx_pth, i + suffix)
        pexec('imp.load_dynamic("wx." + "%s", "%s")'%(i, p), wait=1)
        done.append(p)
    libs = [x for x in libs if x not in done]
    for i in libs:
        p = os.path.splitext(os.path.basename(i))[0]
        pexec('imp.load_dynamic("wx." + "%s", "%s")'%(p, i), wait=1)

    
def wxPython_thread():
    class AttrHolder: pass

    call_holder = AttrHolder()
    call_holder.call = None
    call_holder.main_thread_id = thread.get_ident()
    pexec = ParallelExec()

    # Create wrappers to wxPython extension modules:
    for name in get_extmodules(pexec):
        module = sys.modules[name]
        if module is None: # happens with gui_thread.wxPython
            continue
        sys.modules[name] = wrap_extmodule(module,call_holder)
    
    if sys.modules.has_key('wxPython.wxc'):
        pexec('import wxPython; reload(wxPython)', wait=1)
        import wxPython
    else:
        pexec('import wx; reload(wx)', wait=1)
        import wx

    pexec('from wxBackgroundApp import wxBackgroundApp')
    pexec('call_holder.call=wxBackgroundApp()',wait=1)
    # Start wxPython application in background:
    pexec('call_holder.call.run()')

    atexit.register(pexec.shutdown)
    return
