"""
    This module is separated from gui_thread because it imports
    wxPython.  We can't do that at the top level of gui_thread.
    It must be done in the wxPython thread.  Classes that need
    wxPython classes in there definitions (inheritance) are 
    banished to this module.
    
    If someone comes up with a way to dump this into gui_thread
    that isn't ugly, we may consider doing that.  I'm sure
    there is a way.  It's the 'isn't ugly' part that concerns me...
    
    The file has the following contents:
    
    o An event handler that is chained onto any proxified window 
      to notifies the wrapping proxy object when the window has 
      been closed.
    
    o Proxy Event class used to pass requests from the main thread
      to the wxPython thread
    
    o A base class for all proxy objects called proxy_base

    o A proxied_callable class that handles callable attributes of a
      proxied clas by making the call secondary thread.

    o A proxy_attr class that handles attributes of a proxied object
      correctly by dispatching calls to the secondary thread.
    
    o A wxPython class called event_catcher that lives in the
      secondary thread and catches/dispatches all events requested
      from the primary thread
    
    o The main wxApp that lives in the secondary thread (very simple)
      
"""
import sys
_in_thread = 0
if not sys.modules.has_key('wxPython'):
    _in_thread = 1
    print '<Importing wxPython... ',
    sys.stdout.flush()
from wxPython.wx import *
if _in_thread:
    print 'done.>\n>>> ',
    sys.stdout.flush()
del _in_thread

import thread, threading
import types, traceback
import UserList, UserDict
import main

#################################
# Window Close Event Handler 
#################################        

class CloseEvtHandler(wxEvtHandler):
    """ An Event Handler added to all proxified objects to notify
        the proxy wrapper when the object has been closed.    
    """
    def __init__(self,proxy_obj):
        """ Build an event handler to handle the Close event.
            It needs access to the proxy object to tell it when
            the underlying wxPython window has been closed.
        """
        wxEvtHandler.__init__(self)
        self.proxy_obj = proxy_obj
        EVT_CLOSE(self,self.ProxyOnCloseWindow)
    def ProxyOnCloseWindow(self,evt):        
        """ Tell the proxy wrapper that the object is no longer
             alive.  Logic in the wrapper methods use this info
             to make sure they don't call a dead wxPython object.   
        """
        try:
            self.proxy_obj.kill_proxy()
        except AttributeError:
            self.proxy_obj.proxy_object_alive = 0
        del self.proxy_obj
        evt.Skip()

#################################
# Proxy Event Object
#################################        

wxEVT_PROXY_EVENT = 26001
def PROXY_EVENT(win, func):
    win.Connect(-1, -1, wxEVT_PROXY_EVENT, func)    
        
class proxy_event(wxPyEvent):
    """ Proxy events transport a method code
        object, the arguments to this method
        and a synchronization method from one
        thread to another.
    """
    def __init__(self,method,args,kw,finished):    
        """ Create a proxy_event object.
        
            method -- a code object that is to be called
            args   -- list of arguments.  if the call is
                      to an unbound method, the first entry
                      should be a class object of the 
                      appropriate type.
            kw     -- dictionary of key word arguments.
            finished -- a threading.Event() object used
                        for syncrhonization. 
        """
        wxPyEvent.__init__(self)
        self.SetEventType(wxEVT_PROXY_EVENT)
        self.method = method
        self.args = args
        self.kw = kw            
        self.finished = finished


#################################
# A proxied callable object
#################################

def is_immutable(x):
    """ Checks if object is completely immutable.  A tuple is not
    considered completely immutable since it could contain references
    to objects that could change.  Returns 1 if it is immutable or 0
    if object is mutable."""    
    imm = ()
    try:
        imm = (types.StringType, types.FloatType, types.IntType,
               types.ComplexType, types.NoneType, types.UnicodeType)
    except AttributeError:
        imm = (types.StringType, types.FloatType, types.IntType,
               types.ComplexType, type.NoneType)
    if type(x) in imm:
        return 1
    else:
        return 0


def smart_return(ret, proxy):    
    """ This intelligently returns an appropriately proxied object to
    prevent crashing the interpreter.  If the object is immutable it
    simply returns it.  If it is callable it returns a
    proxied_callable and if it is niether it returns a proxy_attr."""
    if callable(ret):
        return proxied_callable(proxy, ret)
    elif is_immutable(ret):
        return ret
    elif type(ret) in (types.ListType, types.TupleType):
        return proxy_list(proxy, ret)
    elif type(ret) is types.DictType:
        return proxy_dict(proxy, ret)
    else:
        return proxy_attr(proxy, ret)


class proxied_callable:

    """This wraps any callable object so that the call is made via a
    proxy using a proxy_event.  This makes it possible for the user to
    call methods of a proxy's attribute.  """
    
    def __init__(self, proxy, call_obj):        
        """ Create the proxied callable object.

            proxy -- the proxied object whose attribute is being
            proxied.

            call_obj -- the callable attribute that is being proxied.
        """
        self.proxy = proxy
        self.__dont_mess_with_me_unless_you_know_what_youre_doing = call_obj
        try:
            call_obj.__doc__
        except AttributeError:
            pass
        else:
            self.__doc__ = call_obj.__doc__
            
    def __call__(self, *args, **kw):        
        """Performs the call to the proxied callable object by
        dispatching the method to the secondary thread."""        
        obj = self.__dont_mess_with_me_unless_you_know_what_youre_doing
        finished = threading.Event()
        evt = proxy_event(obj, args, kw, finished)
        self.proxy.post(evt)
        finished.wait()
        if finished.exception_info:
            raise finished.exception_info[0], \
                  finished.exception_info[1]

        return smart_return(finished._result, self.proxy)


################################################################
# A proxied list that proxies lists and tuples safely and yet
# behaves very much like them.
################################################################

class proxy_list(UserList.UserList):
    def __init__(self, proxy, val):
        UserList.UserList.__init__(self, val)
        self.proxy = proxy

    def __getitem__(self, i):
        val = self.data[i]
        return smart_return(val, self.proxy)

    def __getslice__(self, i, j):
        i = max(i, 0); j = max(j, 0)
        return self.__class__(self.proxy, self.data[i:j])


################################################################
# A proxied list that proxies dictionaries safely and yet behaves
# very much like one.
################################################################

class proxy_dict(UserDict.UserDict):
    def __init__(self, proxy, val):
        UserDict.UserDict.__init__(self, val)
        self.proxy = proxy

    def __getitem__(self, key):
        val = self.data[key]
        return smart_return(val, self.proxy)

    def keys(self):
        return smart_return(self.data.keys(), self.proxy)

    def items(self):
        return smart_return(self.data.items(), self.proxy)

    def values(self):
        return smart_return(self.data.values(), self.proxy)

    def get(self, key, failobj=None):
        val = self.data.get(key, failobj)
        return smart_return(val, self.proxy)

    def setdefault(self, key, failobj=None):
        if not self.data.has_key(key):
            self.data[key] = failobj
        return smart_return(self.data[key], self.proxy)

    def popitem(self):
        return smart_return(self.data.popitem(), self.proxy)


################################################################
# A proxied attribute that handles callable internal attrubutes
# properly.
################################################################

class proxy_attr:

    """
    Description

      This wraps a proxy's attribute such that a call to the
      attribute's member function returns a proxied_callable and any
      request for a data attribute will return another proxy_attr.
      This way a user can call a method of a proxy's attribute.  Use
      the get_members and get_methods member functions to get a list
      of the data and functions of the proxied object.

    Caveats

      The user can always use
      _proxy_attr__dont_mess_with_me_unless_you_know_what_youre_doing
      and still get around proxying but atleast they were warned. :)
      Also, this isnt a normal Python object in that
      self.__class__.__dict__ will not contain any methods of the
      proxied object.  Only self.__dict__ will contain these.  It
      works but isn't perfect.

    """
    
    __DONT_WRAP = ('__del__', '__init__', '__getattr__')

    def __init__(self, proxy, obj):
        """ Create the proxied attribute.

            proxy -- the proxied object whose attribute is being
            proxied.

            obj -- the attribute that is being proxied.
        """
        self.__proxy = proxy
        self.__dont_mess_with_me_unless_you_know_what_youre_doing = obj
        if hasattr(obj, '__class__'):
            try:
                doc = obj.__class__.__doc__
                self.__dict__['__doc__'] = doc
            except AttributeError:
                pass
            self._set_attrs(obj, main.get_all_methods(obj.__class__))
        else:
            # this is most probably a built in type. but in any case
            self._set_attrs(obj, dir(obj))

    def _get_meth (self, obj, name):
        if name in self.__DONT_WRAP:
            return None
        try:
            ret = getattr(obj, name)
        except AttributeError:
            return None
        if callable(ret):
            return proxied_callable(self.__proxy, ret)
        elif is_immutable(ret):
            return ret
        else:
            return None

    def _set_attrs(self, obj, names):
        for i in names:
            tmp = self._get_meth(obj, i)
            if tmp:
                self.__dict__[i] = tmp

    def get_members(self):
        """Returns a list of the members of the proxied attribute."""
        obj = self.__dont_mess_with_me_unless_you_know_what_youre_doing
        return dir(obj)

    def get_methods(self):
        """Returns a list of the methods of the proxied attribute."""
        obj = self.__dont_mess_with_me_unless_you_know_what_youre_doing
        try:
            return dir(obj.__class__)
        except AttributeError:
            return dir(obj)

    def __getattr__(self, key):
        obj = self.__dont_mess_with_me_unless_you_know_what_youre_doing
        ret = getattr(obj, key)
        return smart_return(ret, self.__proxy)

#################################
# Base class for all automatically generated proxy classes
#################################        

class proxy_base:
    is_proxy = 1
    # all proxy objects share the same event catcher.
    # This is set to be a single instance of the
    # event_catcher class when the second thread is 
    # started.
    catcher = None
    
    def kill_proxy(self):
        """ Used by the CloseEvtHandler to notify
            proxy class that underlying wxPython window
            has been closed.
        """
        self.proxy_object_alive = 0
        
    def post(self,evt):
        """ All proxy classes send their events
            to a single "catcher" window using
            this method.  The catcher is set at 
            startup in the OnInit method of 
            second_thread_app.
        """
        wxPostEvent(self.catcher, evt)        

    def __getattr__(self,key):        
        try:
            ret = self.__dict__[key]
        except KeyError:
            ret = getattr(self.__dict__['wx_obj'],key)
        return smart_return(ret, self)
    
    # This needs a little thought
    #def __setattr__(self,key,val):        
    #    finished = threading.Event()
    #    evt = proxy_event(setattr,(self.wx_obj,key,val),finished)
    #    self.post_event(evt)
    #    finished.wait()
    #    if finished.exception_info:
    #        raise finished.exception_info[0], finished.exception_info[1]
    #    return finished._result

#################################
# Receives/Dispatches all events requested from the main thread.
#################################        

class event_catcher(wxFrame):
    """ The "catcher" frame in the second thread.
        It is invisible.  It's only job is to receive
        proxy_vents from the main thread, and then call
        the method indicated in the event.
    """ 
    def __init__(self):
        wxFrame.__init__(self, NULL,-1,'',size=(300,300))        
        EVT_CLOSE(self,self.OnCloseWindow)
        PROXY_EVENT(self,self.proxy_handler)    
    def proxy_handler(self,evt): 
        try:
            evt.finished._result = apply(evt.method,evt.args,evt.kw)
            evt.finished.exception_info = None
            evt.finished.set()
        except:
            type, value, tb = sys.exc_info()
            info = traceback.extract_tb(tb)
            filename, lineno, function, text = info[-1] # last line only
            tmp = "\nFile \"%(filename)s\", line %(lineno)s, "\
                  "in %(function)s\n%(text)s"%locals()
            value = str(value) + tmp
            evt.finished.exception_info = [type, value]
            type = value = tb = None # clean up
            evt.finished.set()

    def OnCloseWindow(self,evt):
        main.app.ExitMainLoop()

#################################
# The (very simple) wxApp that lives in secondary thread
#################################        
                                
class second_thread_app(wxApp):
    """  wxApp that lives in the second thread.
    """    
    def OnInit(self):
        """ Registers the event_catcher window as 
            the object to whom proxy classes should 
            send their events.
        """    
        proxy_base.catcher = event_catcher()
        return true

