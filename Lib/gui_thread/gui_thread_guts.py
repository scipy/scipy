"""
    This module is separated from gui_thread because it imports
    wxPython.  We can't do that at the top level of gui_thread.
    It must be done in the wxPython thread.  Classes that need
    wxPython classes in there definitions (inheritance) are 
    banished to this module.
    
    If someone comes up with a way to dump this into gui_thread
    that isn't ugly, we may consider doing that.  I'm sure
    there is a way.  It's "isn't ugly" part that concerns me...
    
    The file has the following contents:
    
    o An event handler that is chained onto any proxified window 
      to notifies the wrapping proxy object when the window has 
      been closed.
    
    o  Proxy Event class used to pass requests from the main thread
       to the wxPython thread
    
    o A base class for all proxy objects called proxy_base
    
    o A wxPython class called event_catcher that lives in the
      secondary thread and catches/dispatches all events requested
      from the primary thread
    
    o The main wxApp that lives in the secondary thread (very simple)
      
"""

from wxPython.wx import *
print '<wxPython imported>\n>>> ',
import thread, threading

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
        self.proxy_obj.kill_proxy()
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
        from wxPython.wx import wxPostEvent
        wxPostEvent(self.catcher, evt)        

    def __getattr__(self,key):        
        try:
            return self.__dict__[key]            
        except KeyError:
            return getattr(self.__dict__['wx_obj'],key)

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
            #print 'exception occured in the following call:'
            #print '\tmethod:', evt.method
            #print '\targuments:', evt.args
            #print '\tkey words:', evt.kw            
            import sys
            err_type, err_msg = sys.exc_info()[:2]
            #print '\terror:',err_type
            #print '\tmessage:', err_msg

            import sys
            evt.finished.exception_info = sys.exc_info()[:2]           
            evt.finished.set()

    def OnCloseWindow(self,evt):
        import main
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

