""" Routines executed in main thread by gui_thread.

    o Starting up wxPython in secondary thread.

    o Shutting down the wxPython thread on system exit.

    o Generating proxy classes for wxPython window
      objects.
      
    This file does not contain an import of wxPython.
    Anything requiring the import of wxPython for its 
    definition should be placed in gui_thread_guts to 
    make sure the import occurs in the second thread.    
"""

###########################
# Secondary Thread Start-up
###########################
import thread, threading
import sys, os, new

running_in_second_thread = 0
gui_thread_finished = threading.Event()

def gui_thread(finished):
    """ Indirectly imports wxPython into the second thread
    """
    try:
        # If we can find a module named wxPython.  Odds are (maybe 100%),
        # we don't want to start a new thread with a MainLoop() in it.
        if not sys.modules.has_key('wxPython'):
            #import must be done inside if statement!!!
            from gui_thread_guts import second_thread_app
            # Variable used to see if the wxApp is
            # running in the main or secondary thread
            # Used to determine if proxies should be generated.            
            global running_in_second_thread,app,gui_thread_finished   
            app = second_thread_app(0)
            running_in_second_thread = 1
            try:
                app.MainLoop()
                # when the main loop exits, we need to single the
                # exit_gui_thread function that it is OK to shut down.
            finally:
                gui_thread_finished.set()
    finally: 
        finished.set()             
        
def start():                    
    finished = threading.Event()
    t1 = threading.Thread(target=gui_thread,args=(finished,))
    t1.setDaemon(1)
    t1.start()
    #!  I'd like to block waiting until I'm sure wxPython
    #!  has been initialized in the second thread.  Unfortunately,
    #!  there is some strange relationship between "import"
    #!  and thread synchronization.  Without this method call,
    #!  the wxPython is import fine.  If you uncomment it,
    #!  however, the import (actually any import) in the
    #!  !second! thread blocks indefinitely.  Remove the
    #!  import and the code executes as expected.
    #!  This issue is documented in the thread_tests directory.
    #finished.wait()
   
###########################
# Secondary Thread Clean-up
###########################
oldexitfunc = getattr(sys, 'exitfunc', None)
def exit_gui_thread(last_exit = oldexitfunc):    
    # don't wait on MS platforms -- it hangs.
    # On X11, we have to shut down the secondary thread.
    if running_in_second_thread and os.name != 'nt':
    	import gui_thread_guts    	
    	event_poster = gui_thread_guts.proxy_base()
    	event_catcher = event_poster.catcher
        finished = threading.Event()
        evt = gui_thread_guts.proxy_event(event_catcher.Close,
                                          (),{},finished)
        event_poster.post(evt)
        # wait for event to get handled
        finished.wait()
        # wait for the gui_thread to die.        
        gui_thread_finished.wait()
    if last_exit: last_exit()
       
sys.exitfunc = exit_gui_thread        

##############################################
# Registering classes and generating proxies.
##############################################

def register(wx_class):
    """ Create a gui_thread compatible version of wx_class
    
        Test whether a proxy is necessary.  If so,
        generate and return the proxy class.  if not,
        just return the wx_class unaltered.       
    """    
    if running_in_second_thread:
        #print 'proxy generated'
        return proxify(wx_class)
    else:
        wx_class.init2 = wx_class.__init__
        wx_class.__init__ = plain_class__init__
        return wx_class

def plain_class__init__(self,*args,**kw):
    self.init2(*args,**kw)
    add_close_event_handler(self)
    self.proxy_object_alive = 1
    
def proxify(wx_class):
    """ Create a proxy class for a wx_class.
    
        It searches for all the unique methods in wx_class
        and its bases classes and creates proxies for
        these calls.  The proxies fire events that
        are handled by the real wx_class object living
        in the wxPython thread.        
    """
    # generate proxy
    from gui_thread_guts import proxy_base
    class_name = 'proxy_'+wx_class.__name__
    class_bases = (proxy_base,)
    class_dict = {}
    class_doc = ""
    try:
        class_doc = getattr(wx_class, '__doc__')
    except AttributeError:
        pass
    else:
        class_dict['__doc__'] = class_doc
    class_methods = get_all_methods(wx_class)
    for method in class_methods:
        func = generate_method(method,wx_class)
        if func:
            class_dict[method] = func
    proxy_class = new.classobj(class_name,class_bases,class_dict)
    return proxy_class

def add_close_event_handler(proxy_obj):
    """ Add an event handler to the wxPython object being wrapped by proxy_obj
        that temporarily catches the Close event.  This is necessary to 
        notify the proxy_obj when the wxPython object is dead.
    """
        
    import gui_thread_guts
    close_handler = gui_thread_guts.CloseEvtHandler(proxy_obj)
    try:
        proxy_obj.wx_obj.PushEventHandler(close_handler)
    except AttributeError:
        # its a standard class that just needs to let us
        # know when it dies.    
        proxy_obj.PushEventHandler(close_handler)
        
def generate_method(method,wx_class):
    """ Create a proxy method.
    
        This first creates a text version of the method, 
        accounting for slight differences between __init__
        methods and all other methods.  It converts the
        text to a code object (using exec) and returns the
        code object.  The code is never actually written to
        a file.  
        
        It takes about .4 seconds on a wxFrame object with
        150 methods.  This is a one time cost at start up.
        It might be beneficial if we use the same code over
        and over to actually write the proxy class to a 
        module. (.pyc file?)        
    """
    module_name = wx_class.__module__
    class_name = wx_class.__name__
    import_statement = 'from %s import %s' % (module_name,class_name)    
    documentation = ""
    try:
        documentation = getattr(getattr(wx_class, method), '__doc__')
    except AttributeError:
        pass
    if method == '__init__':
        call_method = class_name
        pre_test = ''
        #pre_test = 'from gui_thread_guts import proxy_base;'\
        #           'proxy_base.__init__(self)'
        arguments = 'arg_list = args'
        results  = 'self.wx_obj = finished._result;' \
                   'add_close_event_handler(self);' \
                   'self.proxy_object_alive = 1;'
    elif method == '__getattr__':
        return None
    else:
        pre_test =  "if not self.proxy_object_alive: proxy_error()"
        call_method = '%s.%s' % (class_name,method)
        arguments = 'arg_list = tuple([self.wx_obj] + list(args))'
        results  = 'return smart_return(finished._result, self)'
    body = """def %(method)s(self,*args,**kw):
                \"\"\"%(documentation)s\"\"\"
                %(pre_test)s
                from gui_thread_guts import proxy_event, print_exception, smart_return
                %(import_statement)s #import statement
                finished = threading.Event()
                # remove proxies if present
                args = dereference_arglist(args)                
                %(arguments)s #arguments
                evt = proxy_event(%(call_method)s,arg_list,kw,finished)
                self.post(evt)
                finished.wait()
                if finished.exception_info:
                    print_exception(finished.exception_info)
                    raise finished.exception_info['type'],finished.exception_info['value']
                %(results)s #results\n""" %locals()
    #if method == '__init__':
    #    print body
    exec(body)
    return eval(method)

#---------------------- utility functions ----------------------
import types
def is_method(x,cls):
    return type(cls.__dict__[x]) == types.FunctionType
    
def class_methods(cls):
    contents = dir(cls)
    methods = filter(lambda x,cls=cls: callable(getattr(cls, x)), contents)
    # using callable since it will handle more than just functions,
    # using getattr since it works for vtk objects too
    # (cls.__dict__[x]) won't work for vtk objects
    return methods

def get_all_methods(wx_class,top=1):
    assert (wx_class) is types.ClassType or hasattr(wx_class, '__bases__')
    methods = []
    for base in wx_class.__bases__:
        methods.extend(get_all_methods(base,0))
    methods.extend(class_methods(wx_class))    
    if top:
        methods = remove_duplicates(methods)
    return methods        

def remove_duplicates(lst):
    """ Discard any methods that have equivalent names
    """
    res = []    
    for x in lst:
        if not x in res: 
            res.append(x)
    return res

def is_proxy(x):
    return hasattr(x,'is_proxy')
    
def dereference_arglist(lst):
    """ Scan for proxy objects and convert to underlying object
    """
    res = []
    for arg in lst:
        if is_proxy(arg): 
            res.append(arg.wx_obj)
            #print 'dereferenced ', arg.wx_obj
        else: res.append(arg)            
    return res
       
def proxy_error():
    raise ValueError, 'This window has been destroyed'
