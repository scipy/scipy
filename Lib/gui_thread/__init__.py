"""
    Allows interpreter and GUI windows to co-exist peacefully.

    gui_thread runs wxPython in a back ground thread and allows you 
    to create and manipulate wxPython windows from the command line.

    NEW WAY:
    An example usage follows::
           >> import gui_thread
           >> gui_thread.start()
           >> from gui_thread.examples import SimpleFrame
           >> a = SimpleFrame()
           >> a.Show(1)
           # a window with a button should pop up...
           >> a.SetTitle('bob')
           >> a.GetTitle()
           'bob'
           >> a.GetSize()
           (350,200)
           >> a.Blue()
           # window button should have turned blue.
           >> b = SimpleFrame()
           >> b.Show(1)
           # a 2nd window appears

    OLD WAY:
    An example usage follows::
           >> import gui_thread
           >> gui_thread.start(use_main=1)
           <wait several seconds while wxPython starts>
           wxPython imported           
           >> from gui_thread.examples import SimpleFrame
           >> prxyMyFrame = gui_thread.register(SimpleFrame)
           >> a = prxyMyFrame()
           >> a.Show(1)
           # a window with a button should pop up...
           >> a.SetTitle('bob')
           >> a.GetTitle()
           'bob'
           >> a.GetSize()
           (350,200)
           >> a.Blue()
           # window button should have turned blue.
           >> b = prxyMyFrame()
           >> b.Show(1)
           # a 2nd window appears
        
        Modules based on gui_thread will also run in wxPython applications
        without modification.  So, by using gui_thread, you can create
        viewers for common objects that are usable from the standard
        python interpreter, PythonWin (on MS), and wxPython applications.
        Since Tkinter and wxPython currently do not coexist well together,
        IDLE is not compatible with gui_thread.        
        
        gui_thread should always be the first module you import into an 
        application, unless it is a wxPython app.  In that case, it can
        be imported at any time.  When using it in interactive mode, make
        wait until the "wxPython Imported" message appears before importing
        any other modules that use wxPython.
"""

# New hooks for importing wxPython to its own thread
from wxPython_thread import wxPython_thread

_use_main = 0
def start(extra=None,use_main=0):
    global _use_main
    _use_main = use_main
    if _use_main:
        import main
        main.start(extra)
    else:
        wxPython_thread()
    return

def start_up(src_modname,dst_modname,extra=None):
    global _use_main
    if _use_main:
        import main
        main.start_up(src_modname,dst_modname,extra)
    return

def register(wx_class):
    global _use_main
    if _use_main:
        import main
        return main.register(wx_class)
    return wx_class
