#
# Title: Defines wxBackgroundApp used by wxPython_thread.
#
# Author: Pearu Peterson <pearu@cens.ioc.ee>
# Created: October, 2003
#

__all__ = ['wxBackgroundApp']

import sys
import thread
import threading
import traceback
import atexit

wxPython_thread_id = thread.get_ident()
from wxPython import wx
wxEVT_PROXY_EVENT = wx.wxNewId()

class ProxyEvent(wx.wxPyEvent):
    def __init__(self, data, finished):
        wx.wxPyEvent.__init__(self)
        self.SetEventType(wxEVT_PROXY_EVENT)
        self.data = data
        self.finished = finished

class EventCatcher(wx.wxFrame):

    def __init__(self, app):
        wx.wxFrame.__init__(self, wx.NULL,-1,'',size=(300,300))
        self.Connect(-1, -1, wxEVT_PROXY_EVENT, self.proxy_handler)
        wx.EVT_CLOSE(self, self.OnCloseWindow)
        self.app = app

    def proxy_handler(self,evt):
        assert wxPython_thread_id==thread.get_ident(),\
               'wrong thread (not wxPython)'
        code,frame = evt.data
        try:
            exec (code,frame.f_globals,frame.f_locals)
        except Exception:
            type, value, tb = sys.exc_info()
            info = traceback.extract_tb(tb)
            filename, lineno, function, text = info[-1]
            tmp = "\nFile \"%(filename)s\", line %(lineno)s, "\
                  "in %(function)s\n%(text)s"% locals()
            evt.finished.exception_info = (type,str(value) + tmp)
            type = value = tb = None
        evt.finished.set()

    def OnCloseWindow(self,evt=None):
        self.app.ExitMainLoop()
        if evt is not None:
            evt.Skip()

class ExecThread(threading.Thread):
    """ Execute command in a new thread."""
    def __init__(self,cmd,globs,locs):
        threading.Thread.__init__(self)
        self.cmd = cmd
        self.globs = globs
        self.locs = locs
        self.start()
    def run(self):
        exec self.cmd in self.globs,self.locs

class wxBackgroundApp(wx.wxPySimpleApp):

    def __init__(self):
        """ Construct wxPython application
        Must call it from wxPython thread.
        """
        assert wxPython_thread_id==thread.get_ident(),\
               'wrong thread (not wxPython)'
        wx.wxPySimpleApp.__init__(self)
        self.run_is_ready = threading.Event()
        atexit.register(self._shutdown)

    def __call__(self, code, frame=None):
        """ Execute code in wxPython thread within a frame.
        Don't call it from wxPython thread.
        """
        assert wxPython_thread_id!=thread.get_ident(),\
               'wrong thread (wxPython)'
        
        if frame is None:
            frame = sys._getframe(1)

        self.run_is_ready.wait()

        finished = threading.Event(1)
        cmd = 'wx.wxPostEvent(self.event_catcher,'\
              'ProxyEvent((code, frame),finished))'
        ExecThread(cmd,globals(),locals())
        finished.wait(0.5)

        exc_info = getattr(finished,'exception_info',None)
        if exc_info is not None:
            raise exc_info[0],exc_info[1]

    def run(self):
        """ Start wxPython application.
        Must call it from wxPython thread.
        """
        assert wxPython_thread_id==thread.get_ident(),\
               'wrong thread (not wxPython)'
        self.event_catcher = EventCatcher(self)
        self.run_is_ready.set()

        self.MainLoop()

        del self.event_catcher
        self.run_is_ready.clear()

    def _shutdown(self):
        # Note that __init__ registers shutdown in atexit.
        if hasattr(self,'event_catcher'):
            self('self.event_catcher.OnCloseWindow()')
