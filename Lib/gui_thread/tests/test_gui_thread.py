""" Testing gui_thread.  Please run this under gui_thread and also
without it.  Right now the testing is extremely rudimentary."""

import unittest, gui_thread, weakref, time, sys

from wxPython.wx import *

def threaded():
    return gui_thread.main.running_in_second_thread

def is_alive(obj):
    if obj() is None:
        return 0
    else:
        return 1

def yield():
    if not threaded():
        # forces the event handlers to finish their work.
        # this also forces deletion of windows.
        wxYield() 
    else:
        time.sleep(0.05) # sync threads
    

class TestPanel(wxPanel):
    def __init__(self, parent):
        wxPanel.__init__(self, parent, -1)
        ID = NewId()
        btn = wxButton(self, ID, "Hello")
        EVT_BUTTON(self, ID, self.OnButton)

    def OnButton(self, evt):
        print "TestPanel.OnButton"


class TestFrame(wxFrame):
    def __init__(self, parent):
        wxFrame.__init__(self, parent, -1, "Hello Test")
        self.panel = TestPanel(self)
        EVT_CLOSE(self, self.OnClose)
        self.Show(1)

    def OnClose(self, evt):
        self.Destroy()


class TestClass:
    def __init__(self):
        self.a = 1

    def test(self):
        print self.a


class test_gui_thread(unittest.TestCase):
    def check_wx_class(self):
        "Checking a wxFrame proxied class"
        for i in range(5):
            f = gui_thread.register(TestFrame)
            a = f(None)
            p = weakref.ref(a)
            a.Close(1)
            del a
            yield()
            # this checks for memory leaks
            self.assertEqual(is_alive(p), 0)
            
    def check_normal_class(self):
        "Checking non-wxWindows proxied class "
        f = gui_thread.register(TestClass)
        a = f()
        p = weakref.ref(a)
        # the reference count has to be 2 since nothing special should
        # be done for these classes.
        self.assertEqual(sys.getrefcount(a), 2)
        del a
        self.assertEqual(is_alive(p), 0)        
        
    def check_exception(self):
        "Checking exception handling"
        f = gui_thread.register(TestFrame)
        a = f(None)
        p = weakref.ref(a)
        self.assertRaises(TypeError, a.Close, 1, 2, 3)
        a.Close()
        del a
        yield()
        # this checks for memory leaks
        self.assertEqual(is_alive(p), 0)

def test_suite():
    suites = []
    suites.append(unittest.makeSuite(test_gui_thread, 'check_'))
    total_suite = unittest.TestSuite(suites)
    return total_suite

def test():
    all_tests = test_suite()
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(all_tests)


############################################################

# Utility clases for running the tests and for beautification.

############################################################

class TesterApp (wxApp):
    def OnInit (self):
        f = TesterFrame(None)
        return true
    
class TesterFrame(wxFrame):
    def __init__(self, parent):
        wxFrame.__init__(self, parent, -1, "Tester")

        self.CreateStatusBar()
        sizer = wxBoxSizer(wxHORIZONTAL)
        ID = NewId()
        btn = wxButton(self, ID, "Start Test")
        EVT_BUTTON(self, ID, self.OnStart)
        msg = "Click to start running tests. "\
              "Tester Output will be shown on the shell."
        btn.SetToolTip(wxToolTip(msg))
        sizer.Add(btn, 1, wxEXPAND)
        ID = NewId()
        btn = wxButton(self, ID, "Close")
        EVT_BUTTON(self, ID, self.OnClose)
        btn.SetToolTip(wxToolTip("Click to close the tester."))
        sizer.Add(btn, 1, wxEXPAND)

        sizer.Fit(self)
        self.SetAutoLayout(true)
        self.SetSizer(sizer)
        
        self.Show(1)        

    def OnStart(self, evt):
        self.SetStatusText("Running Tests")
        test()
        self.SetStatusText("Finished Running Tests")

    def OnClose(self, evt):
        self.Close(1)

if __name__ == "__main__":
    if not threaded():        
        app = TesterApp()
        app.MainLoop()
    else:
        test()

