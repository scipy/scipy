""" Testing gui_thread.  Please run this under gui_thread and also
without it.  Right now the testing is extremely rudimentary."""

import unittest, gui_thread, weakref, time, sys

from wxPython.wx import *

def threaded():
    return gui_thread.main.running_in_second_thread

class TestPanel(wxPanel):
    def __init__(self, parent):
        wxPanel.__init__(self, parent, -1)
        ID = NewId()
        btn = wxButton(self, ID, "Hello")
        EVT_BUTTON(self, ID, self.OnButton)

##    def __del__(self):
##        print "TestPanel.__del__"                

    def OnButton(self, evt):
        print "TestPanel.OnButton"


class TestFrame(wxFrame):
    def __init__(self, parent):
##        print "TestFrame.__init__"
        wxFrame.__init__(self, parent, -1, "Hello Test")
        panel = TestPanel(self)
        EVT_CLOSE(self, self.OnClose)
        self.Show(1)

##    def __del__(self):
##        print "TestFrame.__del__"                

    def OnClose(self, evt):
        self.Destroy()


class TestClass:
    def __init__(self):
        self.a = 1

    def test(self):
        print self.a

def is_alive(obj):
    if obj() is None:
        return 0
    else:
        return 1


class test_gui_thread(unittest.TestCase):
    def check_wx_class(self):
        "Checking a wxFrame proxied class"
        for i in range(5):
            f = gui_thread.register(TestFrame)
            a = f(None)
            p = weakref.ref(a)
            a.Close(1)
            del a
            time.sleep(0.25) # sync threads
            # this checks for memory leaks
            self.assertEqual(is_alive(p), 0)
            
    def check_normal_class(self):
        "Checking non-wxWindows proxied class "
        f = gui_thread.register(TestClass)
        a = f()
        p = weakref.ref(a)
        # the reference count has to be 2.
        self.assertEqual(sys.getrefcount(a), 2)
        del a
        self.assertEqual(is_alive(p), 0)        
        
    def check_exception(self):
        "Checking exception handling"
        f = gui_thread.register(TestFrame)
        a = f(None)
        self.assertRaises(TypeError, a.Close, 1, 2, 3)
        a.Close()
        del a

def test_suite():
    suites = []
    suites.append(unittest.makeSuite(test_gui_thread, 'check_'))
    total_suite = unittest.TestSuite(suites)
    return total_suite

def test():
    all_tests = test_suite()
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(all_tests)


class NoThreadTestFrame(wxFrame):
    def __init__(self, parent):
        wxFrame.__init__(self, parent, -1, "Hello Test")
        test()
        self.Close(1)

if __name__ == "__main__":
    if not threaded():        
        app = wxPySimpleApp()
        frame = NoThreadTestFrame(None)
        app.MainLoop()
    else:
        test()

