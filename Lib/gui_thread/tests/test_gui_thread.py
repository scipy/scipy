""" Testing gui_thread.  Please run this under gui_thread and also
without it.  Right now the testing is extremely rudimentary.  The
proxied return types are tested so are a few simple classes.

To test the code you should run the tests both without and with using
gui_thread in another thread.  Like so:

 $ python
 >>> import test_gui_thread
 >>> test_gui_thread.test()
 # or
 $ python test_gui_thread.py

And then:
 $ python
 >>> import gui_thread
 >>> <Importing wxPython...  done.>
 >>> import test_gui_thread
 
 >>> test_gui_thread.test()

"""

from wxPython.wx import *

import unittest, gui_thread, weakref, time, sys, math

def threaded():
    return gui_thread.main.running_in_second_thread

def is_alive(obj):
    if obj() is None:
        return 0
    else:
        return 1

def Yield():
    if not threaded():
        # forces the event handlers to finish their work.
        # this also forces deletion of windows.
        wxYield() 
    else:
        time.sleep(0.05) # sync threads
    

class TestPanel(wxPanel):
    def __init__(self, parent):
        wxPanel.__init__(self, parent, -1)
        self.attr = None
        ID = NewId()
        btn = wxButton(self, ID, "Hello")
        EVT_BUTTON(self, ID, self.OnButton)

    def OnButton(self, evt):
        print "TestPanel.OnButton"


class TestFrame(wxFrame):
    def __init__(self, parent):
        wxFrame.__init__(self, parent, -1, "Hello Test")
        self.panel = TestPanel(self)
        self.attr = [0, 0]
        EVT_CLOSE(self, self.OnClose)
        self.Show(1)

    def OnClose(self, evt):
        self.Destroy()


class TestClass:
    def __init__(self):
        self.a = 1

    def test(self):
        return self.a

class _TestHang:
    def get_obj(self):
        val = TestHang()
        return val

TestHang = gui_thread.register(_TestHang)

class _A:
    pass

A = gui_thread.register(_A)

class _TestDeref:
    def __init__(self):
        self.a = A()
        
    def test_proxy(self, arg1):
        "Returns 0 if the test fails."
        if (gui_thread.main.is_proxy(arg1) == 1) or \
           (gui_thread.main.is_proxy_attr(arg1) == 1):
            return 0
        else:
            return 1

    def test_proxy_attr(self):
        a = TestDeref()
        return self.test_proxy(a.a)

    def test_kw(self, arg1=None):
        return self.test_proxy(arg1)

TestDeref = gui_thread.register(_TestDeref)

dummy_instance = TestClass()

class TestProxyAttr:
    def __init__(self):
        self.lst = [1,2]
        self.dict = {'a':'a', 'b':'b'}

    def test(self):
        return 1

    def get_int(self):
        return 1

    def get_float(self):
        return math.pi

    def get_complex(self):
        return 1.0 + 5j

    def get_string(self):
        return "This is a string"

    def get_list(self):
        return [[0,1,2, 'foo', dummy_instance], 1,2, "foo", dummy_instance]

    def get_tuple(self):
        return ((1,2,3), 1, 2, 'foo',  dummy_instance)

    def get_dict(self):
        dict = {'int': 1, 'string': "This is a string",
                'list': [1,2,'foo'], 'tuple':(1,2, 'foo'),
                'obj': dummy_instance,
                'dict': {'int':1, 'string': 'This is a string'}}
        return dict

    def get_callable(self):
        return self.test

class test_gui_thread(unittest.TestCase):
    def check_wx_class(self):
        "Checking a wxFrame proxied class"
        for i in range(5):
            f = gui_thread.register(TestFrame)
            a = f(None)
            p = weakref.ref(a)
            a.Close(1)
            del a
            Yield()
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
        Yield()
        # this checks for memory leaks
        self.assertEqual(is_alive(p), 0)

    def check_proxy_proxy(self):
        "Checking proxied function calling another proxied function"
        a = TestHang()
        a.get_obj()

    def check_dereference(self):
        "Checking if dereferencing works."
        a = TestDeref()
        self.assertEqual(a.test_proxy(a), 1)
        self.assertEqual(a.test_kw(arg1=a), 1)
        self.assertEqual(a.test_proxy_attr(), 1)

    def check_setattr(self):
        "Checking if __setattr__ works."
        f = gui_thread.register(TestFrame)
        a = f(None)
        a.attr = [1, "string"]
        self.assertEqual(a.attr[0], 1)
        self.assertEqual(a.attr[1], "string")

        # checking the proxy_attr's __setattr__
        a.panel.attr = [1, "string"]
        self.assertEqual(a.panel.attr[0], 1)
        self.assertEqual(a.panel.attr[1], "string")
        a.Close()


class test_proxy_attribute(unittest.TestCase):
    
    """ This class tests the different types of smart returning used
    by the proxy classes.  The aim is to make the proxy behave very
    much like the un-proxied classes so checking the return types is
    important especially because lists, tuples and dictionaries have
    to be properly proxied."""
    
    def check_int(self):
        "Checking proxying of function returning int"
        f = gui_thread.register(TestProxyAttr)
        a = f()
        self.assertEqual(a.get_int(), 1)

    def check_float(self):
        "Checking proxying of function returning float"
        f = gui_thread.register(TestProxyAttr)
        a = f()
        self.assertEqual(a.get_float(), math.pi)

    def check_complex(self):
        "Checking proxying of function returning complex"
        f = gui_thread.register(TestProxyAttr)
        a = f()
        self.assertEqual(a.get_complex(), 1.0 + 5j)
            
    def check_string(self):
        "Checking proxying of function returning string"
        f = gui_thread.register(TestProxyAttr)
        a = f()
        self.assertEqual(a.get_string(), "This is a string")

    def check_list(self):
        "Checking proxying of function returning list"
        f = gui_thread.register(TestProxyAttr)
        a = f()
        no_proxy = TestProxyAttr()
        l = a.get_list()
        l1 = no_proxy.get_list()
        self.assertEqual(l1, l)

    def check_tuple(self):
        "Checking proxying of function returning tuple"
        f = gui_thread.register(TestProxyAttr)
        a = f()
        no_proxy = TestProxyAttr()
        t = a.get_tuple()
        t1 = no_proxy.get_tuple()
        for i in range(len(t1)):
            if type(t1[i]) is types.TupleType:
                self.assertEqual(list(t1[i]), list(t[i]))
            elif isinstance(t1[i], TestClass):
                self.assertEqual(t[i].test(), t1[i].test())
            else:
                self.assertEqual(t[i], t1[i])
    
    def check_dict(self):
        "Checking proxying of function returning dictionary"
        f = gui_thread.register(TestProxyAttr)
        a = f()
        no_proxy = TestProxyAttr()
        d = a.get_dict()
        d1 = no_proxy.get_dict()
        self.assertEqual(d.keys(), d1.keys())
        for i in d1.keys():
            if i == 'tuple':
                self.assertEqual(list(d[i]), list(d1[i]))
            elif i == 'obj':
                self.assertEqual(d[i].test(), d1[i].test())
            else:
                self.assertEqual(d[i], d1[i])

    def check_callable(self):
        "Checking proxying of function returning callable"
        f = gui_thread.register(TestProxyAttr)
        a = f()
        no_proxy = TestProxyAttr()
        self.assertEqual(a.get_callable()(), no_proxy.get_callable()())

    def check_attr_list(self):
        "Checking if changing list attributes work"
        f = gui_thread.register(TestProxyAttr)
        a = f()
        a.lst.append(3)
        self.assertEqual(a.lst[-1], 3)
        
    def check_attr_dict(self):
        "Checking if changing dictionary attributes work"
        f = gui_thread.register(TestProxyAttr)
        a = f()
        a.dict['c'] = 'c'
        self.assertEqual(a.dict['c'], 'c')

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
