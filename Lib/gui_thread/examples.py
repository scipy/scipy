import gui_thread
from wxPython.wx import *

def modal_sample():
    import gui_thread
    prxyDirDialog = gui_thread.register(wxDirDialog)
    q=prxyDirDialog(None)
    q.ShowModal()
    
class TrackingEvtHandler:
    """ Used occasionly for debugging which handlers are 
        connected
    """
    def Connect( self, *arguments, **namedarguments ):
        print 'arguments', arguments, namedarguments
        if not hasattr( self, "__ids" ):
            self.__ids = []
        self.__ids.append( arguments[0] )
        return apply( wxEvtHandler.Connect, (self,)+arguments,namedarguments )
    def ReleaseEventHandlers( self ):
        while self.__ids:
            self.Disconnect( self.__ids[-1] )
            del self.__ids[-1]

#class SimpleFrame(TrackingEvtHandler, wxFrame):
class SimpleFrame( wxFrame):
    def __init__(self, parent=None, ID=-1, title='no title', 
                 pos=wxDefaultPosition, size=wxSize(350,200)):
        wxFrame.__init__(self, parent, ID, title, pos, size)
        panel = wxPanel(self, -1)

        button = wxButton(panel, 1003, "Close Me")
        button.SetPosition(wxPoint(15, 15))
        self.button = button       
        EVT_BUTTON(self, 1003, self.OnCloseMe)
        EVT_CLOSE(self, self.OnCloseWindow)
        # Ah! The next line is ignored.  It would seem you must Disconnect
        # events before adding new ones.
        EVT_CLOSE(self, self.OnCloseWindow2)
        self.name = 'bubba'
        
    def Blue(self):
        self.button.SetBackgroundColour(wxBLUE)

    def Red(self):
        self.button.SetBackgroundColour(wxRED)

    def OnCloseMe(self, event):
        print 'hit'
        self.Close(true)

    def OnCloseWindow(self, event):
        print 'first called'
        self.Destroy()

    def OnCloseWindow2(self, event):
        print '2nd called'
        self.Destroy()

def non_modal_sample():
    prxySimpleFrame = gui_thread.register(SimpleFrame)
    win = prxySimpleFrame(None, -1, "This is a wxFrame", wxDefaultPosition, 
                      wxSize(350, 200))
    win.Show(true)
    return win