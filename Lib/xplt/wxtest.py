#!/usr/bin/env python
#  $Id$
#  --------------------------------------------------------------------
#
#  NAME:     wxtest.py
#
#  PURPOSE:  Test interaction when combining wxPython and PyGist
#
#  NOTES:
#  Dave Grote (LBNL) was able to develop a GUI with wxWindows and PyGist
#  because wxWindows has an idle event.  Because in wxPython, the idle
#  event is only called once each time the window goes idle, Dave came
#  up with a workaround available in wxPython (to request that another
#  idle event be sent).  His solution is illustrated below.
#  
#  While we can get the GUI and the Gist window to each respond to 
#  events, the solution is not perfect.  If the GUI window covers the 
#  Gist window, the Gist window is not redrawn when exposed until the
#  mouse enters GUI window. 
#
#  CHANGES:
#  03/06/03 llc Originated.
#  --------------------------------------------------------------------

from wxPython.wx import *
import sys
import gist

ID_GIST = wxNewId()
ID_QUIT = wxNewId()

class MainFrame (wxFrame):
   def __init__ ( self, parent, id ):
      wxFrame.__init__ ( self, parent, id, "Idle Test" )
      wxButton ( self, ID_GIST, "Gist", pos=(0,0) )
      wxButton ( self, ID_QUIT, "Quit", pos=(0,50) )
      EVT_BUTTON ( self, ID_GIST, self.onGist )
      EVT_BUTTON ( self, ID_QUIT, self.onQuit )
      EVT_IDLE ( self, self.onIdle )
      self.isGistWindowOn = 0
   def onGist ( self, event ):
#     print "... Gist button pressed"
      self.isGistWindowOn = 1
      gist.plg([0,1]) 
   def onQuit ( self, event ):
#     print "... Quit button pressed"
      self.isGistWindowOn = 0
      sys.exit()
   def onIdle ( self, event ):
#     print "... Idle event"
      if self.isGistWindowOn:
         gist.pyg_pending() # Handle any new events
         gist.pyg_idler()   # Redraw window after changes from events
         event.RequestMore (1)

class MainApp ( wxApp ):
   def OnInit ( self ):
      self.frame = MainFrame ( NULL, -1 )
      self.frame.Show ( true )
      self.SetTopWindow ( self.frame )
      return true

if (__name__=="__main__"):
   app = MainApp (0)
   app.MainLoop()
