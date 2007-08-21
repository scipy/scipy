## Automatically adapted for scipy Oct 31, 2005 by

#!/usr/bin/env python
#  $Id$
#  -----------------------------------------------------------------
#
#  NAME:     tktest.py
#
#  PURPOSE:  Test interaction when combining Tkinter and PyGist
#
#  NOTES:
#  Tkinter and PyGist do not work together.
#  Problem:  Tkinter does not have a proper idle event.
#  Responds to Tk events and Gist mouse events, but:
#  - Gist window does not update properly
#  - Python command line input is no longer accepted
#
#  -----------------------------------------------------------------

import Tkinter
import sys
import tkgist
from gist import *
import gistdemolow

def DrawOval(Event):
    # Event.widget will be the main canvas:
    Event.widget.create_oval(Event.x-5,Event.y-5, Event.x+5,Event.y+5)
def DrawRectangle(Event):
    Event.widget.create_rectangle(Event.x-5,Event.y-5, Event.x+5,Event.y+5)
def MoveButton(Side):
    # The methods pack_forget() and grid_forget() unpack
    # a widget, but (unlike the destroy() method)
    # do not destroy it; it can be re-displayed later.
    QuitButton.pack_forget()
    QuitButton.pack(side=Side)

def plot():
    print ".. Gist button pressed"
    plg([0,1])

def main():
    root=Tkinter.Tk()
    MainCanvas=Tkinter.Canvas(root)
    MainCanvas.bind("<Button-1>",DrawOval)
    MainCanvas.bind("<Button-3>",DrawRectangle)
    MainCanvas.pack(fill=Tkinter.BOTH,expand=Tkinter.YES)

    gistButton=Tkinter.Button(MainCanvas,text="Gist",
     command=plot)
    quitButton=Tkinter.Button(MainCanvas,text="Quit",
     command=sys.exit)
    gistButton.pack(side=Tkinter.TOP)
    quitButton.pack(side=Tkinter.BOTTOM)

    root.bind("<Up>",lambda e:MoveButton(Tkinter.TOP))
    root.bind("<Down>",lambda e:MoveButton(Tkinter.BOTTOM))
    root.bind("<Left>",lambda e:MoveButton(Tkinter.LEFT))
    root.bind("<Right>",lambda e:MoveButton(Tkinter.RIGHT))
    root.geometry("300x300") # Set minimum window size

    root.mainloop()

if (__name__=="__main__"):
    main()
