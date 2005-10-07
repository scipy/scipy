#!/usr/bin/env python
#  $Id$
#  ---------------------------------------------------------------------
#
#  NAME:  tkgist.py
#
#  CHANGES:
#  12/06/02 dhm Originated.
#  01/21/03 llc gistC should be gist.
#  
#  ---------------------------------------------------------------------
"""
Make gistC module compatible with _tkinter module
(The work is done when this module is imported.)

Need to ensure two things:
(1) pyg_idler gets called by tcl/tk as an "after idle" method
    - this allows gist to complete deferred actions
(2) pyg_pending gets called by tcl/tk whenever events arrive for
    a gist window
    - this is unnecessary under win32 because the OS automatically
      sends window messages to the gist winproc
"""

__version__ = "$Id$"

import sys
import Tkinter
import _tkinter
import scipy.xplt.gist as gist

_root = None
_refresh = 1

# problem: how can you reset an "after idle" to make it into a true
#          idle event??

def do_pending(file, mask):
   global _root
   global _refresh
   print "python do_pending"
   if _root and _refresh:
      _refresh = 0
   gist.pyg_pending()

#def do_idle(*a):
def do_idle(dis, fd):
   global _root
   global _refresh
   print "python do_idle; dis=", dis, "fd=", fd
   gist.pyg_idler()
   _refresh = 1

   if _root:
      print _root
      if dis == 0:
         print dis
         _root.after_idle(do_idle)
#        _root.after(100, do_idle)

def do_connector(dis, fd):
    global _root
    print "python do_connector; dis=", dis, "fd=", fd
    if not _root:
        _root = Tkinter._default_root
        print "set _root", _root
    if fd >= 0:
        if dis == 0:
            _tkinter.createfilehandler(fd, _tkinter.READABLE, do_pending)
            print "file handle created"
        else:
            _tkinter.deletefilehandler(fd)
            print "file handle deleted"

def do_keyhandler(msg):
    print "python do_keyhandler"
    print msg
    # should evaluate msg as python command here

""" 
Initialization for Tkinter cooperating with Gist
"""

#  .. Make very sure gist has not stolen PyOS_InputHook
#     -- probably possible for this to be insufficient

gist.pyg_unhook()

#  .. can set this explicitly if Tkinter._default_root is not set
#     e.g.- if running IDLE (or anything else that uses PyShell):
#   import tkgist
#   tkgist._root = PyShell.root
#   from gist import *

print "_root, _refresh", _root, _refresh
gist.pyg_register(do_connector, do_keyhandler);
