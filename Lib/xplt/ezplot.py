# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

# The following more or less permanent defaults can be overruled
# by assignment, provided that this file is execfile'd rather than import'ed.
defgridx = "off"
defgridy = "off"
defthick = 1.0
defmark = " "
defstyle = "solid"
deflabels = 1
deflabel = None
defscale = "linlin"
defcolor = "fg"
deftop = ""
defbot = ""
defleft = ""
defright = ""
defvsc = .5
deflev = 8
defcgm = "yes"
defps = "no"

_null_object_list_ = [ [], [], [], [], [], [], [], []]

def _set_frame_defaults_ ( ) :
   global _titles_, _object_list_, _grid_type_, _axis_labels_, _axis_limits_, \
          _axis_scales_, _gridx_, _gridy_, _thick_, _mark_, _style_, \
          _labels_enabled_, _axis_scales_, _color_, _label_list_, _text_, \
          _text_pos_, _text_size_, _text_color_, _tosys_, _cmin_, \
          _cmax_, _vsc_, _lev_, _krange_, _lrange_, _kstyle_, _lstyle_, \
          _region_, _default_label_, _color_bar_
   _region_ = "all"
   _lev_ = deflev
   _krange_ = None
   _lrange_ = None
   _kstyle_ = "solid"
   _lstyle_ = "solid"
   _vsc_ = defvsc
   _titles_ = [defbot, deftop,
            defleft, defright]  # bottom, top, left, right titles
   _object_list_ = [ [], [], [], [], [], [], [], []] # (ensure copy!!)
      # list of objects to plot (one for each possible device)
   _label_list_ = []
   _text_ = [" "]
   _text_pos_ = [ [0., 0.]]
   _text_size_ = [0]
   _text_color_ = ["bg"]
   _tosys_ = [1]
   _gridx_ = defgridx
   _gridy_ = defgridy
   _grid_type_ = "axes"         # grid type for current graph
   _axis_labels_ = ["X axis",
                    "Y axis",
                    "YR axis"]  # axis labels
   _axis_limits_ = [["d","d"],
                    ["d","d"],
                    ["d","d"]]  # axis limits; "d" implies default.
   _axis_scales_ = ["lin",
                    "lin",
                    "lin"]      # axis scales
   _thick_ = defthick
   _mark_ = defmark
   _style_ = defstyle
   _labels_enabled_ = deflabels
   _default_label_ = " "
   _axis_scales_ = defscale
   _color_ = defcolor
   _cmin_ = _cmax_ = 0.0
   _color_bar_ = 0

def _set_initial_defaults_ ( ) :
   # Set defaults
   global _ezcshow_, _ezcreset_, _window_, _cgm_, _ps_, _displays_, \
          _ps_file_, _cgm_plotter_, _ps_plotter_, _win_plotters_, \
          _current_graph_, _cgm_file_, _graphics_, _win_on_, _zt_, \
          _rt_, _ireg_, _pvar_, _cindex_, _ut_, _vt_
   _ezcshow_ = "true"           # no plot until nf
   _ezcreset_ = "true"          # return to initial values after nf
   _window_ = "no"              # no window
   _cgm_ = defcgm               # create cgm file
   _ps_ = defps                 # create postscript file
   _displays_ = [" ", " ", " ",
                 " ", " ", " ",
                 " ", " "]      # hosts where display occurs
   _graphics_ = ["", "", "",
                 "", "", "",
                 "", ""]        # type of graphics to use
   _ps_file_ = "Aa00.ps"        # name of postscipt file for display
   _cgm_file_ = "Aa00.cgm"      # name of cgm file for display
   _cgm_plotter_ = -1           # plotter which writes to cgm file
                                # (subscript into _win_plotters_)
   _ps_plotter_ = -1            # plotter which writes to postscript file
                                # (subscript into _win_plotters_)
   _win_plotters_ = [None,
                     None,
                     None,
                     None,
                     None,
                     None,
                     None,
                     None]      # plotters which plot to windows (max 8)
                                # _cgm_plotter_ and _ps_plotter_ are
                                # on this list, if they exist.
   _win_on_ = [0, 0, 0, 0,
               0, 0, 0, 0]      # windows actually activated
   _current_graph_ = [None,
                      None,
                      None,
                      None,
                      None,
                      None,
                      None,
                      None]      # current graph2d object
   _set_frame_defaults_ ( )
   _zt_ = None
   _rt_ = None
   _ireg_ = None
   _pvar_ = None
   _cindex_ = None
   _ut_ = None
   _vt_ = None

_set_initial_defaults_ ( )

GraphicsError = "GraphicsError"
import os
try:
   graphics = os.environ["PYGRAPH"]
except KeyError:
   graphics = "Gist"

from gist import * # to get mfit function


if graphics [0:3] == "Nar" :
   import NarPlotter
   GistPlotter = None
   Graphics = NarPlotter
   print "Opening Narcisse graphics."
elif graphics == "Gist" :
   import GistPlotter
   NarPlotter = None
   Graphics = GistPlotter
   print "Opening Gist graphics."
else :
   raise GraphicsError , \
      graphics + " is an unknown graphics package. Check PYGRAPH " + \
         "environment variable."

from region import *
from quadmesh import *
from curve import *
from graph2d import *
from cellarray import *

import shapetest
from Numeric import *
from scipy_base.fastumath import *

_ezdict_ = {'t': 'true' , 'T': 'true', 'y': 'true', 'Y': 'true',
            'f': 'false', 'F': 'false', 'n': 'false', 'N': 'false'}

def ezcshow (val) :
   global _ezcshow_, _ezdict_
   EzcShowErr = "EzcShowErr"
   if not _ezdict_.has_key (val [0][0]) :
      raise EzcShowErr , \
         "Illegal argument! must be 'true' or 'false'."
   _ezcshow_ = _ezdict_ [val [0][0]]

def ezcreset (val) :
   global _ezcreset_, _ezdict_
   EzcResetErr = "EzcResetErr"
   if not _ezdict_.has_key (val [0][0]) :
      raise EzcResetErr , \
         "Illegal argument! must be 'true' or 'false'."
   _ezcreset_ = _ezdict_ [val [0][0]]

from gistC import bytscl

def ezcfmc (val) :
   """returns an array of indices into the current palette."""
   return bytscl (val)

_WinSpecError_ = "WinSpecError"

def _get_free_win_number_ ( ) :
   # Gets the lowest number not associated with an open device, else -1.
   global _win_on_
   try :
      return _win_on_.index (0)
   except:
      return -1

def _get_first_open_win_number () :
   # Gets the number of the lowest numbered open window, if any, else -1..
   global _win_on_, _cgm_plotter_, _ps_plotter_
   for i in range (8) :
      if _win_on_ [i] != 0 and i != _cgm_plotter_ and i != _ps_plotter_:
         return i
   return -1

def win (val , *n, **keywords) :
   """win (str [, n] [, display = <string1>] [, graphics = <string2>])
         str = "on" : reserves a graphics window for display. n, if
            present, must be between 0 and 7 and specifies which window.
            if n is absent, it will default to the lowest-numbered window
            available. display, if present, specifies via <string1>
            where the window will be displayed, e. g., "icf.llnl.gov:0.0".
            The default is $DISPLAY. graphics, if present, specifies
            via <string2> which type of graphics is desired, "Nar"
            for Narcisse, and "Gist" for Gist being the currently
            allowed values. The default is $PYGRAPH if set, otherwise
            "Gist".
         str = "close" : closes a graphics window; wndow n if
            specified, otherwise window 0.
   """
   global _win_on_, _displays_, _ps_plotter_, _cgm_plotter_, _graphics_, \
      _win_plotters_, _object_list_
   if len (n) > 0 :
      _window_number_ = n [0]
   else :
      if keywords.has_key ("window") :
         _window_number_ = keywords ["window"]
      else :
         _window_number_ = "all"
   if val == "on" or val == "ON" or val == "On" \
      or val == "open" or val == "OPEN" or val == "Open" :
      _window_number_ = _check_win_ (_window_number_)
      if _window_number_ == []:
         _window_number_ = [_get_free_win_number_ ( )]
         if _window_number_ [0] < 0 :
            raise _WinSpecError_, "win ('open'): all windows are in use."
      for i in _window_number_:
         if i < 0 or i > 7 :
            raise _WinSpecError_, "Window number (" + `i` + \
               ") is not between 0 and 7."
         _win_on_ [i] = 1
         if keywords.has_key ("display") :
            _displays_ [i] = keywords ["display"]
         if keywords.has_key ("graphics") :
            _graphics_ [i] = keywords ["graphics"]
   if val == "off" or val == "OFF" or val == "Off" or \
      val == "close" or val == "CLOSE" or val == "Close" :
      _window_number_ = _check_win_ (_window_number_, minmax = 1)
      if _window_number_ == []:
         _window_number_ = [0]
      for i in _window_number_:
         if i < 0 or i > 7 :
            raise _WinSpecError_, "win ('close'): window number (" + `i` + \
               ") is not between 0 and 7."
         _displays_ [i] = " "
         _graphics_ [i] = ""
         _win_on_ [i] = 0
         _object_list_ [i] = []
         if _win_plotters_ [i] is not None :
            _win_plotters_ [i].close ()
            _win_plotters_ [i] = None
         if i == _ps_plotter_ :
            _ps_plotter_ = -1
         elif i == _cgm_plotter_ :
            _cgm_plotter_ = -1

tv = win

_CgmError_ = "CgmError"

def new_plot_file ( filename ) :
   """new_plot_file ( ) will eventually return the name of the next
   available file in the sequence "Aa00.ps/cgm", Ab00.ps/cgm", etc.
   """
   if not os.path.isfile ("./" + filename) :
      return filename
   n = string.find (string.letters, filename [1])
   if n == -1 or n == len (string.letters) - 1 :
      raise _PsError_ , "No more ps filenames available."
   letters = string.letters [n + 1:]
   OK = 0
   for l in string.letters :
      tmpfile = filename[0:1] + l + filename[2:]
      if not os.path.isfile ("./" + tmpfile) :
         OK = 1
         break
   if OK == 1 :
      return tmpfile
   else :
      raise _PsError_ , "No more ps or cgm filenames available."

def _min_ ( x ) :
   if is_scalar (x) :
      return x
   else :
      return min (ravel (x))

def _max_ ( x ) :
   if is_scalar (x) :
      return x
   else :
      return max (ravel (x))

def _get_objects_ () :
   global _object_list_, _null_object_list_
   _this_object_list_ = []
   if _object_list_ != _null_object_list_ :
      for i in range (8):
         if _object_list_ [i] != [] :
            _this_object_list_.append (i)
   return _this_object_list_

def _get_open_win_list ():
   global _win_on_
   openl = []
   for i in range (len (_win_on_)) :
      if _win_on_ [i]:
         openl.append (i)
   return openl


_BadWindow = "BadWindow"

def _check_win_ (window, minmax = 0) :
   """_check_win_ (window) checks the integrity of the window specification,
   and if OK, returns a list of windows. If minmax = 0, treats "min" and "max"
   as errors, otherwise allows them."""
   global _win_on_
   if window == "all" :
      return _get_open_win_list ()
   elif minmax == 1 and window == "min" :
      for i in range (len (_win_on_)) :
         if _win_on_ [i]:
            return [i]
   elif minmax == 1 and window == "max" :
      for i in range (len (_win_on_)) :
         if _win_on_ [7 - i]:
            return [7 - i]
   elif type (window) == IntType and 0 <= window <= 7:
      return [window]
   elif type (window) == ListType:
      err = 0
      for i in range (len (window)):
         if type (window [i]) == IntType and 0 <= window [i] <= 7:
            continue
         else:
            err = 1
            break
      if err == 0:
         return window
   raise _BadWindow, "Incorrect window specification: " + `window` + "."

def _set_axis_limits_ (window = "all") :
   """_set_axis_limits_ ( ) makes sure legal values are used for the defaults."""
   global _axis_limits_, _object_list_, _null_object_list_
   ## Someday this will be improved to have a list of object lists, one for each device.
   _object_numbers_ = _check_win_ (window, minmax = 1)
   if _object_numbers_ == []:
      _axis_limits_ = [ [0., 0.], [0., 0.], [0., 0.]]
      return
   min_ = max_ = None
   for j in range (len (_object_numbers_)):
      if _object_list_ [j] == [] :
         continue
      for i in range (3) :
         if _axis_limits_ [i] == ["d", "d"] :
            _axis_limits_ [i] = [0., 0.]
            continue
         if _axis_limits_ [i] [0] == "d" :
            if i == 0 :
               for k in range (len (_object_list_ [j])) :
                  new_min = _min_ (_object_list_ [j] [k].x)
                  if min_ is None or new_min < min_ :
                     min_ = new_min
            elif i == 1 :
               for k in range (len (_object_list_ [j])) :
                  new_min = _min_ (_object_list_ [j] [k].y)
                  if min_ is None or new_min < min_ :
                     min_ = new_min
            elif i == 2 :
               min_ = 0.
            _axis_limits_ [i] [0] = min_
         if _axis_limits_ [i] [1] == "d" :
            if i == 0 :
               for k in range (len (_object_list_ [j])) :
                  new_max = _max_ (_object_list_ [j] [k].x)
                  if max_ is None or new_max > max_ :
                     max_ = new_max
            elif i == 1 :
               for k in range (len (_object_list_ [j])) :
                  new_max = _max_ (_object_list_ [j] [k].y)
                  if max_ is None or new_max > max_ :
                     max_ = new_max
            elif i == 2 :
               max_ = 0.
            _axis_limits_ [i] [1] = max_

_CgmError_ = "CgmError"

def _get_win_ (window, error):
   """_get_win_ (window) returns the number of a window which has a
   non-empty _object_list_. Has an exception if there is none--nothing to plot."""
   global _object_list_
   win_num = -1
   if window == "min":
      # Find the earliest window with a nonempty list (if any)
      for i in range (8):
         if _object_list_ [i] != []:
            win_num = i
            break

   elif window == "max":
      # Find the last window with a nonempty list (if any)
      for i in range (8):
         if _object_list_ [7 - i] != []:
            win_num = 7 - i
            break

   elif 0 <= window <= 7:
      # Check for a nonempty object list:
      if _object_list_ [window] != [] :
         win_num = window

   else:
      # Find the earliest window with a nonempty list (if any)
      for i in range (8):
         if _object_list_ [i] != []:
            win_num = i
            break
      
   if win_num == -1:
      raise error, "There seems to be nothing there!"
   return win_num
   
def cgm (val, *n, **kw) :
   global _cgm_, _cgm_plotter_, _current_graph_, _cgm_file_, _win_on_, \
          _win_plotters, _cmin_, _cmax_, _titles_, _axis_limits_, \
          _axis_labels_, _axis_scales_, _text_, _text_pos_, \
          _text_size_, _text_color_, _tosys_, _grid_type_, _object_list_
   if kw.has_key ("new_frame") :
      new_frame = kw ["new_frame"]
   elif _ezcshow_ == "false" :
      new_frame = "no"
   else :
      new_frame = "yes"
   if val == "on" or val == "ON" or val == "On" \
      or val == "open" or val == "OPEN" or val == "Open" :
      _cgm_ = "yes"
      if len (n) == 0 :
         if 7 >= _cgm_plotter_ >= 0 :
            if _win_on_ [_cgm_plotter_] :
               return # a window is already assigned to cgm and is open
            else :
               _win_on_ [_cgm_plotter_] = 1
               return # a window is already assigned to cgm
         else :
            _window_number_ = _get_free_win_number_ ()
            if _window_number_ < 0 :
               raise _CgmError_, "No window available for cgm use."
      else :
         _window_number_ = n [0]
         if _window_number_ < 0 or _window_number_ > 7 :
            raise _CgmError_, "Window number (" + `_window_number_` + \
               ") is not between 0 and 7."
         if 0 <= _cgm_plotter_ <= 7 :
            if _cgm_plotter_ == _window_number_ :
               return # this window is already assigned to cgm
            else :
               # another cgm window was assigned. unassign it.
               _win_on_ [_cgm_plotter_] = 0
               if _win_plotters_ [_cgm_plotter_] is not None :
                  _win_plotters_ [_cgm_plotter_].close ()
                  _win_plotters_ [_cgm_plotter_] = None
      _cgm_plotter_ = _window_number_
      _win_on_ [_cgm_plotter_] = 1
      return
   if val == "close" or val == "CLOSE" or val == "Close" :
      _cgm_ = "no"
      if 0 <= _cgm_plotter_ <= 7 :
         print "Closing cgm file " + _cgm_file_ + "."
         _win_on_ [_cgm_plotter_] = 0
         if _win_plotters_ [_cgm_plotter_] is not None :
            _win_plotters_ [_cgm_plotter_].close ()
            _win_plotters_ [_cgm_plotter_] = None
         _cgm_plotter_ = -1
      return
   if val == "off" or val == "OFF" or val == "Off" :
      # Turn off transmissions to window but don't close
      _cgm_ = "no"
      return
   if val == "send" or val == "SEND" or val == "Send" or val == "plot" :
      # Figure out which object list is to be sent
      if _cgm_plotter_ < 0 or _object_list_ [_cgm_plotter_] == 0 and \
         not kw.has_key ("window"):
         kw ["window"] = "min"
      if kw.has_key ("window") :
         # User has specified a window whose object list is to be plotted
         # or else the cgfm window doesn't have a plotting list of its own
         # -- guaranteed to return a valid window number; exception otherwise.
         win_no = _get_win_ (kw ["window"], _CgmError_)
      # time to send a graph to a cgm file
      if _cgm_plotter_ < 0 :
         _cgm_plotter_ = _get_free_win_number_ ()
         if _cgm_plotter_ < 0 :
            raise _CgmError_, "No window available for cgm use."
         _win_on_ [_cgm_plotter_] = 1
      if _win_plotters_ [_cgm_plotter_] is None :
         _cgm_file_ = new_plot_file (_cgm_file_)
         print "Opening plot file " + _cgm_file_ + "."
         _win_plotters_ [_cgm_plotter_] = \
            GistPlotter.Plotter ( "none" , n = _cgm_plotter_, hcp = _cgm_file_)
      _win_plotters_ [_cgm_plotter_].set_bytscl (_cmin_, _cmax_)
      _set_axis_limits_ ( )
      if _cgm_plotter_ != win_no :
         _object_list_ [_cgm_plotter_] = _object_list_ [win_no]
      if _current_graph_ [_cgm_plotter_] is None :
         _current_graph_ [_cgm_plotter_] = Graph2d (_object_list_ [win_no], \
            plotter = _win_plotters_ [_cgm_plotter_], \
            titles = _titles_, axis_limits = _axis_limits_, axis_labels = \
            _axis_labels_, axis_scales = _axis_scales_,  \
            text = _text_, text_pos = _text_pos_, text_size = \
            _text_size_, text_color = _text_color_, tosys = _tosys_, \
            grid_type = _grid_type_)
      elif val == "plot" :
          _current_graph_ [_cgm_plotter_].set (plotter = _win_plotters_ [_cgm_plotter_],
              axis_limits = _axis_limits_)
      else :
         _current_graph_ [_cgm_plotter_].new (_object_list_ [win_no], \
            plotter = _win_plotters_ [_cgm_plotter_], \
            text = _text_, text_pos = _text_pos_, text_size = \
            _text_size_, text_color = _text_color_, tosys = _tosys_, \
            titles = _titles_, axis_limits = _axis_limits_, axis_labels = \
            _axis_labels_, axis_scales = _axis_scales_, grid_type = _grid_type_)
      if new_frame == "no" :
         _current_graph_ [_cgm_plotter_].plot ( )
         _win_plotters_ [_cgm_plotter_].plot_text ( )
      if val != "plot" : # Got to turn off the file.
         _cgm_ = "no"

import string

_PsError_ = "PsError"

def ps (val, *n, **kw) :
   global _ps_, _ps_plotter_, _current_graph_, _ps_file_, _win_on_, \
          _win_plotters, _cmin_, _cmax_, _titles_, _axis_limits_, \
          _axis_labels_, _axis_scales_, _text_, _text_pos_, \
          _text_size_, _text_color_, _tosys_, _grid_type_, _object_list_
   if kw.has_key ("new_frame") :
      new_frame = kw ["new_frame"]
   elif _ezcshow_ == "false" :
      new_frame = "no"
   if val == "on" or val == "ON" or val == "On" \
      or val == "open" or val == "OPEN" or val == "Open" :
      _ps_ = "yes"
      if len (n) == 0 :
         if 7 >= _ps_plotter_ >= 0 :
            if _win_on_ [_ps_plotter_] :
               return # a window is already assigned to ps
            else :
               _win_on_ [_ps_plotter_] = 1
               return
         else :
            _window_number_ = _get_free_win_number_ ()
            if _window_number_ < 0 :
               raise _PsError_, "No window available for ps use."
      else :
         _window_number_ = n [0]
         if _window_number_ < 0 or _window_number_ > 7 :
            raise _PsError_, "Window number (" + `_window_number_` + \
               ") is not between 0 and 7."
         if 0 <= _ps_plotter_ <= 7 :
            if _ps_plotter_ == _window_number_ :
               return # this window is already assigned to ps
            else :
               # another ps window was assigned. unassign it.
               _win_on_ [_ps_plotter_] = 0
               if _win_plotters_ [_ps_plotter_] is not None :
                  _win_plotters_ [_ps_plotter_].close ()
                  _win_plotters_ [_ps_plotter_] = None
      _ps_plotter_ = _window_number_
      _win_on_ [_ps_plotter_] = 1
      return
   if val == "off" or val == "OFF" or val == "Off" :
      _ps_ = "no"
      return
   if val == "close" or val == "CLOSE" or val == "Close" :
      _ps_ = "no"
      if 0 <= _ps_plotter_ <= 7 :
         _win_on_ [_ps_plotter_] = 0
         if _win_plotters_ [_ps_plotter_] is not None :
            _win_plotters_ [_ps_plotter_].close ()
            _win_plotters_ [_ps_plotter_] = None
         _ps_plotter_ = -1
         print "Closing plot file " + _ps_file_ + "."
      return
   if val == "send" or val == "SEND" or val == "Send" or val == "plot" :
      # Figure out which object list is to be sent
      if _ps_plotter_ < 0 or _object_list_ [_ps_plotter_] == 0 and \
         not kw.has_key ("window"):
         kw ["window"] = "min"
      if kw.has_key ("window") :
         # User has specified a window whose object list is to be plotted
         # or else the cgfm window doesn't have a plotting list of its own
         # -- guaranteed to return a valid window number; exception otherwise.
         win_no = _get_win_ (kw ["window"], _PsError_)
      # time to send a graph to a ps file
      if _ps_plotter_ < 0 :
         _ps_plotter_ = _get_free_win_number_ ()
         if _ps_plotter_ < 0 :
            raise _PsError_, "No window available for ps use."
         _win_on_ [_ps_plotter_] = 1
      if _win_plotters_ [_ps_plotter_] is None :
         _ps_file_ = new_plot_file (_ps_file_)
         print "Opening plot file " + _ps_file_ + "."
         _win_plotters_ [_ps_plotter_] = \
            GistPlotter.Plotter ( "none" , n = _ps_plotter_, hcp = _ps_file_)
      _win_plotters_ [_ps_plotter_].set_bytscl (_cmin_, _cmax_)
      _set_axis_limits_ ( )
      if _ps_plotter_ != win_no :
         _object_list_ [_ps_plotter_] = _object_list_ [win_no]
      if _current_graph_ [_ps_plotter_] is None :
         _current_graph_ [_ps_plotter_] = Graph2d (_object_list_ [win_no], \
            plotter = _win_plotters_ [_ps_plotter_], \
            titles = _titles_, axis_limits = _axis_limits_, axis_labels = \
            _axis_labels_, axis_scales = _axis_scales_,  \
            text = _text_, text_pos = _text_pos_, text_size = \
            _text_size_, text_color = _text_color_, tosys = _tosys_, \
            grid_type = _grid_type_)
      elif val == "plot" :
          _current_graph_ [_ps_plotter_].set (plotter = _win_plotters_ [_ps_plotter_],
              axis_limits = _axis_limits_)
      else :
         _current_graph_ [_ps_plotter_].new (_object_list_ [win_no], \
            plotter = _win_plotters_ [_ps_plotter_], \
            text = _text_, text_pos = _text_pos_, text_size = \
            _text_size_, text_color = _text_color_, tosys = _tosys_, \
            titles = _titles_, axis_limits = _axis_limits_, axis_labels = \
            _axis_labels_, axis_scales = _axis_scales_, grid_type = _grid_type_)
      if new_frame == "no" :
         _current_graph_ [_ps_plotter_].plot ( )
         _win_plotters_ [_ps_plotter_].plot_text ( )
      if val != "plot" : # Got to turn off the file.
         _ps_ = "no"


_NfError_ = "NfError"
_UndoError_ = "UndoError"

def undo (item = None, window = "min") :
   """undo (item, window) removes the item-th object from window's
   object list, if it can. item counts starting at 1."""
   global _object_list_
   win_no = _get_win_ (window, _UndoError_)
   lo = len (_object_list_ [win_no])
   if item is None :
      n = lo
   if n > lo or n < 1:
      raise _UndoError_, "There aren't " + `n` + " objects in the list."
   _object_list_ [win_no] [n - 1:n] = []
   if _object_list_ [win_no] != [] and _ezcshow_ == "true" :
      sf (win_no)
   
def nf ( new_frame = "yes" , window = "all" ) :
   global _win_plotters_, _current_graph_, _object_list_, _label_list_, \
          _cmin_, _cmax_, _cgm_, _ps_, _labels_enabled_, _color_bar_, \
          _win_on_, _cgm_plotter_, _ps_plotter_, _titles_, _axis_limits_, \
          _axis_labels_, _axis_scales_, _text_, _text_pos_, _text_size_, \
          _text_color_, _tosys_, _grid_type_, _object_list_, _graphics_, \
          _current_graph_, _labels_temporarily_enabled_
   enablen = zeros (8, Int)
   enable_cgm = 0
   enable_ps = 0
   if _object_list_ == _null_object_list_ :
      # Nothing to plot
      return
   if _cgm_ == "yes" and _cgm_plotter_ < 0 :
      # Open default cgm plotter
      print "Opening default cgm plotter."
      cgm ("on")
   if window == "cgm" :
      enable_cgm = 1
   elif window == "ps" :
      enable_ps = 1
   else :
      windows = _check_win_ (window, minmax = 1)
      for i in windows:
         enablen [i] = 1
   _set_axis_limits_ ( )
   for i in range (len (_win_plotters_)) :
      for j in range (len (_object_list_ [i])) :
         if _labels_temporarily_enabled_ and len (_label_list_) > j :
            _object_list_ [i] [j].set (label = _label_list_ [j])
            _object_list_ [i] [j].set (marks = 1)
         else :
            _object_list_ [i] [j].set (label = " ")
            _object_list_ [i] [j].set (marks = 0)
#           if _object_list_ [i] [j].line_type != "none" :
#              _object_list_ [i] [j].set (marks = 0)
      if i == _cgm_plotter_ and (enable_cgm or enablen [i]) and _cgm_ :
         cgm ("plot", new_frame = new_frame, window = _cgm_plotter_)
      elif i == _ps_plotter_ and (enable_ps or enablen [i]) and _ps_ :
         ps ("plot", new_frame = new_frame, window = _ps_plotter_)
      elif enablen [i] :
         if _current_graph_ [i] is None :
            _current_graph_ [i] = Graph2d (_object_list_ [i], titles = _titles_,
               text = _text_, text_pos = _text_pos_, text_size = \
               _text_size_, text_color = _text_color_, tosys = _tosys_, \
               axis_limits = _axis_limits_, axis_labels =_axis_labels_,
               axis_scales = _axis_scales_, grid_type = _grid_type_,
               color_bar = _color_bar_)
         else :
            _current_graph_ [i].new (_object_list_ [i], titles = _titles_,
               text = _text_, text_pos = _text_pos_, text_size = \
               _text_size_, text_color = _text_color_, tosys = _tosys_, \
               axis_limits = _axis_limits_, axis_labels =_axis_labels_,
               axis_scales = _axis_scales_, grid_type = _grid_type_,
               color_bar = _color_bar_)
         if _win_on_ [i] == 1 :
            if _win_plotters_ [i] is None :
               if _graphics_ [i] [0:3] == "Nar" :
                  gr = NarPlotter
               elif _graphics_ [i] == "Gist" :
                  gr = GistPlotter
               else :
                  gr = Graphics
               _win_plotters_ [i] = gr.Plotter (n = i, display = _displays_ [i])
            if new_frame == "no" :
               _win_plotters_ [i].set_bytscl (_cmin_, _cmax_)
               _current_graph_ [i].set (plotter = _win_plotters_ [i])
               _current_graph_ [i].plot ( )
               _win_plotters_ [i].plot_text ( )
         if new_frame == "yes" :
            _object_list_ [i] = []
            _label_list_ = []
            if _ezcreset_ == 'true' :
               _set_frame_defaults_ ( )

def sf ( window = "all") :
   nf ("no", window = window)

def display (vals) :
   global _displays_
   if shapetest.is_scalar (vals) :
      _displays_ = [vals]
   else :
      _displays_ = vals

def titles (* vals) :
   global _titles_
   if len (vals) == 0 :
      _titles_ = [defbot, deftop, defleft, defright]
   elif len (vals) == 1 :
      _titles_ = [defbot, vals [0], defleft, defright]
   elif len (vals) >= 2 :
      _titles_ [0] = vals [1]
      _titles_ [1] = vals [0]
      if len (vals) >= 3 :
         _titles_ [2] = vals [2]
         if len (vals) >= 4 :
            _titles_ [3] = vals [3]
         else :
            _titles_ [3] = defright
      else :
         _titles_ [2] = defleft
         _titles_ [3] = defright

def titleb (val) :
   global _titles_
   _titles_ [0] = val

def titlet (val) :
   global _titles_
   _titles_ [1] = val

def titlel (val) :
   global _titles_
   _titles_ [2] = val

def titler (val) :
   global _titles_
   _titles_ [3] = val

_mark_dict_ = {"dot": ".", "circle": "o", "cross": "x", "plus": "+", 
              "asterisk": "*"}

def attr ( *kw, ** keywords ) :
   """attr ( <keyword arguments> ) changes the values of certain
   keywords for the next graphics commands, until changed by
   keyword arguments to a plot command (changed for that command
   only), another attr command, or an nf (if ezcreset == "true",
   then attributes revert to their default values).
   If a keyword is a frame attribute (grid, scale, labels) then
   a new graph will be plotted regardless of the setting of _ezcshow_.
   """
   global _grid_type_, _axis_scales_, _thick_, _style_, _mark_, \
          _color_, _lev_, _krange_, _lrange_, _kstyle_, _lstyle_, \
          _region_, _labels_enabled_, _default_label_, _object_list_, \
          _null_object_list_, _win_on_
   if len (kw) > 0 :
      keywords = kw [0]
   if keywords.has_key ("grid") or keywords.has_key ("scale") or \
      keywords.has_key ("labels"):
      _nf_ = 1
   else :
      _nf_ = 0
   if keywords.has_key ("region") :
      _region_ = keywords ["region"]
   if keywords.has_key ("krange") :
      _krange_ = keywords ["krange"]
   if keywords.has_key ("lrange") :
      _lrange_ = keywords ["lrange"]
   if keywords.has_key ("kstyle") :
      _kstyle_ = keywords ["kstyle"]
   if keywords.has_key ("lstyle") :
      _lstyle_ = keywords ["lstyle"]
   if keywords.has_key ("grid"):
      _grid_type_ = keywords ["grid"]
      if _grid_type_ [0:2] == "no" :
         _grid_type_ = "none"
      elif _grid_type_ == "tickonly" :
         _grid_type_ = "axes"
      elif _grid_type_ == "xy" :
         _grid_type_ = "wide"
      else : # someday will implement x and y only
         _grid_type_ = "wide"
   if keywords.has_key ("scale"):
      _axis_scales_ = keywords ["scale"]
   if keywords.has_key ("thick"):
      _thick_ = keywords ["thick"]
   if keywords.has_key ("style"):
      _style_ = keywords ["style"]
   if keywords.has_key ("mark"):
      _mark_ = keywords ["mark"]
      if _mark_dict_.has_key (_mark_) :
         _mark_ = _mark_dict_ [_mark_]
   if keywords.has_key ("labels"):
      # Ignored if not either yes or no
      lbl = keywords ["labels"]
      if type (lbl) == StringType :
         if lbl [0] == 'y' or lbl [0] == 'Y' :
            _labels_enabled_ = 1
         elif lbl [0] == 'n' or lbl [0] == 'N' :
            _labels_enabled_ = 0
   if keywords.has_key ("label"):
      _default_label_ = keywords ["label"]
   if keywords.has_key ("color"):
      _color_ = keywords ["color"]
   if keywords.has_key ("lev") :
      _lev_ = keywords ["lev"]
   # Don't do the following if no windows are open yet, sf() will die.
   if _nf_ and _win_on_ != [0, 0, 0, 0, 0, 0, 0, 0]:
      # a frame trait was changed; plot the graph.
      # all other attributes are active only for later objects
      # introduced into the graph.
      for i in _object_list_ :
         if i == [] :
            continue # (there is nothing to graph yet)
         for j in range (len (i)) :
            # either set or reset object labels
            if _labels_enabled_ == 1 :
               i [j].set (label = _label_list_ [j])
            else :
               i [j].set (label = " ")
      if _object_list_ != _null_object_list_ :
         sf ( )
 
_style_dict_ = {"solid": "solid", "dashed": "dash", "dotted" : "dot",
                "pm" : "solid", "none": "none", "dotdash" : "dashdot"}
      
def plot (y = None, x = None, ** keywords) :
   global _grid_type_, _axis_scales_, _style_, _label_list_, \
      _labels_enabled_, _default_label_, _object_list_, _style_dict_, \
      _mark_dict_, _color_, _thick_, _style_, _ezcshow_, _labels_temporarily_enabled_
   if y is None :
      attr (keywords)
      return
   if keywords.has_key ("window") :
      window = keywords ["window"]
   else :
      window = "all"
   window = _check_win_ (window, minmax = 1)
   mark = " "
   if keywords.has_key ("grid") or keywords.has_key ("scale") :
      # Force a newframe if either of these keywords is specified
      _nf_ = 1
   else :
      _nf_ = 0
   if keywords.has_key ("grid"):
      _grid_type_ = keywords ["grid"]
      if _grid_type_ [0:2] == "no" :
         _grid_type_ = "none"
      elif _grid_type_ == "tickonly" :
         _grid_type_ = "axes"
      elif _grid_type_ == "xy" :
         _grid_type_ = "wide"
      else : # someday will implement x and y only
         _grid_type_ = "wide"
   if keywords.has_key ("scale"):
      _axis_scales_ = keywords ["scale"]
   _labels_temporarily_enabled_ = _labels_enabled_
   if keywords.has_key ("labels") :
      lbl = keywords ["labels"]
      if type (lbl) == StringType :
         if lbl [0] == 'y' or lbl [0] == 'Y' :
            _labels_temporarily_enabled_ = 1
         elif lbl [0] == 'n' or lbl [0] == 'N' :
            _labels_temporarily_enabled_ = 0
   if keywords.has_key ("label") :
      lbl = keywords ["label"]
   else :
      lbl = _default_label_
   if keywords.has_key ("color") :
      col = keywords ["color"]
   else :
      col = _color_
   # if y is a vector of coordinates lbl, x, and col should be converted to
   # objects of the same size.
   if len (shape (y)) == 1 :
      if no_of_dims (lbl) == 1 :
         lbl = lbl [0]
      if no_of_dims (col) == 1 :
         col = col [0]
      if x is not None and no_of_dims (x) == 2 :
         x = x [0]
      no_of_coords = 1
   else :
      no_of_coords = shape (y) [0]
      if is_scalar (lbl) :
         lbl = [lbl] * no_of_coords
      elif len (lbl) < no_of_coords :
         lbl = lbl + [" "] * (no_of_coords - len (lbl))
      if is_scalar (col) :
         col = [col] * no_of_coords
      elif len (col) < no_of_coords :
         col = col + [_color_] * (no_of_coords - len (col))
      if x is None or type (x) == ArrayType and len (shape (x)) == 1 :
         x = [x] * no_of_coords
      elif shape (x) [0] < no_of_coords :
         x = x + [None] * no_of_coords - shape (x) [0]
   if keywords.has_key ("thick") :
      thk = keywords ["thick"]
   else :
      thk = _thick_
   if keywords.has_key ("style") :
      stl = keywords ["style"]
      if stl == "dotted" :
         mark = "."
         stl = "dot"
      elif _style_dict_.has_key (stl) :
         stl = _style_dict_ [stl]
      # non-blank mark forced "none"
      hid = 0
   else :
      if _style_dict_.has_key (_style_) :
         stl = _style_dict_ [_style_]
      else :
         stl = _style_
      hid = 0
   if keywords.has_key ("mark"):
      mark = keywords ["mark"]
      if _mark_dict_.has_key (mark) :
         mark = _mark_dict_ [mark]
      # non-blank mark forced "none"
      if mark != " " :
         hid = 0
         stl = "none"
   if (_labels_enabled_ == 1 or _labels_temporarily_enabled_) and mark == " " \
      and lbl != " " and no_of_coords == 1 and stl == "none":
      mark = lbl [0]
   for j in window:
      if no_of_coords == 1 :
         if mark == " "  and not _labels_enabled_ and not _labels_temporarily_enabled_:
            _object_list_ [j].append (Curve (y = y, x = x, color = col, # label = lbl,
               type = stl, width = thk, hide = hid))
         elif  mark == " "  and (_labels_enabled_ or _labels_temporarily_enabled_):
            _object_list_ [j].append (Curve (y = y, x = x, color = col, label = lbl,
               type = stl, marks = 1, width = thk, hide = hid))
         else :
            _object_list_ [j].append (Curve (y = y, x = x, color = col, label = lbl,
               type = stl, marks = 1, marker = mark, width = thk, hide = hid))
         _label_list_.append (lbl)
      else :
         for i in range (no_of_coords) :
            if mark == " " and not (_labels_enabled_ or _labels_temporarily_enabled_):
               _object_list_ [j].append (Curve (y = y [i], x = x [i], color = col [i], 
                   # label = lbl [i],
                   type = stl, width = thk, hide = hid))
            elif  mark == " "  and (_labels_enabled_ or _labels_temporarily_enabled_):
               _object_list_ [j].append (Curve (y = y [i], x = x [i], color = col [i],
                   label = lbl [i], type = stl, marks = 1, width = thk,
                   hide = hid))
            else :
               _object_list_ [j].append (Curve (y = y [i], x = x [i], color = col [i],
                   label = lbl [i], type = stl, marks = 1, marker = mark,
                   width = thk, hide = hid))
         _label_list_ = _label_list_ + lbl
   if _nf_ or _ezcshow_ == 'true' :
      sf ( window = window )

_bnd_dict_ = {"Yes" : 1 , "yes" : 1 , "YES" : 1 , "y" : 1 , "Y" : 1 ,
              "ON" : 1 , "On" : 1 , "on" : 1 , 1 : 1,
              "No" : 0 , "no" : 0 , "NO" : 0 , "n" : 0 , "N" : 0 ,
              "None" : 0 , "none" : 0 , "NONE" : 0 ,
              "OFF" : 0 , "off" : 0 , "Off" : 0 , 0 : 0 }

_SetMeshError_ = "SetMeshError"

def set_mesh ( ** keywords ) :
   """set_mesh (rt = <array1>, zt = <array2>, ireg = <array3>,
   ut = <array3>, vt = <array4> , pvar = <array5>, cindex = <array6>)
   defines a two-dimensional mesh for plotting. rt and zt are real,
   two-dimensional arrays of the same shape defining the mesh. ireg
   is a two-dimensional integer array defining which region of the
   mesh each quadrilateral in it belongs to. It should be the same
   shape as rt and zt, but the first row and first column are
   constrained to be 0. ut and vt are velocity components used to
   plot vector fields. pvar and cindex are mutually exclusive. pvar
   is a real array used to color a filled mesh, while cindex is a
   character array whose components specify an index into a color
   table. All arrays must be the same shape. Once set,
   these variables will define the mesh until the next set_mesh
   command.  Any variable undefined by set_mesh must be supplied
   as a keyword argument to a plot command or must have been
   previously defined by a set_mesh command. A variable defined by
   set_mesh can be overruled temporarily by a keyword argument to
   a plot command.  Variables set by set_mesh will remain defaults
   unless they are overruled by keyword arguments to plot commands
   or until they are changed permanently by another set_mesh or by
   a clear_mesh.
   """
   global _zt_, _rt_, _ut_, _vt_, _ireg_, _pvar_, _cindex_
   if keywords.has_key ("zt") :
      _zt_ = keywords ["zt"]
   else :
      _zt_ = None
   if keywords.has_key ("rt") :
      _rt_ = keywords ["rt"]
   else :
      _rt_ = None
   if keywords.has_key ("ut") :
      _ut_ = keywords ["ut"]
   else :
      _ut_ = None
   if keywords.has_key ("vt") :
      _vt_ = keywords ["vt"]
   else :
      _vt_ = None
   if keywords.has_key ("ireg") :
      _ireg_ = keywords ["ireg"]
   else :
      _ireg_ = None
   if keywords.has_key ("pvar") and keywords.has_key ("cindex") :
      raise _SetMeshError_, "Only one of 'pvar' and 'cindex' is allowed."
   if keywords.has_key ("pvar") :
      _pvar_ = keywords ["pvar"]
      _cindex_ = None
   elif keywords.has_key ("cindex") :
      _cindex_ = keywords ["cindex"]
      _pvar_ = None

def clear_mesh ( ) :
   global _zt_, _rt_, _ireg_, _pvar_, _cindex_,  _ut_, _vt_
   _zt_ = None
   _rt_ = None
   _ireg_ = None
   _pvar_ = None
   _cindex_ = None
   _ut_ = None
   _vt_ = None

def ezcx (val) :
   global _zt_
   _zt_ = val

def ezcy (val) :
   global _rt_
   _rt_ = val

def ezcpvar (val) :
   global _pvar_
   _pvar_ = val

def ezccindex (val) :
   global _cindex_
   _cindex_ = val

def ezcireg (val) :
   global _ireg_
   _ireg_ = val

def ezcu (val ) :
   global _ut_
   _ut_ = val

def ezcv (val ) :
   global _vt_
   _vt_ = val

_MeshDefError_ = "MeshDefError"
_FillError_ = "FillError"

def _minm_ (x) :
   """_minm_ (x) computes and returns something guaranteed smaller than
   the minimum element in the array x."""
   min_ = min (ravel (x))
   if min_ == 0 :
      return -1
   elif min_ > 0 :
      return min_ - .01 * min_
   else :
      return min_ + .01 * min_

_PlotfError_ = "PlotfError"
_RangeError_ = "RangeError"

def plotm (*kwds, ** keywords) :
   global _grid_type_, _axis_scales_, _cmin_, _cmax_, _krange_, _lrange_, \
          _region_, _color_bar_, _labels_enabled_, _zt_, _rt_, _ut_, _vt_, \
          _ireg_, _pvar_, _cindex_, _object_list_, _labels_temporarily_enabled_
   _labels_temporarily_enabled_ = _labels_enabled_
   if len (kwds) != 0 :
      keywords = kwds [0]
   if keywords.has_key ("window") :
      window = keywords ["window"]
   else :
      window = "all"
   window = _check_win_ (window, minmax = 1)
   if keywords.has_key ("color_bar") :
      _color_bar_ = keywords ["color_bar"]
   if keywords.has_key ("contour_plot") and keywords ["contour_plot"] != 0 :
      cp = 1
   else :
      cp = 0
   if keywords.has_key ("lev") :
      __lev_ = keywords ["lev"]
   else :
      __lev_ = _lev_
   if keywords.has_key ("fill") :
      if _bnd_dict_.has_key (keywords ["fill"]) :
         _fill_ = _bnd_dict_ [keywords ["fill"]]
      else :
         raise _FillError_, "Unknown value <" + `keywords ["fill"]` + \
            "> for keyword 'fill.'"
   else :
      _fill_ = 0
   if keywords.has_key ("zt") :
      __zt_ = keywords ["zt"]
   else :
      __zt_ = _zt_
   if keywords.has_key ("rt") :
      __rt_ = keywords ["rt"]
   else :
      __rt_ = _rt_
   if keywords.has_key ("ut") :
      __ut_ = keywords ["ut"]
   else:
      __ut_ = _ut_
   if keywords.has_key ("vt") :
      __vt_ = keywords ["vt"]
   else:
      __vt_ = _vt_
   if keywords.has_key ("ireg") :
      __ireg_ = keywords ["ireg"]
   else :
      __ireg_ = _ireg_
   if __zt_ is None and __rt_ is not None and cp == 0 :
      raise _MeshDefError_, "zt is missing from the mesh.\n" + \
         "Use set_mesh command or zt keyword in plotm command."
   elif __zt_ is not None and __rt_ is None and cp == 0 :
      raise _MeshDefError_, "rt is missing from the mesh.\n" + \
         "Use set_mesh command or rt keyword in plotm command."
   elif __zt_ is None and __rt_ is None and cp == 0 :
      raise _MeshDefError_, "Both zt and rt are missing from the mesh.\n" + \
         "Use set_mesh command or rt and zt keywords in plotm command."
   if _fill_ == 0 and cp == 0 :
      __pvar_ = _pvar_
      __cindex_ = _cindex_
   if _fill_ == 1 or cp == 1 : # look for either pvar or cindex.
      __pvar_ = None
      __cindex_ = None
      if keywords.has_key ("pvar") :
         __pvar_ = keywords ["pvar"]
      elif keywords.has_key ("cindex") :
         __cindex_ = keywords ["cindex"]
      else :
         __pvar_ = _pvar_
         __cindex_ = _cindex_
      if __pvar_ is not None :
         k = shape (__pvar_) [0]
         l = shape (__pvar_) [1]
      elif __cindex_ is not None :
         k = shape (__cindex_) [0]
         l = shape (__cindex_) [1]
      else :
         raise _MeshDefError_, "Neither pvar nor cindex has been specified."
   if __ireg_ is None and __zt_ is not None : # Default to a single region
      if _fill_ == 0 :
         k = shape (__zt_) [0]
         l = shape (__zt_) [1]
      __ireg_ = ones ( (k, l), Int)
      __ireg_ [:, 0] = 0
      __ireg_ [0, :] = 0
   _z_ = None

   __krange_ = _krange_
   __lrange_ = _lrange_
   if keywords.has_key ("krange"):
      __krange_ = keywords ["krange"]
      if len (__krange_) >= 1 :
         klo = __krange_ [0]
         if klo < 0 or klo > shape (__ireg_) [0] :
            raise _RangeError_, `klo` + " is an illegal low subscript."
         if len (__krange_) >= 2 :
            khi = __krange_ [1]
            if khi < 0 or khi > shape (__ireg_) [0] :
               raise _RangeError_, `khi` + " is an illegal high subscript."
            if len (__krange_) >= 3 :
               kstr = __krange_ [2]
               if kstr < 0 or kstr > shape (__ireg_) [0] :
                  raise _RangeError_, `kstr` + " is an illegal stride."
            else :
               kstr = 1
         else :
            khi = shape (__ireg_) [0]
      else :
         raise _RangeError_, "<" + `__krange_` + "> is an illegal range."
   if keywords.has_key ("lrange"):
      __lrange_ = keywords ["lrange"]
      if len (__lrange_) >= 1 :
         llo = __lrange_ [0]
         if llo < 0 or llo > shape (__ireg_) [1] :
            raise _RangeError_, `llo` + " is an illegal low subscript."
         if len (__lrange_) >= 2 :
            lhi = __lrange_ [1]
            if lhi < 0 or lhi > shape (__ireg_) [1] :
               raise _RangeError_, `lhi` + " is an illegal high subscript."
            if len (__lrange_) >= 3 :
               lstr = __lrange_ [2]
               if lstr < 0 or lstr > shape (__ireg_) [1] :
                  raise _RangeError_, `lstr` + " is an illegal stride."
            else :
               lstr = 1
         else :
            lhi = shape (__ireg_) [1]
      else :
         raise _RangeError_, "<" + `__lrange_` + "> is an illegal range."
   if __krange_ is not None and __lrange_ is None :
      llo = 0
      lhi = shape (__ireg_) [1]
      lstr = 1
      __lrange_ = (llo, lhi, lstr)
   elif __krange_ is None and __lrange_ is not None :
      klo = 0
      khi = shape (__ireg_) [0]
      kstr = 1
      __krange_ = (klo, khi, kstr)
   # These attributes are sticky, so set the globals
   _krange_ = __krange_
   _lrange_ = __lrange_

   if __krange_ is not None : 
      # All the data must be resized
      __ireg_ = __ireg_ [klo:khi + 1:kstr, llo:lhi + 1:lstr]
      if __zt_ is not None :
         __zt_ = __zt_ [klo:khi + 1:kstr, llo:lhi + 1:lstr]
      if __rt_ is not None :
         __rt_ = __rt_ [klo:khi + 1:kstr, llo:lhi + 1:lstr]
      if __ut_ is not None :
         __ut_ = __ut_ [klo:khi + 1:kstr, llo:lhi + 1:lstr]
         __vt_ = __vt_ [klo:khi + 1:kstr, llo:lhi + 1:lstr]
      if __pvar_ is not None :
         __pvar_ = __pvar_ [klo:khi + 1:kstr, llo:lhi + 1:lstr]
      elif __cindex_ is not None:
         __cindex_ = __cindex_ [klo:khi + 1:kstr, llo:lhi + 1:lstr]

   _cscale_ = None
   if _fill_ == 1 or cp == 1 : # look for either pvar or cindex.
      if keywords.has_key ("cscale") :
         _cscale_ = keywords ["cscale"]
      else :
         _cscale_ = "lin"
      if _cscale_ != "lin" and __pvar_ is None :
         raise _PlotfError_, \
            "pvar keyword must be present if cscale is not 'lin.'"
      if __pvar_ is not None :
         _z_ = __pvar_
      elif __cindex_ is not None :
         _z_ = __cindex_
      if keywords.has_key ("edges") :
         _edges_ = keywords ["edges"]
         if keywords.has_key ("ewidth") :
            _ewidth_ = keywords ["ewidth"]
         else :
            _ewidth_ = 1
         if keywords.has_key ("ecolor") :
            _ecolor_ = keywords ["ecolor"]
         else :
            _ecolor_ = "fg"
      else :
         _edges_ = 0
         _ewidth_ = 0
         _ecolor_ = "bg"
   
   mark = " "
   if keywords.has_key ("grid") or keywords.has_key ("scale") :
      # Force a newframe if either of these keywords is specified
      _nf_ = 1
   else :
      _nf_ = 0
   if keywords.has_key ("grid"):
      _grid_type_ = keywords ["grid"]
      if _grid_type_ [0:2] == "no" :
         _grid_type_ = "none"
      elif _grid_type_ == "tickonly" :
         _grid_type_ = "axes"
      elif _grid_type_ == "xy" :
         _grid_type_ = "wide"
      else : # someday will implement x and y only
         _grid_type_ = "wide"
   if keywords.has_key ("scale"):
      _axis_scales_ = keywords ["scale"]
   else :
      _axis_scales_ = "linlin"
   if keywords.has_key ("style") :
      stl = keywords ["style"]
      if stl == "dotted" :
         mark = "."
      elif _style_dict_.has_key (stl) :
         stl = _style_dict_ [stl]
      hid = stl == "none"
   else :
      stl = "solid"
      hid = 0
   if keywords.has_key ("bnd") and _bnd_dict_.has_key (keywords ["bnd"]) :
      bnd = _bnd_dict_ [keywords ["bnd"]]
   else :
      bnd = 0
   # default kstyle and lstyle depend on setting of bnd.
   if bnd == 1 :
      _kstyle_ = "none"
      _lstyle_ = "none"
   else :
      _kstyle_ = stl
      _lstyle_ = stl
   if keywords.has_key ("kstyle"):
      _kstyle_ = keywords ["kstyle"]
      if _style_dict_.has_key (_kstyle_) :
         _kstyle_ = _style_dict_ [_kstyle_]
   if keywords.has_key ("lstyle"):
      _lstyle_ = keywords ["lstyle"]
      if _style_dict_.has_key (_lstyle_) :
         _lstyle_ = _style_dict_ [_lstyle_]
   _inhibit_ = 0
   if _kstyle_ == "none" :
      _inhibit_ = _inhibit_ + 1
   if _lstyle_ == "none" :
      _inhibit_ = _inhibit_ + 2
   if keywords.has_key ("thick") :
      thk = keywords ["thick"]
   else :
      thk = _thick_
   if keywords.has_key ("color") :
      col = keywords ["color"]
   else :
      col = _color_
   if keywords.has_key ("labels") :
      lbl = keywords ["labels"]
   else :
      lbl = " "
   if keywords.has_key ("mark"):
      _mark_ = keywords ["mark"]
   else :
      _mark_ = " "
   if keywords.has_key ("marksize"):
      _marksize_ = keywords ["marksize"]
   if keywords.has_key ("region"):
      __region_ = keywords ["region"]
      _region_ = __region_             # sticky
   else :
      __region_ = _region_
   if keywords.has_key ("vsc"):
      __vsc_ = keywords ["vsc"]
   else :
      __vsc_ = _vsc_
   if is_scalar (__region_) and __region_ != "all" :
      __region_ = [__region_]
   elif __region_ == "all" and bnd == 1 :
      # We need to make up a list of valid region numbers.
      # This is only for the boundary case--Gist only does
      # the outside boundary if region == 0, contrary to its
      # specifications.
      __region_ = []
      n = shape (__ireg_) [0]
      m = shape (__ireg_) [1]
      for i in range (n) :
         for j in range (m) :
            if __ireg_ [i,j] != 0 and not __ireg_ [i, j] in __region_ :
               __region_.append (__ireg_ [i, j])
   filled = 0
   contours = 0
   if cp == 1 and col [0:4] == "fill" :
      contours = filled = 1
      col = "fg"
   elif cp == 1 and col == "fillnl" :
      filled = 1
      col = "fg"
   elif cp == 1 :
      contours = 1
   if cp == 1 :
      _fill_ = filled
   if __region_ != "all" :
      # Plot the indicated regions
      _region_list_ = []
      for i in range (len (__region_)) :
         if _mark_ == " " and _labels_enabled_ == 0 :
            _region_list_.append ( Region (number = __region_ [i],
               boundary = bnd, boundary_type = stl, boundary_color = col,
               inhibit = _inhibit_, type = stl, color = col,
               ktype = _kstyle_, ltype = _lstyle_, filled = _fill_ ,
               width = thk, label = lbl, hide = hid, contours = contours,
               levels = __lev_))
         elif _mark_ == " " and _labels_enabled_ == 1 :
            _region_list_.append ( Region (number = __region_ [i],
               boundary = bnd, boundary_type = stl, boundary_color = col,
               inhibit = _inhibit_, type = stl, color = col,
               width = thk, label = lbl, hide = hid, marks = 1,
               ktype = _kstyle_, ltype = _lstyle_, filled = _fill_ ,
               contours = contours, levels = __lev_))
         else :
            _region_list_.append ( Region (number = __region_ [i],
               boundary = bnd, boundary_type = stl, boundary_color = col,
               inhibit = _inhibit_, type = stl, color = col,
               width = thk, label = lbl, hide = hid, marks = 1,
               ktype = _kstyle_, ltype = _lstyle_, filled = _fill_ ,
               marker = _mark_, contours = contours,
               levels = __lev_))
      for j in window :
         _object_list_ [j].append (QuadMesh (x = __zt_, y = __rt_, ireg = __ireg_,
            vx = __vt_, vy = __ut_, scale = __vsc_, z_scale = _cscale_,
            regions = _region_list_, filled = _fill_, contours = contours,
            levels = __lev_, z = _z_))
   else :
      # Plot the whole mesh
      for j in window :
         if _mark_ == " " and _labels_enabled_ == 0 :
            _object_list_ [j].append (QuadMesh (x = __zt_, y = __rt_, ireg = __ireg_,
               vx = __vt_, vy = __ut_, scale = __vsc_, regions = __region_,
               boundary = bnd, boundary_type = stl, boundary_color = col,
               inhibit = _inhibit_, type = stl, color = col,
               ktype = _kstyle_, ltype = _lstyle_, filled = _fill_ ,
               width = thk, label = lbl, hide = hid, z = _z_,
               contours = contours, levels = __lev_, z_scale = _cscale_))
         elif _mark_ == " " and _labels_enabled_ == 1 :
            _object_list_ [j].append (QuadMesh (x = __zt_, y = __rt_, ireg = __ireg_,
               vx = __vt_, vy = __ut_, scale = __vsc_, regions = __region_,
               boundary = bnd, boundary_type = stl, boundary_color = col,
               inhibit = _inhibit_, type = stl, color = col, z_scale = _cscale_,
               width = thk, label = lbl, hide = hid, marks = 1,
               ktype = _kstyle_, ltype = _lstyle_, filled = _fill_ ,
               z = _z_, contours = contours, levels = __lev_))
         else :
            _object_list_ [j].append (QuadMesh (x = __zt_, y = __rt_, ireg = __ireg_,
               vx = __vt_, vy = __ut_, scale = __vsc_, regions = __region_,
               boundary = bnd, boundary_type = stl, boundary_color = col,
               inhibit = _inhibit_, type = stl, color = col, z_scale = _cscale_,
               width = thk, label = lbl, hide = hid, marks = 1,
               ktype = _kstyle_, ltype = _lstyle_, filled = _fill_ ,
               marker = _mark_, z = _z_, contours = contours, levels = __lev_))
   _label_list_.append (lbl)
   if _nf_ or _ezcshow_ == 'true' :
      sf ( window = window )

_PlotvError_ = "PlotvError"

def plotv (*args, **keywords) :
   """plotv ( [zt [, rt [, vt [,ut [, ireg]]]]], <keyword arguments> )
   plots vectors on a mesh."""
   rt = _rt_
   zt = _zt_
   ut = _ut_
   vt = _vt_
   ireg = _ireg_
   if not keywords.has_key ("kstyle") :
      keywords ["kstyle"] = "none"
   if not keywords.has_key ("lstyle") :
      keywords ["lstyle"] = "none"
   if len (args) > 0 :
      zt = args [0]
      keywords ["zt"] = zt
      if len (args) > 1 :
         rt = args [1]
         keywords ["rt"] = rt
         if len (args) > 2 :
            vt = args [2]
            keywords ["vt"] = vt
            if len (args) > 3 :
               ut = args [3]
               keywords ["ut"] = ut
               if len (args) > 4 :
                  ireg = args [4]
                  keywords ["ireg"] = ireg
               elif keywords.has_key ("ireg") :
                  ireg = keywords ["ireg"]
            elif keywords.has_key ("ut") :
               ut = keywords ["ut"]
         elif keywords.has_key ("vt") :
            vt = keywords ["vt"]
      elif keywords.has_key ("rt") :
         rt = keywords ["rt"]
   else :
      if keywords.has_key ("zt") :
         zt = keywords ["zt"]
      if keywords.has_key ("rt") :
         rt = keywords ["rt"]
      if keywords.has_key ("ut") :
         ut = keywords ["ut"]
      if keywords.has_key ("vt") :
         vt = keywords ["vt"]
   if rt is None :
      raise _PlotvError_, "rt has not been defined. Use set_mesh, or pass\n" + \
         "as second argument or keyword argument to plotv."
   if zt is None :
      raise _PlotvError_, "zt has not been defined. Use set_mesh, or pass\n" + \
         "as first argument or keyword argument to plotv."
   if ut is None :
      raise _PlotvError_, "ut has not been defined. Use set_mesh, or pass\n" + \
         "as fourth argument or keyword argument to plotv."
   if vt is None :
      raise _PlotvError_, "vt has not been defined. Use set_mesh, or pass\n" + \
         "as third argument or keyword argument to plotv."
   plotm ( keywords )

def plotb (** keywords) :
   """plotb ( <keyword arguments> ) is the same as plotm with the
   keyword bnd = 1, i. e., plot the boundary of the mesh. """
   keywords ["bnd"] = 1
   plotm ( keywords )

_PlotcError_ = "PlotcError"

def plotc (** keywords) :
   """plotc ( <keyword arguments> ) is similar to plotm; however, it
   adds to the keywords in such a way that plotm will do a contour plot.
   As in plotf, zt and rt may have been specified as one dimensional,
   and if so, they need to be expanded to two dimensions.
   """
   global _lev_, _pvar_
   if keywords.has_key ("pvar") :
      (k, l) = shape (keywords ["pvar"])
   elif _pvar_ is not None :
      (k, l) = shape (_pvar_)
   else :
      raise _PlotcError_ , "pvar keyword must be specified or " + \
         "pvar must be defined by set_mesh."
   if keywords.has_key ("zt") :
      x = keywords ["zt"]
      if no_of_dims (x) == 1 :
         if len (x)  != k :
            raise _PlotcError_, "Length of zt <" + `len (x)` + \
               "> does not match first dimension of pvar <" + `k` + ">."
         keywords ["zt"] = multiply.outer (x, ones (l))
      elif shape (y) != (k, l) :
         raise _PlotcError_, "zt has a funny shape: " + `shape (y)` + \
            "."
   if keywords.has_key ("rt") :
      y = keywords ["rt"]
      if no_of_dims (y) == 1 :
         if len (y)  != l : 
            raise _PlotcError_, "Length of rt <" + `len (y)` + \
               "> does not match first dimension of pvar <" + `l` + ">."
         keywords ["rt"] = multiply.outer (ones (k), y)
      elif shape (x) != (k, l) :
         raise _PlotcError_, "rt has a funny shape: " + `shape (x)` + \
            "."
   keywords ["kstyle"] = "none"
   keywords ["lstyle"] = "none"
   keywords ["contour_plot"] = 1
   if keywords.has_key ("lev") :
      lev = keywords ["lev"]
   else :
      lev = _lev_
   if lev == "lin" or lev == "linear" :
      lev = deflev
   if lev == "log" or type (lev) == IntType and lev < 0 :
      keywords ["cscale"] = "log"
      if lev == "log" :
         lev = deflev
      else :
         lev = - lev
   keywords ["lev"] = lev
   plotm (keywords)

_PlotzError_ = "PlotzError"

def plotz (* args, **keywords) :
   """plotz (z [,x [,y]], <keyword arguments> ) is the non mesh-oriented
   contour plot. keyword pvar will be set to z, rt will be set to x
   (if present), and zt will be set to y (if present). If x or y or
   both are 1d, they will be expanded to a 2d mesh. This is so that
   plotm can be used."""
   global _lev_
   keywords ["kstyle"] = "none"
   keywords ["lstyle"] = "none"
   keywords ["contour_plot"] = 1
   if len (args) != 0 :
      z = args [0]
   elif keywords.has_key ("z") :
      z = keywords ["z"]
      del keywords ["z"]
   else :
      raise _PlotzError_, "No z argument was given for plotz!!"
   if len (args) > 1 :
      x = args [1]
      if len (args) > 2 :
         y = args [2]
      else :
         y = None
   else :
      x = None
      y = None
   if x is not None and y is not None and \
      len (x.shape) == len (y.shape) == len (z.shape) == 1 and \
      x.shape == y.shape == z.shape:
      # random data; turn into a regular mesh.
      xmax = max (x)
      xmin = min (x)
      ymax = max (y)
      ymin = min (y)
      nx = min (int (sqrt (len (x))) + 1, 50)
      ny = nx
      dx = (xmax - xmin) / nx
      dy = (ymax - ymin) / ny
      xcplot = span (xmin, xmax - dx, nx)
      ycplot = span (ymin, ymax - dy, ny)
      del xmin, xmax, ymin, ymax, dx, dy
      xt = subtract.outer (x, x)
      yt = subtract.outer (y, y)
      aa = sqrt (xt * xt + yt * yt + rsq)
      del xt, yt
      alpha = array (z, copy = 1)
      alpha = solve_linear_equations (aa, alpha)
      z = mfit (alpha, x, xcplot, y,
         ycplot, rsq)
      del aa, alpha, rsq
      # Expand coordinates to 2d arrays to match zcplot
      x = multiply.outer (xcplot, ones (ny, Float))
      y = multiply.outer (ones (nx, Float), ycplot)
      del xcplot, ycplot, nx, ny
   keywords ["pvar"] = z
   (k, l) = shape (z)
   if x is not None and len (shape (x)) == 1 :
      if len (x) != k :
         raise _PlotzError_, "length of x <" + `len (x)` + \
           "> does not match first dimension of z <" + `k` + ">."
      x = multiply.outer (x, ones (l, Float))
   elif x is None :
      x = multiply.outer (arange (1, k + 1, typecode = Float), ones (l, Float))
   if y is not None and len (shape (y)) == 1 :
      if len (y) != l :
         raise _PlotzError_, "length of y <" + `len (y)` + \
           "> does not match second dimension of z <" + `l` + ">."
      y = multiply.outer (ones (k, Float), y)
   elif y is None :
      y = multiply.outer (ones (k, Float), arange (1, l + 1, typecode = Float))
   keywords ["zt"] = y
   keywords ["rt"] = x
   if keywords.has_key ("lev") :
      lev = keywords ["lev"]
   else :
      lev = _lev_
   if lev == "lin" or lev == "linear" :
      lev = deflev
   if lev == "log" or type (lev) == IntType and lev < 0 :
      keywords ["cscale"] = "log"
      if lev == "log" :
         lev = deflev
      else :
         lev = - lev
   keywords ["lev"] = lev
   plotm (keywords)

_PlotiError_ = "PlotiError"

def ploti (* pvar, ** keywords) :
   """ploti (pvar, <keyword arguments> ) is the cell array plot."""
   global _object_list_, _ezcshow_, _pvar_, _labels_enabled_, \
      _labels_temporarily_enabled_
   _labels_temporarily_enabled_ = _labels_enabled_
   if _object_list_ != _null_object_list_ :
      raise _PlotiError_, "Cell array plot can't be superposed with others."
   if len (pvar) > 0 :
      pvar = pvar [0]
   elif keywords.has_key ("pvar") :
      pvar = keywords ["pvar"]
   elif _pvar_ is not None :
      pvar = _pvar_
   else :
      pvar = _cindex_
   if keywords.has_key ("window") :
      window = keywords ["window"]
   else :
      window = "all"
   window = _check_win_ (window, minmax = 1)
   if keywords.has_key ("style") :
      hid = keywords ["style"] == "none"
   else :
      hid = 0
   if keywords.has_key ("labels") :
      lbl = keywords ["labels"]
   else :
      lbl = " "
   for j in window:
      _object_list_ [j].append ( CellArray (z = pvar, hide = hid, label = lbl))
   if _ezcshow_ == 'true' :
      sf ( window = window )

def plotf (*args, ** keywords) :
   """plotf ( arg1 [,arg2 [,arg3 [,arg4 ]]] [,<keyword arguments>] )
   is the same as plotm with the keyword fill = 1, i. e., plot a
   filled mesh. arg1 if integer is an array which gives color
   numbers; if real, gives the values of a real variable which
   will be used to interpolate into the current palette. arg2,
   if present, is zt. arg3, if present, is rt. arg4, if present,
   is ireg."""
   global _zt_, _rt_, _pvar_
   n = len (args)
   __zt_ = _zt_
   __rt_ = _rt_
   if n != 0 :
      if type (args[0][0][0]) == IntType :
         plvar = keywords ["cindex"] = array(args [0],'b') # Gist knows chars
      else :
         plvar = keywords ["pvar"] = args [0]
         if n > 1 :
            __zt_ = keywords ["zt"] = args [1]
         if n > 2 :
            __rt_ = keywords ["rt"] = args [2]
         if n > 3 :
            keywords ["ireg"] = args [3]
   elif keywords is not None and not keywords.has_key ("cindex") \
      and not keywords.has_key ("pvar") and _pvar_ is None :
      raise _PlotfError_, "plotf requires at least one non-keyword " + \
         "argument\n or else one of the keywords 'cindex' or 'pvar.'"
   elif not keywords.has_key ("cindex") and not keywords.has_key ("pvar") \
        and _pvar_ is not None :
      plvar = keywords ["pvar"] = _pvar_
   elif keywords.has_key ("cindex") :
      plvar = keywords ["cindex"]
   else :
      plvar = keywords ["pvar"]
   # zt and rt might not be specified, or else thay are vectors.
   # They must be 2d and match plvar in shape.
   (k, l) = shape (plvar)
   if shape (__zt_) != shape (plvar) :
      if __zt_ is None:
         __zt_ = multiply.outer (arange (1, k + 1, typecode = Float), ones (l, Float))
      elif no_of_dims (__zt_) == 1 :
         if len (__zt_) != k :
            raise _PlotfError_, "Length of zt <" + `len (__zt_)` + \
               "> does not match first dimension of pvar <" + `k` + ">."
         __zt_ = multiply.outer (__zt_, ones (l, Float))
      else :
         raise _PlotfError_, "zt has a funny shape: " + `shape (__zt_)` + \
            "."
      keywords ["zt"] = __zt_
   if shape (__rt_) != shape (plvar) :
      if __rt_ is None:
         __rt_ = multiply.outer (ones (k, Float), arange (1, l + 1, typecode = Float))
      elif no_of_dims (__rt_) == 1 :
         if len (__rt_) != l :
            raise _PlotfError_, "Length of rt <" + `len (__rt_)` + \
               "> does not match second dimension of pvar <" + `l` + ">."
         __rt_ = multiply.outer (ones (k, Float), __rt_)
      else :
         raise _PlotfError_, "rt has a funny shape: " + `shape (__rt_)` + \
            "."
      keywords ["rt"] = __rt_
   keywords ["fill"] = 1
   plotm ( keywords )

def frame ( *args , ** keywords ) :
   """frame ( [xmin [, xmax [, ymin [, ymax]]]] [, <keyword srguments> )
   has the effect of changing the axis limits of the current graph.
   The limits may be expressed by positional arguments or by keyword
   arguments. The limits are effective immediately on the current
   graph and disappear on the next one. Any limits not set by this
   command will default to the maximum and the minimum of the data.
   """
   global _axis_limits_
   if keywords.has_key ("window") :
      window = keywords ["window"]
   else :
      window = "all"
   window = _check_win_ (window, minmax = 1)
   if keywords.has_key ("xmin") :
      _axis_limits_ [0] [0] = keywords ["xmin"]
   elif len (args) > 0 :
      _axis_limits_ [0] [0] = args [0]
   else :
      _axis_limits_ [0] [0] = "d"
   if keywords.has_key ("xmax") :
      _axis_limits_ [0] [1] = keywords ["xmax"]
   elif len (args) > 1 :
      _axis_limits_ [0] [1] = args [1]
   else :
      _axis_limits_ [0] [1] = "d"
   if keywords.has_key ("ymin") :
      _axis_limits_ [1] [0] = keywords ["ymin"]
   elif len (args) > 2 :
      _axis_limits_ [1] [0] = args [2]
   else :
      _axis_limits_ [1] [0] = "d"
   if keywords.has_key ("ymax") :
      _axis_limits_ [1] [1] = keywords ["ymax"]
   elif len (args) > 3 :
      _axis_limits_ [1] [1] = args [3]
   else :
      _axis_limits_ [1] [1] = "d"
   sf ( window = window )

def fr ( *args , ** keywords ) :
   """fr ( [xmin [, xmax [, ymin [, ymax]]]] [, <keyword srguments> )
   is the equivalent of nf ( ) followed by frame with the same
   set of arguments."""
   global _axis_limits_
   if keywords.has_key ("window") :
      window = keywords ["window"]
   else :
      window = "all"
   window = _check_win_ (window, minmax = 1)
   nf ( window = window )
   if keywords.has_key ("xmin") :
      _axis_limits_ [0] [0] = keywords ["xmin"]
   elif len (args) > 0 :
      _axis_limits_ [0] [0] = args [0]
   else :
      _axis_limits_ [0] [0] = "d"
   if keywords.has_key ("xmax") :
      _axis_limits_ [0] [1] = keywords ["xmax"]
   elif len (args) > 1 :
      _axis_limits_ [0] [1] = args [1]
   else :
      _axis_limits_ [0] [1] = "d"
   if keywords.has_key ("ymin") :
      _axis_limits_ [1] [0] = keywords ["ymin"]
   elif len (args) > 2 :
      _axis_limits_ [1] [0] = args [2]
   else :
      _axis_limits_ [1] [0] = "d"
   if keywords.has_key ("ymax") :
      _axis_limits_ [1] [1] = keywords ["ymax"]
   elif len (args) > 3 :
      _axis_limits_ [1] [1] = args [3]
   else :
      _axis_limits_ [1] [1] = "d"
   sf ( window = window )

def text (str, rt, zt, size, tosys = 1, color = "fg", window = "all") :
   """text (str, zt, rt, size [, tosys]) plots the specified string starting at
   (zt, rt) in point size 'size.' If tosys is specified and nonzero,
   then coordinates are not window-relative."""
   global _text_, _text_pos_, _text_color_, _text_size_, _tosys_, _ezcshow_
   _text_.append (str)
   _text_pos_.append ([rt, zt])
   _text_color_.append (color)
   _text_size_.append (size)
   _tosys_.append (tosys)
   if _ezcshow_ == 'true' :
      sf ( window = window )

def list_devices ( ) :
   global _win_on_, _win_plotters_, _cgm_file_, _cgm_, _ps_file_, _ps_, \
      _cgm_plotter_, _ps_plotter_, _graphics_, _displays_
   if _win_on_ == [0, 0, 0, 0, 0, 0, 0, 0] :
      print "No devices are currently open."
      return
   for i in range (len (_win_plotters_) ) :
      if i == _cgm_plotter_ :
         if _cgm_ == "yes" :
            op = "on"
         else :
            op = "off"
         print "device %d: cgm file '%s', currently %s." % \
            (i, _cgm_file_, op)
      elif i == _ps_plotter_ :
         if _ps_ == "yes" :
            op = "on"
         else :
            op = "off"
         print "device %d: ps file '%s', currently %s." % \
            (i, _ps_file_, op)
      elif _win_on_ [i] :
         if _graphics_ [i] [0:3] == "Nar" or Graphics == NarPlotter :
            gr = "Narcisse"
         elif _graphics_ [i] == "Gist" or Graphics == GistPlotter :
            gr = "Gist"
         else :
            gr = Graphics
         print "device %d: %s Xwindow for display '%s'." % \
            (i, gr, _displays_ [i])
   return
