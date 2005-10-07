# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

from gist import *
from shapetest import *
from Numeric import *
from scipy_base.fastumath import *
from graftypes import *
from movie import *
from pl3d import *
from plwf import *
from slice3 import *
from colorbar import *
from mesh3d import *
from arrayfns import *
 
def minmax ( x ) :
   """minmax (x) where x is a two-dimensional array computes the minimum and
   maximum values in the array and returns them as a list [min, max]. I use
   this routine because there are apparently some circumstances in which Gist
   fails to calculate default axis limits correctly.
   """
   return [floor (min (ravel (x))), ceil (max (ravel (x)))]

def true_minmax ( x ) :
   return [min (ravel (x)), max (ravel (x))]

class Plotter :
   """Plotter ( <name> , <keyword arguments> ) creates a Gist Plotter
   object. The <name>, if non-blank, specifies where the display is to
   occur; if blank, it will be your default display; and if null (""),
   there will be no window at all. The keywords
   and their defaults are as follows:
     n (0) -- the number of the graphics window (0 to 7 are allowed).
              each plotter object corresponds to a separate window.
     dpi (100) -- the size of the window wanted. 100 and 75 are allowed;
                  100 is the larger size.
     wait (1) -- used to make sure everything is plotted before
                 changing frames.
     private (0) -- use a common colormap.
     hcp -- if not present, use default hardcopy file used by all
            windows. If present, names a file unique to this window.
     dump (0) -- if 1, dumps the color palette at the beginning of
                 each page of hardcopy output, otherwise converts
                 to grey scale.
     legends (0) -- controls whether (1) or not (0) curve legends
                    are dumped to the hardcopy.
     style ("work.gs") -- name of a Gist style sheet.
   """

   def type (self) :
      return GistType
 
   def _set_defaults ( self ) :
      self._n = -1
      self._dpi = 100
      self._wait = 1
      self._private = 0
      self._dump = 0
      self._legends = 0
      self._style = "work.gs"
      self._xscale = 0 # 0 for linear, 1 for log
      self._yscale = 0 # 0 for linear, 1 for log
      self._xmin = "e"
      self._xmax = "e"
      self._ymin = "e"
      self._ymax = "e"
      self._xmindefault = 1
      self._ymindefault = 1
      self._xmaxdefault = 1
      self._ymaxdefault = 1
      self._text = []
      self._text_color = []
      self._text_size = []
      self._text_pos = []
      self._tosys = []
      self._xyequal = 0
      self._cmin = 0.0
      self._cmax = 0.0
      self._z_c_switch = None
      # (ZCM 3/26/97) Added this so that plmesh need be called
      # only once for a given quadmesh.
      self._quadmesh = None
      self._current_palette = None
      self.making_movie = 0
      self.equal_axes = 0

   from os import environ

   _DisplayError = "DisplayError"

   def __init__ ( self , filename = ' ' , ** kw ) :
      self._set_defaults ( )
      for k in kw.keys ( ) :
         setattr ( self , "_" + k , kw [k] )
      # if there is a window with this number, kill it and its hard copy.
      if (self._n == -1) : # no number specified so use 0
         winkill (0)
         self._n = 0
      else :
         winkill ( self._n )
      if filename == " " and hasattr (self, "_display") :
         filename = self._display
      if filename == " " :
         try :
            filename = self.environ ["DISPLAY"]
         except :
            raise self._DisplayError, \
               "Gist does not know where to display your plots.\n" + \
               "You need to exit python and set your DISPLAY variable,\n" + \
               "or else instantiate your Graph object using the 'display'\n" + \
               "keyword argument."
      if (filename == "none" or filename == "" or filename is None) and \
         hasattr ( self , "_hcp" ) :
         window ( self._n , display = "" , hcp = self._hcp ,  dpi = self._dpi ,
                  wait = self._wait , legends = self._legends ,
                  style = self._style )
         hcpon ()
      elif hasattr ( self , "_hcp" ) and self._hcp != "" and filename != " " :
         hcp_file (self._hcp)
         window ( self._n , display = filename , dpi = self._dpi ,
                  wait = self._wait , private = self._private ,
                  hcp = self._hcp , dump = self._dump ,
                  legends = self._legends , style = self._style )
         hcpon ()
      elif hasattr( self , "_hcp" ) and self._hcp != "" and filename == " " :
         hcp_file (self._hcp)
         window ( self._n , dpi = self._dpi , wait = self._wait ,
                  private = self._private , hcp = self._hcp ,
                  dump = self._dump , legends = self._legends ,
                  style = self._style )
         hcpon ()
      elif filename != " " :
         window ( self._n ,  display = filename , dpi = self._dpi ,
                  wait = self._wait , private = self._private ,
                  legends = self._legends , style = self._style )
      else :
         window ( self._n , dpi = self._dpi , wait = self._wait ,
                  private = self._private ,
                  legends = self._legends , style = self._style )
      self._filename = filename
      if filename != "none" and filename != "" :
         print "Opening graphics window " + `self._n` + "."
      self._open = 1
      self._GistError = "GistError"
 
   def close ( self ) :
      "close () closes the connection to Gist."
      if self._open :
         if (self._filename == "none" or self._filename == "") \
            and hasattr ( self , "_hcp" ) :
            if window () != self._n :
               window ( self._n )
            hcp ( ) # Make sure last frame goes out
            fma ( )
            print "Closing hard copy file " + hcp_finish (self._n) + "."
            winkill ( self._n )
         else :
            winkill ( self._n )
            print "Closing graphics window " + `self._n` + "."
         self._open = 0
 
   def __del__ ( self ) :
      self.close ( )
 
   def query ( self ) :
       if not self._open :
          return -1
       else :
          return 1

   narcisse_to_gist_cc_dict = { "absolute" : "gray.gp" , "binary" : "gray.gp" ,
       "bluegreen" : "stern.gp" , "default" : "rainbow.gp" ,
       "negative" : "yarg.gp" , "positive" : "gray.gp" ,
       "rainbow" : "rainbow.gp" , "rainbowhls" : "rainbow.gp" ,
       "random" : "earth.gp" , "redblue" : "heat.gp" , "redgreen" : "heat.gp" ,
       "shifted" : "earth.gp" , 0 : "gray.gp" , 1 : "gray.gp" , 2 : "stern.gp" ,
        3 : "gray.gp" , 4 : "yarg.gp" , 5 : "gray.gp" , 6 : "rainbow.gp" ,
        7 : "rainbow.gp" , 8 : "earth.gp" , 9 : "heat.gp" , 10 : "heat.gp" ,
        11 : "earth.gp" }

   _legal_card = { 'earth.gp' : 1, 'stern.gp' : 1, 'rainbow.gp' : 1,
                   'heat.gp' : 1, 'gray.gp' : 1, 'yarg.gp' : 1 }

   def set_style ( self, st ) :
       """set_style ( st ) sets the style sheet to a new value"""
       self._style = st
       window (style = st)

   def set_distance ( self, * val ) :
       """set_distance (d) sets the viewing distance depending on d,
       a number between 0 and 20. 0 signifies infinite distance;
       20 signals very close. These need to be translated to true
       distance values for Gist.
       """
       if len (val) == 0 or val [0] <= 0 or val [0] >= 20 :
          setz3 (None)
       else :
          if val [0] < 0.1 :
             setz3 (1.e6)
          else :
             setz3 ( - 50251.206 * val [0] + 1005025.1206)

   def set_bytscl ( self, cmin, cmax ) :
       """set_bytscl ( cmin, cmax ) sets _cmin and _cmax so that
       bytscl will be called for the next plf command.
       """
       self._cmin = cmin
       self._cmax = cmax

   def set_color_card ( self , cd , now = 0) :
       """set_color_card ( palette ) indicates a predefined color card
       for a plot. Standard ones are 'earth.gp', 'stern.gp',
       'rainbow.gp', 'heat.gp', gray.gp', and 'yarg.gp'. If a Narcisse
       color card is specified, we do a translation. If an illegal one
       is specified, use the default.
       """
       hh = cd
       if self.narcisse_to_gist_cc_dict.has_key (cd) :
          hh = self.narcisse_to_gist_cc_dict [cd]
       elif not self._legal_card.has_key (cd) :
          print "Gist warning:" , cd , "is not a legal Gist color card."
          hh = "rainbow.gp"
       self._palette = hh
       if now == 1 :
          palette ( self._palette )
          del self._palette

   def send_color_card ( self ) :
       if hasattr (self, "_palette") :
          if self._current_palette != self._palette :
             palette ( self._palette )
             self._current_palette = self._palette
             del self._palette

   def set_titles ( self , * vals ) :
       """set_titles ('bottom', 'top', 'left', 'right')
       All arguments are optional. Missing ones default to ' '.
       Gist uses plt to plot the titles; values must be saved till
       the graph is plotted.
       """
       if len ( vals ) == 0 or len ( vals [0] ) == 0 :
          self._titles = [" ", " ", " ", " "]
       elif is_scalar ( vals [0] ) :
          if len (vals) == 1 :
             self._titles = [vals [0]] + [" ", " ", " "]
          elif len (vals) == 2 :
             self._titles = [vals [0], vals [1]] + [" ", " "]
          elif len (vals) == 3 :
             self._titles = [vals [0], vals [1], vals [2]] + [" "]
          else :
             self._titles = [vals [0], vals [1], vals [2], vals [3]]
       elif len ( vals[0] ) == 1 :
          self._titles = vals[0] + [" ", " ", " "]
       elif len ( vals[0] ) == 2 :
          self._titles = vals[0] + [" ", " "]
       elif len ( vals[0] ) == 3 :
          self._titles = vals[0] + [" "]
       else :
          self._titles = vals[0]
 
   def set_title_colors ( self , *vals ) :
       """set_title_colors (bottom_color, top_color, left_color, right_color)
         All arguments are optional, integers from 0 to 63 representing
         a color in some color map. Missing arguments default
         to foreground."""
       if len ( vals ) == 0 :
          self._title_colors = ["fg", "fg", "fg", "fg"]
       elif is_scalar ( vals[0] ) :
          if len (vals) == 1 :
             self._title_colors = [vals [0]] + ["fg", "fg", "fg"]
          elif len (vals) == 2 :
             self._title_colors = [vals [0], vals [1]] + ["fg", "fg"]
          elif len (vals) == 3 :
             self._title_colors = [vals [0], vals [1], vals [2]] + ["fg"]
          else :
             self._title_colors = [vals [0], vals [1], vals [2],vals [3]]
       elif len ( vals[0] ) == 1 :
          self._title_colors = vals[0] + ["fg", "fg", "fg"]
       elif len ( vals[0] ) == 2 :
          self._title_colors = vals[0] + ["fg", "fg"]
       elif len ( vals[0] ) == 3 :
          self._title_colors = vals[0] + ["fg"]
       else :
          self._title_colors = vals[0]

 
   def _lengthen_if_necessary ( self , n ) :
       while len (self._text) <= n :
          self._text.append ( " " )
          self._text_color.append ( "fg" )
          self._text_size.append ( 14.0 )
          self._text_pos.append ( [0., 0.] )
          self._tosys.append ( 0 )
 
   # The following five functions set various aspects of the text
   # array. In Gist, you must call plot_text in order to have the
   # text appear on the graph.
   def set_text ( self, t , n ) :
       "set_text (str, ix) sets the ix'th text to str."
       self._lengthen_if_necessary ( n )
       self._text [n] = t
       return

   def set_tosys ( self, v, n ) :
       "set_tosys ( v, n ) sets the n'th tosys to v."
       self._lengthen_if_necessary ( n )
       self._tosys [n] = v
       return

   def add_text (self, str, x, y, size, color="fg", tosys = 1) :
       """add_text (str, x, y, size [, color]) adds a text to a graph."""
       self._text.append (str)
       self._text_pos.append ([x, y])
       self._text_size.append (size)
       self._text_color.append (self._figure_color (color))
       self._tosys.append (tosys)

   def clear_text ( self ) :
       "clear_text ( ) makes sure that there is no text to be put out."
       self._text = []
       self._text_pos = []
       self._text_size = []
       self._text_color = []
       self._tosys = []
 
   def set_text_color ( self, n1, n2 ) :
       """set_text_color (col, ix) sets the ix'th text color to col"""
       self._lengthen_if_necessary ( n2 )
       self._text_color [n2] = self._figure_color (n1)
       return
 
   def set_text_size ( self , n , m ) :
       """set_text_size (sz, ix) sets the ix'th text size to sz.
       Text sizes mean different things in Gist and Narcisse, so
       the user has to be aware of this."""
       self._lengthen_if_necessary ( m )
       self._text_size [m] = n
       return
 
   def set_text_pos ( self , n, m, mm , tosys = 0) :
       """set_text_pos (x, y, ix) positions the ix'th text at (x, y),
       which are real numbers between 0 and 1 giving relative
       position in the graphics window."""
       self._lengthen_if_necessary ( mm )
       self._text_pos [mm] = [n, m]
       return

   def new_frame ( self ) :
       fma ( )
 
   def plot_text ( self ) :
       """plot_text ( ) sends the accumulated texts out to the
       graph."""
       if len ( self._text ) == 0 :
           return
       for n in range (len ( self._text )) :
          if self._text [n] != '' and self._text [n] != ' ' :
             plt ( self._text [n] , self._text_pos [n][0] ,
                   self._text_pos [n][1] , color = self._text_color [n] ,
                   height = self._text_size [n] , tosys = self._tosys [n])
       return
 
   def freeze_graph ( self ) :
       pass
 
   def set_no_concat ( self ) :
       pass
 
   def set_axis_labels ( self , x ) :
       pass

   def set_3d_grid_type (self, t) :
       pass

   def set_grid_type ( self , * val ) :
       """set_grid_type ( string ) determines how intrusive the axes
       and grids are. The legal arguments are:
       'none'--no axes and grids are drawn.
       'axes'--axes with tick marks.
       'wide'--widely spaced grid in x and y (2d or 3d).
       'full'--narrowly spaced grid in x and y (2d or 3d).
           'full' is not currently implemented.
       If no argument is specified, the default is 'axes'."""
 
       if len ( val ) > 1 :
          raise  self._GistError, "Too many arguments to set_grid_type."
       if len ( val ) == 0 or val [0] == "axes" :
          self._gridtype = 2
       elif val [0] == "none" :
          self._gridtype = 0
       elif val [0] == "wide" :
          self._gridtype = 1
 
   def set_axis_lin ( self , ax ) :
       """set_axis_lin (ax) where ax can be 'x', 'y', or 'all'.
       The specified axis will have a linear scale."""
       if ( ax == "x" ) :
          self._xscale = 0
       elif ( ax == "y" ) :
          self._yscale = 0
       elif ( ax == "all" ) :
          self._xscale = 0
          self._yscale = 0
       else :
          return
       #logxy ( self._xscale, self._yscale )
 
   def set_axis_log ( self , ax ) :
       """set_axis_log (ax) where ax can be 'x', 'y', or 'all'.
       The specified axis will have a logarithmic scale."""
       if ( ax == "x" ) :
          self._xscale = 1
       elif ( ax == "y" ) :
          self._yscale = 1
       elif ( ax == "all" ) :
          self._xscale = 1
          self._yscale = 1
       else :
          return
       #logxy ( self._xscale, self._yscale )
 
   def set_linlin ( self ) :
       "set_linlin () sets both x and y axes to linear scale."
       self._xscale = 0
       self._yscale = 0
       #logxy ( 0 , 0 )
 
   def set_linlog ( self ) :
       'set_linlog () sets x axis to linear, y axis to logarithmic.'
       self._xscale = 0
       self._yscale = 1
       #logxy ( 0 , 1 )
 
   def set_loglin ( self ) :
       'set_loglin () sets x axis to logarithmic, y axis to linear.'
       self._xscale = 1
       self._yscale = 0
       #logxy ( 1 , 0 )
 
   def set_loglog ( self ) :
       'set_loglog () sets both x and y axes to logarithmic scale.'
       self._xscale = 1
       self._yscale = 1
       #logxy ( 1 , 1 )
 
   def set_axis_max ( self , ax , * val1 ) :
       """set_axis_max (ax, val) where ax is 'x' or 'y'.
          The maximum of the specified axis will be set to val.
          val should be a PyFloat object, or, if omitted, the
          default will be chosen.
       """
       # (6/13/96) In all the axis limits routines, I make note whether
       # axis defaults have been asked for, and set flags accordingly.
       # For flags which are set, I will calculate defaults myself.
       # This is because I have found various instances in which
       # Gist fails to do so correctly. The calculation, unfortunately,
       # has to wait until curves (etc.) are sent to be plotted.
       if ax == "x" :
          if len ( val1 ) == 0 :
             self._xmax = "e"
          else :
             self._xmax = val1 [0]
          if self._xmin == 0 and ( self._xmax == 0 or self._xmax == "e" ) :
             self._xmin = "e"
             self._xmax = "e"
          if self._xmin == "e" : self._xmindefault = 1
          else :                 self._xmindefault = 0
          if self._xmax == "e" : self._xmaxdefault = 1
          else :                 self._xmaxdefault = 0
       elif ax == "y" :
          if len ( val1 ) == 0 :
             self._ymax = "e"
          else :
             self._ymax = val1 [0]
          if self._ymin == 0 and ( self._ymax == 0 or self._ymax == "e" ) :
             self._ymin = "e"
             self._ymax = "e"
          if self._ymin == "e" : self._ymindefault = 1
          else :                 self._ymindefault = 0
          if self._ymax == "e" : self._ymaxdefault = 1
          else :                 self._ymaxdefault = 0
 
   def set_axis_min ( self , ax , * val1 ) :
       """set_axis_min (ax, val) where ax is 'x' or 'y'.
          The minimum of the specified axis will be set to val.
          val should be a PyFloat object, or, if omitted, the
          default will be chosen.
       """
       if ax == "x" :
          if len ( val1 ) == 0 :
             self._xmin = "e"
          else :
             self._xmin = val1 [0]
          if self._xmin == 0 and self._xmax == 0 :
             self._xmin = "e"
             self._xmax = "e"
          if self._xmin == "e" : self._xmindefault = 1
          else :                 self._xmindefault = 0
          if self._xmax == "e" : self._xmaxdefault = 1
          else :                 self._xmaxdefault = 0
       elif ax == "y" :
          if len ( val1 ) == 0 :
             self._ymin = "e"
          else :
             self._ymin = val1 [0]
          if self._ymin == 0 and self._ymax == 0 :
             self._ymin = "e"
             self._ymax = "e"
          if self._ymin == "e" : self._ymindefault = 1
          else :                 self._ymindefault = 0
          if self._ymax == "e" : self._ymaxdefault = 1
          else :                 self._ymaxdefault = 0
 
   def set_default_axes_limits ( self , * h ) :
       '''set_default_axes_limits () sets Gist to compute the maximum
       and minimum of the axes depending on the data.'''
 
       self._xmin = "e"
       self._xmax = "e"
       self._ymin = "e"
       self._ymax = "e"
       self._xmindefault = 1
       self._ymindefault = 1
       self._xmaxdefault = 1
       self._ymaxdefault = 1
 
   def set_x_axis_limits ( self , val1 , * val2i ) :
       '''set_x_axis_limits (min, max) sets the limits on the x axis to
       the specified (pyFloat) sizes. If max is missing, take the default.'''
       self._xmin = val1
       if len (val2i) == 0:
          self._xmax = "e"
       else :
          self._xmax = val2i [0]
       if self._xmin == 0 and ( self._xmax == 0 or self._xmax == "e" ) :
          self._xmin = "e"
          self._xmax = "e"
       if self._xmin == "e" : self._xmindefault = 1
       else :                 self._xmindefault = 0
       if self._xmax == "e" : self._xmaxdefault = 1
       else :                 self._xmaxdefault = 0
 
   def set_y_axis_limits ( self , val1 , * val2i ) :
       '''set_y_axis_limits (min, max) sets the limits on the x axis to
       the specified (pyFloat) sizes. If max is missing, take the default.'''
       self._ymin = val1
       if len (val2i) == 0:
          self._ymax = "e"
       else :
          self._ymax = val2i [0]
       if self._ymin == 0 and ( self._ymax == 0 or self._ymax == "e" ) :
          self._ymin = "e"
          self._ymax = "e"
       if self._ymin == "e" : self._ymindefault = 1
       else :                 self._ymindefault = 0
       if self._ymax == "e" : self._ymaxdefault = 1
       else :                 self._ymaxdefault = 0

   def set_z_axis_limits ( self, val1, val2 ) :
       '''set_z_axis_limits (min, max) sets the limits on the z axis to
       the specified (pyFloat) sizes. If val1 and val2 are both 0.0,
       we assume default calculations; otherwise, graph only the part
       of the surface whose z coordinate is between val1 and val2.'''
       if val1 == val2 == 0.0 :
          self._zmin = self._zmax = None
       else :
          self._zmin = val1
          self._zmax = val2

   def set_c_axis_limits ( self, val1, val2 ) :
       '''set_c_axis_limits (min, max) sets the limits on the c values to
       the specified (pyFloat) sizes. If val1 and val2 are both 0.0,
       we assume default calculations; otherwise, graph only the part
       of the surface whose c value is between val1 and val2.'''
       if val1 == val2 == 0.0 :
          self._cmin = self._cmax = None
       else :
          self._cmin = val1
          self._cmax = val2

   def set_xyequal ( self ) :
       self._xyequal = 1

   def reset_xyequal ( self ) :
       self._xyequal = 0
 
   # Translation table from Narcisse colors to gist
   narcisse_to_gist_col = {
                   "foreground" : "fg",
                   "background" : "bg" , "orange" : "yellow" ,
                   "purple" : "magenta" }
   legal_gist_colors = { "bg" : 1 , "fg" : 1 , "black" : 1 , "white" : 1 ,
                         "red" : 1 , "green" : 1 , "blue" : 1 , "cyan" : 1 ,
                         "magenta" : 1 , "yellow" : 1 }

   def _figure_color ( self, col ) :
       """_figure_color ( col ) does the best job it can to return
       a correct color. If the value is legal for Gist (even though
       it may mean something else in another system) then it is
       returned unchanged. If it is a Narcisse value, it is converted
       to Gist if possible. In all other cases, return 'fg.'
       """
       if type (col) == IntType and (-10 <= col <= 63 ) \
          or self.legal_gist_colors.has_key (col) : # probably OK
          return col
       if self.narcisse_to_gist_col.has_key (col) :
          return self.narcisse_to_gist_col [col]
       print "Gist warning:" , col , "is an unknown color."
       return "fg"

   legal_gist_types = { "solid" : 1 , "dash" : 1 , "dot" : 1 ,
                        "dashdot" : 1 , "dashdotdot" : 1 , "none" : 1 }
   narcisse_to_gist_tp = { "line" : "solid" , "step" : "solid" ,
                           "unspecified" : "solid" , 6 : "dot" ,
                           "dotted" : "dot"}

   def _figure_type ( self , tp ) :
       """_figure_type ( tp ) does the best job it can to return
       a correct type. If the value is legal for Gist (even though
       it may mean something else in another system) then it is
       returned unchanged. If it is a Narcisse value, it is converted
       to Gist if possible. In all other cases, return 'solid.'
       """
       if type (tp) == IntType and (0 <= tp <= 5) \
          or self.legal_gist_types.has_key (tp) :
          return tp
       if self.narcisse_to_gist_tp.has_key ( tp ) :
          return self.narcisse_to_gist_tp [tp]
       print "Gist warning:" , tp , "is an unknown curve type."
       return "solid"

   _title_coords = [[.3935, .3655], [.3935, .8963], 
                    [.1278, .6483], [.6473, .6483]]
   _just = ["CB", "CT", "LH", "RH"]

   def _plot_titles ( self ) :
       """_plot_titles ( ) plots the four titles (if any)."""
       for i in range (4) :
          if self._titles [i] != '' and self._titles [i] != ' ' :
             plt (self._titles [i], self._title_coords [i][0],
                  self._title_coords [i][1],
                  color = self._figure_color (self._title_colors [i]),
                  justify = self._just [i]) ## , path = self._paths [i])

   def plot_object ( self , crv ) :

       """plot_object (crv) is a general purpose plotting routine. It should
       be called with one argument, a curve, quadmesh, polymap, cellarray,
       or lines object. In the case of multiple objects on one graph, the
       first call only should be to this routine, subsequent calls to
       add_object. plot_object does some one-time things, such as
       plotting the titles.
       """
       # Due to Gist problems calculating default axis limits correctly,
       # I'm doing it here.
       if crv.type () == QuadMeshType :
          if self._quadmesh != crv : # plmesh it
             self._quadmesh = crv
             # (ZCM 3/26/97) Do a plmesh first, otherwise boundary
             # gets all screwed up.
             plmesh (array (crv.y).astype (Float),
                     array (crv.x).astype (Float),
                     array (crv.ireg).astype (Int),
                     triangle = array (crv.tri).astype (Float))
          if self._xmindefault or self._xmaxdefault :
             new_limits = minmax ( crv.x )
             if self._xmindefault : self._xmin = new_limits [0]
             if self._xmaxdefault : self._xmax = new_limits [1]
          if self._ymindefault or self._ymaxdefault :
             new_limits = minmax ( crv.y )
             if self._ymindefault : self._ymin = new_limits [0]
             if self._ymaxdefault : self._ymax = new_limits [1]
       self._first_curve = 1
       n = current_window ( )
       if (n == -1 or n != self._n) :
          window (self._n)
       fma ( )
       self.add_object ( crv )

   def add_object ( self , crv ) :
       """add_object (crv) simply adds a curve to an existing graph.
       In Gist graphics, this is the same as plot_object except there
       is no frame advance."""
       global _color_bar_
       if crv.type () == QuadMeshType :
          if self._first_curve == 0 :
             # May need to compute new x and y axis limits
             if self._xmindefault or self._xmaxdefault :
                new_limits = minmax ( crv.x )
                if self._xmindefault and new_limits [0] < self._xmin :
                   self._xmin = new_limits [0]
                if self._xmaxdefault and new_limits [1] > self._xmax :
                   self._xmax = new_limits [1]
             if self._ymindefault or self._ymaxdefault :
                new_limits = minmax ( crv.y )
                if self._ymindefault and new_limits [0] < self._ymin :
                   self._ymin = new_limits [0]
                if self._ymaxdefault and new_limits [1] > self._ymax :
                   self._ymax = new_limits [1]
          if self._xyequal :
             # need to revise limits in light of this requirement
             xdiff = self._xmax - self._xmin
             ydiff = self._ymax - self._ymin
             if xdiff > ydiff :
                self._ymax = self._ymax + xdiff - ydiff
             elif ydiff > xdiff :
                self._xmax = self._xmax + ydiff - xdiff
       if self._xmin == 'd':
          self._xmin = 'e'
       if self._xmax == 'd':
          self._xmax = 'e'
       if self._ymin == 'd':
          self._ymin = 'e'
       if self._ymax == 'd':
          self._ymax = 'e'
       limits ( self._xmin, self._xmax )
       ylimits ( self._ymin, self._ymax)
       if hasattr (crv, "color") :
          crv.color = self._figure_color (crv.color)
       if hasattr (crv, "line_type") :
          crv.line_type = self._figure_type (crv.line_type)
       if crv.type () == CurveType :
          if self._first_curve:
             self._first_curve = 0
             self.send_color_card ( )
             gridxy (self._gridtype)
             logxy ( self._xscale, self._yscale )
             self._plot_titles ( )
          if (crv.marks == 0 or crv.marker is None) and crv.label != " " :
             marks = 1
             marker = crv.label [0]
          else :
             marks = crv.marks
             if crv.line_type == "dot" and crv.marker == "." \
                and crv.label != " " :
                marker = crv.label [0]
             else :
                marker = crv.marker
          if marks != 0 and (marker is None or marker == " ") :
              # Give user default markers
              plg ( array (crv.y).astype (Float),
                array (crv.x).astype (Float),
                hide = crv.hide, width = crv.width, type = crv.line_type,
                color = crv.color, marks = marks)
          else :
             plg ( array (crv.y).astype (Float),
                array (crv.x).astype (Float),
                hide = crv.hide, width = crv.width, type = crv.line_type,
                color = crv.color, marks = marks,
                marker = marker )
       elif crv.type () == LinesType :
          if self._first_curve:
             self._first_curve = 0
             self.send_color_card ( )
             gridxy (self._gridtype)
             logxy ( self._xscale, self._yscale )
             self._plot_titles ( )
          pldj ( array (crv.x0).astype (Float), array (crv.y0).astype (Float),
                 array (crv.x1).astype (Float), array (crv.y1).astype (Float),
                 hide = crv.hide, width = crv.width,
                 type = crv.line_type, color = crv.color )
       elif crv.type () == PolyMapType :
          if self._first_curve:
             self._first_curve = 0
             self.send_color_card ( )
             gridxy (self._gridtype)
             logxy ( self._xscale, self._yscale )
             self._plot_titles ( )
          if _color_bar_ :
             self.set_style ("z_work.gs")
             color_bar (min (ravel (crv.z)), max (ravel (crv.z)), 0)
          if type (crv.z[0]) == type (1.0) :
             plfp ( array (crv.z).astype (Float), array (crv.y).astype (Float),
                    array (crv.x).astype (Float), array (crv.n).astype (Int),
                    hide = crv.hide )
          else :
             plfp ( crv.z , array (crv.y).astype (Float),
                    array (crv.x).astype (Float), array (crv.n).astype (Int),
                    hide = crv.hide )
       elif crv.type () == CellArrayType :
          if self._first_curve:
             self._first_curve = 0
             self.send_color_card ( )
             gridxy (self._gridtype)
             logxy ( self._xscale, self._yscale )
          if _color_bar_ :
             self.set_style ("z_work.gs")
             color_bar (min (ravel (crv.z)), max (ravel (crv.z)), 0)
          if crv.x0 is None and crv.x1 is None :
             pli (crv.z, hide = crv.hide )
          elif crv.x0 is None :
             pli (crv.z, crv.x1, crv.y1,
                  hide = crv.hide )
          else :
             pli (crv.z, crv.x0, crv.y0, crv.x1, crv.y1,
                  hide = crv.hide )
          self._plot_titles ( )
       elif crv.type () == QuadMeshType :
          if self._first_curve:
             self._first_curve = 0
             self.send_color_card ( )
             gridxy (self._gridtype)
             logxy ( self._xscale, self._yscale )
          self._plot_quadmesh (crv)
          self._plot_titles ( )
       elif crv.type () == Animation2dType :
          if self._first_curve:
             self._first_curve = 0
             self.send_color_card ( )
             gridxy (self._gridtype)
             logxy ( self._xscale, self._yscale )
             self._plot_titles ( )
          self._do_animation ( crv ) 
       else :
          raise self._GistError , \
             "Unknown object type '" + crv.type () + "'."

   _BadZScale = "BadZScale"

   def _do_animation ( self, anim ) :
       """_do_animation ( anim ) calls the functions supplied in
       object anim in order to do an animation.
       """
       # anim initializes its variables
       anim.initialize ( anim )

       if anim.animation == 1 :
          animate (1)
       else :
          animate (0)

       fma ( )
       for i in range (anim.nsteps) :
           for j in range (len (anim.calculations)) :
               anim.calculations [j] (anim)
               plg (anim.y, anim.x, color = anim.color)
           fma ( )
           anim.update (anim)

       if anim.animation == 1 :
           animate (0)
           for j in range (len (anim.calculations)) :
               anim.calculations [j] (anim)
               plg (anim.y, anim.x, color = anim.color)
          
       
   def _compute_levels ( self , crv , z) :
       """_compute_levels ( crv ) computes the z coordinates of levels
       of contours, if necessary. if crv.levels is an integer, it
       interpolates that many values between zmin and zmax, and
       puts those in crv.levels.
       """
       if type (crv.levels) == IntType and crv.levels >= 2 :
          if crv.z_scale == "normal" :
             zlin = ravel (z)
             lzlin = len (zlin)
             zbar = add.reduce (zlin) / lzlin
             zs = sqrt ( (add.reduce (zlin ** 2) - lzlin * zbar ** 2) /
                 (lzlin - 1))
             z1 = zbar - 2. * zs
             z2 = zbar + 2. * zs
             diff = (z2 - z1) / (crv.levels - 1)
             levs = z1 + arange (crv.levels) * diff
          elif crv.z_scale == "lin" or crv.z_scale == "log" :
             try: 
                [z1, z2] = zmin_zmax (z, crv.ireg)
             except:
                [z1, z2] = true_minmax (z)
          else :
             raise self._BadZScale , \
                `crv.z_scale` + " is an invalid z scale specifier."
          if crv.z_scale == "lin" :
             levs = z1 + arange (1, crv.levels, typecode = Float) * \
                (z2 - z1) / (crv.levels + 1)
          elif crv.z_scale == "log" :
             levs = z1 + exp (arange (1, crv.levels, typecode = Float) * \
                log (z2 - z1) / (crv.levels + 1))
          crv.levels = levs

   _QuadMeshError_ = "QuadMeshError"
 
   def _plot_quadmesh ( self , crv ) :
       """_plot_quadmesh ( <mesh object> ) is called internally by
       plot_object/add_object to plot meshes. This is because plotting
       a mesh has so many more options; this will keep add_object
       from getting way too long.
       """
       global _color_bar_
       if hasattr (crv, "regions") and crv.regions is not None and \
          crv.regions != [] and crv.regions != "all" :
             self._plot_by_regions ( crv )
             return
       if hasattr (crv, "regions") and crv.regions != "all" and \
          crv.regions != [] :
          raise self._QuadMeshError_, "<" + `crv.regions` + \
             "> is not a valid value for the regions keyword."
       if crv.z is not None and crv.filled == 0 and crv.contours == 1 :
          # Must be a contour plot
          if crv.boundary == 1 :
             plm ( boundary = 1, type = self._figure_type (crv.boundary_type) ,
                   color = self._figure_color (crv.boundary_color))
          if crv.levels is None :
             # Accept 8 default levels
             if crv.marks == 0 or crv.marker is None and crv.label != " " :
                marks = 1
                marker = crv.label [0]
             else :
                marks = crv.marks
                marker = crv.marker
             plc (array (crv.z).astype (Float),
                  hide = crv.hide, type = crv.line_type,
                  width = crv.width, color = crv.color, marks = marks,
                  marker = marker, region = crv.region)
          else :
             # Accept user-specified levels
             self._compute_levels (crv, crv.z)
             if crv.marks == 0 or crv.marker is None and crv.label != " " :
                marks = 1
                marker = crv.label [0]
             else :
                marks = crv.marks
                marker = crv.marker
             plc (array (crv.z).astype (Float),
                  hide = crv.hide, type = crv.line_type,
                  width = crv.width, color = crv.color, marks = marks,
                  marker = marker, region = crv.region, levs = crv.levels)
       elif crv.filled == 1 :
          # Must be a filled mesh or filled contours
          if _color_bar_ :
             self.set_style ("z_work.gs")
          if crv.contours == 0 :
             # Plot mesh filled with background color (wire frame)
             # (crv.z is None) or colored by crv.z
             if crv.levels is None :
                if _color_bar_ :
                   color_bar (min (ravel (crv.z)), max (ravel (crv.z)), 0, ncol = 8)
                plfc (crv.z, crv.y, crv.x, crv.ireg,
                     contours = 8,
                     region = crv.region, scale = crv.z_scale)
             else :
                if _color_bar_ :
                   color_bar (min (ravel (crv.z)), max (ravel (crv.z)), 0,
                      ncol=crv.levels)
                plfc (crv.z, crv.y, crv.x, crv.ireg,
                     contours = crv.levels,
                     region = crv.region, scale = crv.z_scale)
             return
          if self._cmax != 0.0 or self._cmin != 0.0 :
             # Use byte scale
             (k, l) = shape (crv.z)
             _z_ = reshape (bytscl (reshape (crv.z, (k * l,)),
                cmin = self._cmin, cmax = self._cmax), (k, l)).astype('b')
          else : 
             _z_ = crv.z
          if _color_bar_ :
             self.set_style ("z_work.gs")
          if crv.contours == 1 :
             # eventually this will plot filled contours
             if crv.marks == 0 or crv.marker is None and crv.label != " " :
                marks = 1
                marker = crv.label [0]
             else :
                marks = crv.marks
                marker = crv.marker
             if crv.levels is None :
                if _color_bar_ :
                   color_bar (min (ravel (_z_)), max (ravel (_z_)), 0, ncol = 8)
                plfc (_z_.astype (Float), crv.y.astype (Float),
                   crv.x.astype (Float), crv.ireg.astype (Int),
                   contours = 8, region = crv.region, 
                   scale = crv.z_scale)
                plc (crv.z.astype (Float),
                   hide = crv.hide, type = crv.line_type,
                   width = crv.width, color = crv.color, marks = marks,
                   marker = marker, region = crv.region)
             else :
                if _color_bar_ :
                   color_bar (min (ravel (_z_)), max (ravel (_z_)), 0,
                      ncol = crv.levels)
                plfc (_z_.astype (Float), crv.y.astype (Float),
                   crv.x.astype (Float), crv.ireg.astype (Int),
                   contours = crv.levels, region = crv.region, 
                   scale = crv.z_scale)
                self._compute_levels (crv, crv.z)
                plc (crv.z.astype (Float),
                   hide = crv.hide, type = crv.line_type, levs = crv.levels,
                   width = crv.width, color = crv.color, marks = marks,
                   marker = marker, region = crv.region)
          self._plot_titles ( )
       else :
          # Must be a plain, ordinary, every day, garden-variety mesh
          # If boundary = 1, we need to do a boundary plot.
          # if either ktype or ltype (or both) is not "none",
          # we need to plot the k lines or l lines (or both).
          if crv.boundary == 1 :
             plm (
                  boundary = 1, hide = crv.hide, 
                  type = self._figure_type (crv.boundary_type),
                  color = self._figure_color (crv.boundary_color), 
                  width = crv.width)#, region = crv.region)
          if crv.ktype != "none" or crv.ltype != "none" :
             if crv.ktype == crv.ltype :
                plm (
                     inhibit = crv.inhibit,
                     hide = crv.hide, type = crv.ktype,
                     width = crv.width, color = crv.color)#, region = crv.region)
                return
             if crv.ktype != "none" :
                plm (
                     inhibit = 2,
                     hide = crv.hide, type = crv.ktype,
                     width = crv.width, color = crv.color)#, region = crv.region)
             if crv.ltype != "none" :
                plm (
                     inhibit = 1,
                     hide = crv.hide, type = crv.ltype,
                     width = crv.width, color = crv.color)#, region = crv.region)
       if crv.vx is not None :
          # Must be a vector field plot on top of everything else
          if crv.scale is not None :
             plv (array (crv.vy).astype (Float), array (crv.vx).astype (Float),
                  hide = crv.hide,
                  width = crv.width, color = crv.color,
                  region = crv.region, scale = crv.scale)
          else :
             plv (array (crv.vy).astype (Float), array (crv.vx).astype (Float),
                  hide = crv.hide,
                  width = crv.width, color = crv.color,
                  region = crv.region)
 
   def _plot_by_regions ( self , crv ) :
       """
       _plot_by_regions (crv) allows one to plot some subset of the
       regions of a mesh with various differing options.
       """
       plmesh (array (crv.y).astype (Float), array (crv.x).astype (Float),
               array (crv.ireg).astype (Int),
               triangle = array (crv.tri).astype (Float))
       # Loop through the regions, plotting each as specified:
       for reg in crv.regions :
          if crv.z is not None and reg.filled == 0 and reg.contours == 1 :
             # Must be a contour plot
             if reg.boundary == 1 :
                plm ( boundary = 1,
                      type = self._figure_type (reg.boundary_type) ,
                      color = self._figure_color (reg.boundary_color))
             if reg.levels is None :
                # Accept 8 default levels
                if reg.marks == 0 or reg.marker is None and reg.label != " " :
                   marks = 1
                   marker = reg.label [0]
                else :
                   marks = reg.marks
                   marker = reg.marker
                plc (array (crv.z).astype (Float), 
                     hide = reg.hide, type = reg.line_type,
                     width = reg.width, color = self._figure_color (reg.color),
                     marks = marks, marker = marker,
                     region = reg.number)
             else :
                # Accept user-specified levels
                if reg.marks == 0 or reg.marker is None and reg.label != " " :
                   marks = 1
                   marker = reg.label [0]
                else :
                   marks = reg.marks
                   marker = reg.marker
                self._compute_levels (reg, crv.z)
                plc (array (crv.z).astype (Float), 
                     hide = reg.hide, type = reg.line_type,
                     width = reg.width, color = self._figure_color (reg.color),
                     marks = marks, marker = marker,
                     region = reg.number, levs = reg.levels)
          elif reg.filled == 1 :
             # Must be a filled mesh
             if crv.z is None :
                # Plot mesh filled with background color (wire frame)
                plf (crv.z, 
                   hide = reg.hide,region = reg.number,edges = reg.edges,
                   ecolor = self._figure_color (reg.ecolor),
                   ewidth = reg.ewidth)
                return
             if self._cmax != 0.0 or self._cmin != 0.0 :
                # Use byte scale
                # crv.z needs to be reshaped to linear then shaped
                # back to 2d.
                (k, l) = shape (crv.z)
                _z_ = reshape (bytscl (reshape (crv.z,  (k * l, )),
                   cmin = self._cmin, cmax = self._cmax), (k, l)).astype('b')
             else : 
                _z_ = crv.z
             plf (array (_z_).astype (Float), 
                  hide = reg.hide,
                  region = reg.number, edges = reg.edges, ecolor = 
                  self._figure_color (reg.ecolor),
                  ewidth = reg.ewidth)
          else :
             # Must be a plain, ordinary, every day, garden-variety mesh
             if reg.boundary == 1 :
                plm (boundary = 1,
                     hide = reg.hide, type = 
                     self._figure_type (reg.boundary_type),
                     width = reg.width,
                     color = self._figure_color (reg.boundary_color),
                     region = reg.number)
             if reg.ktype != "none" or reg.ltype != "none" :
                if reg.ktype == reg.ltype :
                   plm (inhibit = reg.inhibit, 
                        type = self._figure_type (reg.ktype),
                        width = reg.width,
                        color = self._figure_color (reg.color),
                        region = reg.number)
                else :
                   if reg.ktype != "none" :
                      plm (inhibit = 2,  
                           type = self._figure_type (reg.ktype),
                           width = reg.width,
                           color = self._figure_color (reg.color),
                           region = reg.number)
                   if reg.ltype != "none" :
                      plm (inhibit = 1,  
                           type = self._figure_type (reg.ltype),
                           width = reg.width,
                           color = self._figure_color (reg.color),
                           region = reg.number)
          if crv.vx is not None and reg.vectors != 0 :
             # Must be a vector field plot on top of anything else
             if ( crv.scale is not None ) :
                plv (array (crv.vy).astype (Float),
                     array (crv.vx).astype (Float),
                     scale = crv.scale,
                     hide = reg.hide,
                     width = reg.width, color = self._figure_color (reg.color),
                     region = reg.number)
             else :
                plv (array (crv.vy).astype (Float),
                     array (crv.vx).astype (Float),
                     hide = reg.hide,
                     width = reg.width, color = self._figure_color (reg.color),
                     region = reg.number)
 
   def send_graph ( self, graf ) :
       if graf.type () == Graph3dType :
          lims = draw3 (1)
          if lims is not None:
             limits (lims [0], lims [1], lims [2], lims [3])
          if hasattr (self, "_hcp") and self._hcp != "" :
             hcp ()
 
   def set_yr_axis_limits ( self , v, *w ) :
       pass
 
   def synchronize ( self ) :
       pass
 
   def set_label_type ( self, t ) :
       pass
 
   def set_y_axis ( self , lr , *lm ) :
       pass
 
   def set_freeze_each ( self , n ) :
       pass

   def set_connect (self, val) :
       pass

   def set_link (self, val) :
       pass

   def set_mask (self, val) :
       pass

   def set_z_c_switch (self, val) :
       self._z_c_switch = val

   def set_z_contours (self, val) :
       pass

   def set_c_contours (self, val) :
       pass

   def send_generics (self, graf) :
       """send_generics ( graf ) sets up all the things that are generic to any
       graph. It does not actually do any plotting yet.
       """
       self.set_titles ( graf._titles )
       self.set_title_colors ( graf._title_colors )
       # The following line is redundant for Gist. Not sure about Narcisse.
       # self._plot_titles ( )
       if self._dims == 2 :
          self.set_grid_type ( graf._grid_type )
       self.clear_text ( )
       if is_scalar ( graf._text ) :
          if ( graf._text != "" and graf._text != " ") :
             self.set_text ( graf._text , 0 )
       else :
          for i in range ( len ( graf._text ) ) :
             self.set_text ( graf._text [i] , i )
       if is_scalar ( graf._text_color ) :
          self.set_text_color ( graf._text_color , 0 )
       else :
          for i in range ( len ( graf._text_color ) ) :
             self.set_text_color ( graf._text_color [i] , i )
       if is_scalar ( graf._text_size ) :
          self.set_text_size ( graf._text_size , 0 )
       else :
          for i in range ( len ( graf._text_size ) ) :
              self.set_text_size ( graf._text_size [i] , i )
       if is_scalar ( graf._tosys ) :
          self.set_tosys ( graf._tosys , 0 )
       else :
          for i in range ( len ( graf._tosys ) ) :
              self.set_tosys ( graf._tosys [i] , i )
       if is_scalar ( graf._text_pos ) :
          raise graf._GraphSpecError , \
             "Text position must be a point or an array of points."
       if len ( shape ( graf._text_pos )) == 1:
          self.set_text_pos ( graf._text_pos [0] , graf._text_pos [1] , 0 )
       else :
          for i in range (shape (graf._text_pos) [0] ) :
             self.set_text_pos ( graf._text_pos [i][0] ,
                               graf._text_pos [i][1] , i )
       if hasattr (graf, "n") and self._dims == 3 :
          if graf.mask_change :
             self.set_mask ( graf._s [graf.n - 1].mask )
       try:
          no_color = os.environ["NO_COLOR"]
       except KeyError:
          no_color = 0
       if no_color == 0 or no_color == "no" or no_color == "n" :
          self.set_color_card ( graf._color_card )

   def do_generic (self, graf) :
       self.send_generics ( graf )
       self.equal_axes = 0
       self.set_x_axis_limits (graf._axis_limits [0][0],
                               graf._axis_limits [0][1])
       self.set_y_axis_limits (graf._axis_limits [1][0],
                               graf._axis_limits [1][1])
       if graf._no_of_axes > 2 :
          self.set_z_axis_limits (graf._axis_limits [2][0],
                                  graf._axis_limits [2][1])
       if graf._no_of_axes > 3 :
          self.set_c_axis_limits (graf._axis_limits [3][0],
                                  graf._axis_limits [3][1])
       for i in range (graf._no_of_axes) :
          if graf._axis_scales [i] == "lin" :
             self.set_axis_lin (graf._axes [i])
          elif graf._axis_scales [i] == "log" :
             self.set_axis_log (graf._axes [i])
          elif graf._axis_scales [i] == "equal" :
             self.equal_axes = 1
          else :
             raise graf._AxisSpecError , \
                graf._axis_scales [i] + " is not a valid axis scale."
       return

   def quick_plot (self, graf) :
       "quick_plot (graf) plots without recomputing."
       ## Note: I am using phi and theta as the standard polar
       ## angles, so I have put in conversions to the Gist model.
       if graf.type () == Graph2dType :
          self._dims = 2
       else :
          self._dims = 3
       self.do_generic (graf)
       if self._dims == 3 :
          phi = graf._theta
          theta = graf._phi
          if phi is not None:
             phi = pi / 2. - phi  * pi / 180.
          if theta is not None:
             theta = pi / 2. - theta * pi / 180.
          orient3 (phi = phi, theta = theta)
          light3 (ambient = graf.lighting_dict ["ambient"],
                  diffuse = graf.lighting_dict ["diffuse"],
                  specular = graf.lighting_dict ["specular"],
                  spower = graf.lighting_dict ["spower"],
                  sdir = graf.lighting_dict ["sdir"])
          lims = draw3 (1)
          if lims is not None:
             limits (lims [0], lims [1], lims [2], lims [3])
          if hasattr (self, "_hcp") and self._hcp != "" :
             hcp ()
          self._plot_titles ()
       return

   def plot2d (self, graf) :
       """A Graph2d object calls plot2d with itself as argument.
       plot2d sorts out everything for the graph and then does the plot.
       The bulk of this work used to be done in Graph and Graph2d,
       but I decided it was too graphics-dependent.
       """
       global _color_bar_
       _color_bar_ = graf._color_bar
       self._dims = 2
       # (1) Do graph-generic stuff first
       self.do_generic (graf)
       # (2) Do the specifically 2d stuff
       sc = graf._axis_scales
       if graf._xyequal :
          self.set_xyequal ( )
       else :
          self.reset_xyequal ( )
       if is_scalar (sc) :
          if sc == "linlin" or sc == "lin" :
             self.set_linlin ( )
          elif sc == "linlog" :
             self.set_linlog ( )
          elif sc == "loglin" or sc == "log" :
             self.set_loglin ( )
          elif sc == "loglog" :
             self.set_loglog ( )
       else :
          if len (sc) == 1 :
             sc = sc + ["lin", "lin"]
          elif len (sc) == 2 :
             sc = sc + ["lin"]
          for n in range (3) :
             if sc [n] == "log" :
                self.set_axis_log (graf._axes [n])
             else :
                self.set_axis_lin (graf._axes [n])
       for i in range ( graf._c_ln ) :
          if i == 0 :
             self.plot_object ( graf._c [i] )
             if self.equal_axes :
                maxv = max (max (graf._c [i].x), max (graf._c [i].y))
                minv = min (min (graf._c [i].x), min (graf._c [i].y))
          else :
             self.add_object ( graf._c [i] )
             if self.equal_axes :
                maxv = max (max (graf._c [i].x), max (graf._c [i].y), maxv)
                minv = min (min (graf._c [i].x), min (graf._c [i].y), minv)
       if self.equal_axes :
          limits (minv, maxv, minv, maxv)
       self.plot_text ( )

   def plot3d (self, graf, draw = 1) :
       """plot3d (graf) plots a 3d graph object.
       If draw = 0, it only does the calculations, then returns.
       """
       ## Note: I am using phi and theta as the standard polar
       ## angles, so I have put in conversions to the Gist model.
       self._dims = 3
       # (1) Do graph-generic stuff first
       ## For the moment I am ignoring all the generic stuff
       ## except for text and titles. Axis scales and limits
       ## are gilding the lily.
       self.do_generic (graf)

       set_draw3_ (0) # do not draw prematurely
       # 3d graph generics:
       self.set_distance (graf._distance)
       phi = graf._theta
       theta = graf._phi
       if phi is not None:
          phi = - pi / 2. - phi  * pi / 180.
       if theta is not None:
          theta = pi / 2. - theta * pi / 180.
       orient3 (phi = phi, theta = theta)
       if graf._gnomon is not None :
          gnomon (graf._gnomon, chr = [graf._axis_labels [0] [0],
                                       graf._axis_labels [1] [0],
                                       graf._axis_labels [2] [0]])
       light3 (ambient = graf.lighting_dict ["ambient"],
               diffuse = graf.lighting_dict ["diffuse"],
               specular = graf.lighting_dict ["specular"],
               spower = graf.lighting_dict ["spower"],
               sdir = graf.lighting_dict ["sdir"])
       if graf._x_factor != 1. or graf._y_factor != 1. :
          limits_ (xfactor = graf._x_factor, yfactor = graf._y_factor)
       else :
          limits_ (square = graf._xyequal) # set defaults
       if graf._s_ln == 1 and graf._s [0].type () == SurfaceType :
          # A single surface will be plotted by plwf. The
          # allowed options are any one of "wm", "f3", "f4", "s3";
          # or "wm" together with either "f3" or "f4". An
          # intelligent attempt will be made to convert other
          # specifications to one of these.
          #
          # If surface has a _color_card, let it overrule any set
          # by the graph to which it belongs.
          if (hasattr (graf._s [0], "_color_card")) :
             self.set_color_card (graf._s [0]._color_card)
          s = graf._s [0]
          z = s.z
          y = s.y
          x = s.x
          c = s.c
          quad_mesh = 0
          if x is not None and y is not None :
             if len (x.shape) == 2 and len (y.shape) == 2 :
                # Two dimensional quadrilateral mesh (ZCM 7/2/97)
                if x.shape [0] != y.shape [0] or x.shape [1] != y.shape [1] :
                   raise self._GistError, \
                      "shape of x '" + `x.shape` + "' and shape of y '" + \
                      `y.shape` + "' do not match."
             elif len (x.shape) != 1 or len (y.shape) != 1 :
                raise self._GistError, "x and y must be one dimensional."
             elif len (x) != shape (z) [0] or len (y) != shape (z) [1] :
                raise self._GistError, \
                   "lengths of x and y must match shape of z."
             else : # Both are 1; expand to 2d
                y = multiply.outer (ones (len (x), Float), y)
                x = multiply.outer (x, ones (len (y), Float))
          opt = s.opt_3d
          if is_scalar (opt) : opt = [opt]
          if "wm" in opt or "w3" in opt or "w4" in opt :
             edges = 1
          else :
             edges = 0
          c_color_bar = 0
          z_color_bar = 0
          if s.z_c_switch :
             z_axis_max = self._cmax
             z_axis_min = self._cmin
             c_axis_max = self._zmax
             c_axis_min = self._zmin
             dum = z
             z = c
             c = dum
             scale = s.c_contours_scale
             if s.c_contours_array is not None :
                contours = s.c_contours_array
                ncol = len (contours + 1)
             elif s.number_of_c_contours is not None :
                contours = s.number_of_c_contours
                ncol = contours + 1
             else :
                contours = 20
                ncol = contours + 1
          else :
             z_axis_max = self._zmax
             z_axis_min = self._zmin
             c_axis_max = self._cmax
             c_axis_min = self._cmin
             scale = s.z_contours_scale
             if s.z_contours_array is not None :
                contours = s.z_contours_array
                ncol = len (contours + 1)
             elif s.number_of_z_contours is not None :
                contours = s.number_of_z_contours
                ncol = contours + 1
             else :
                contours = 20
                ncol = contours + 1
          if ("f3" in opt or "s3" in opt or "w3" in opt or "i3" in opt) :
             fill = z
             z_color_bar = graf._color_bar
          elif ("f4" in opt or "s4" in opt or "w4" in opt or "i4" in opt) :
             fill = c
             c_color_bar = graf._color_bar
          else :
             fill = None
          if s.mask == "none" :
             cull = 0
          else :
             cull = 1
          if c_color_bar or z_color_bar :
             color_bar (min (ravel (fill)), max (ravel (fill)), graf._split, ncol)
          if "s3" in opt or "i3" in opt :
             # Plot contours in z (or c) direction
             if x is None and y is None :
                x = arange (fill.shape [0], typecode = Float)
                y = arange (fill.shape [1], typecode = Float)
             [nv, xyzv, dum] = slice3mesh (x, y, fill)
             if "s3" in opt :
                plzcont (nv, xyzv, contours = contours, scale = scale,
                   edges = edges, zaxis_min = z_axis_min, zaxis_max = z_axis_max)
             else :
                plzcont (nv, xyzv, contours = contours, scale = scale,
                   edges = edges, color = "bg", zaxis_min = z_axis_min,
                   zaxis_max = z_axis_max)
          elif "s4" in opt or "i4" in opt :
             # Plot contours in c (or z) direction
             if x is None and y is None :
                x = arange (fill.shape [0], typecode = Float)
                y = arange (fill.shape [1], typecode = Float)
             [nv, xyzv, col] = slice3mesh (x, y, z, color = fill, smooth = 1)
             if "s4" in opt :
                pl4cont (nv, xyzv, col, contours = contours, scale = scale,
                   edges = edges, caxis_min = c_axis_min, caxis_max = c_axis_max)
             else :
                pl4cont (nv, xyzv, col, contours = contours, scale = scale,
                   edges = edges, color = "bg", caxis_min = c_axis_min,
                   caxis_max = c_axis_max)
          else :
             plwf (z, y, x, fill = fill, shade = s.shade, edges = edges,
                ecolor = s.ecolor, ewidth = 1, cull = cull, scale =
                graf.z_scale)
          if c_color_bar or z_color_bar :
             self.set_style ("z_nobox.gs")
          else :
             self.set_style (graf._style)
          if draw == 0 :
             self.set_color_card (graf._color_card, now = 1)
             return
          self.lims = draw3 (1)
          if hasattr (self, "_hcp") and self._hcp != "" :
             hcp ()
          limits (self.lims [0], self.lims [1], self.lims [2], self.lims [3])
          self.send_color_card ()
       else :
          # Either there is more than one surface, or else at
          # least one is a mesh3d or a slice, so in any event must
          # be handled differently.  We have to accumulate the pieces
          # of the graph in a pl3tree.
          # N. B. I haven't implemented the Narcisse-style mesh3d
          # plot at this time. It is unclear to me how useful this is.
          clear3 ()
          isosurfaces_present = 0
          # If a split palette is called for, compute separate maxima and
          # minima for the whole graph (used to color plane sections)
          # and for isosurfaces (to effect maximum contrast in shading).
          zmax = None
          zmin = None
          for i in range (graf._s_ln) :
             s = graf._s [i]
             #
             # If some object has a _color_card, let it overrule any set
             # by the graph to which it belongs (last occurring will be used).
             if (hasattr (graf._s [i], "_color_card")) :
                self.set_color_card (graf._s [i]._color_card)
             if s.type () == Slice3dType :
                if s.plane is None :
                   isosurfaces_present = 1
                if s.val is not None : # Get max and min of values
                   if zmax is None :
                      zmax = max (s.val)
                      zmin = min (s.val)
                   else :
                      zmax = max (zmax, max (s.val))
                      zmin = min (zmin, min (s.val))
             elif s.type () == SurfaceType :
                isosurfaces_present = 1
          self.send_color_card ()
          if isosurfaces_present and graf._split != 0 :
             if self._current_palette is None:
                split_palette ( )
             else:
                split_palette (self._current_palette)
             spl = 1
          else:
             spl = 0
          for i in range (graf._s_ln) :
             ncol = None
             if i == 0 :
                clear3 ()
             s = graf._s [i]
             if s.type () == SurfaceType :
                z = s.z
                y = s.y
                x = s.x
                c = s.c
                opt = s.opt_3d
                if is_scalar (opt) : opt = [opt]
                if "wm" in opt :
                   edges = 1
                else :
                   edges = 0
                if s.z_c_switch :
                   z_axis_max = self._cmax
                   z_axis_min = self._cmin
                   c_axis_max = self._zmax
                   c_axis_min = self._zmin
                   dum = z
                   z = c
                   c = dum
                   scale = s.c_contours_scale
                   if s.c_contours_array is not None :
                      contours = s.c_contours_array
                   elif s.number_of_c_contours is not None :
                      contours = s.number_of_c_contours
                   else :
                      contours = 20
                else :
                   z_axis_max = self._zmax
                   z_axis_min = self._zmin
                   c_axis_max = self._cmax
                   c_axis_min = self._cmin
                   scale = s.z_contours_scale
                   if s.z_contours_array is not None :
                      contours = s.z_contours_array
                   elif s.number_of_z_contours is not None :
                      contours = s.number_of_z_contours
                   else :
                      contours = 20
                if ("f3" in opt or "s3" in opt or "w3" in opt) :
                   fill = z
                   z_color_bar = graf._color_bar
                elif ("f4" in opt or "s4" in opt or "w4" in opt) :
                   fill = c
                   c_color_bar = graf._color_bar
                else :
                   fill = None
                if "s3" in opt or "i3" in opt :
                   # Plot contours in z (or c) direction
                   if x is None and y is None :
                      x = arange (fill.shape [0], typecode = Float)
                      y = arange (fill.shape [1], typecode = Float)
                   [nv, xyzv, dum] = slice3mesh (x, y, fill)
                   if "s3" in opt :
                      plzcont (nv, xyzv, contours = contours, scale = scale,
                         clear = 0, edges = edges, zaxis_min = z_axis_min,
                         cmin = zmin, cmax = zmax, zaxis_max = z_axis_max)
                         # Do not clear multiple curves
                   else :
                      plzcont (nv, xyzv, contours = contours, scale = scale,
                         clear = 0, edges = edges, color = "bg", zaxis_min = z_axis_min,
                         cmin = zmin, cmax = zmax, 
                         zaxis_max = z_axis_max)
                elif "s4" in opt or "i4" in opt :
                   # Plot contours in c (or z) direction
                   if x is None and y is None :
                      x = arange (fill.shape [0], typecode = Float)
                      y = arange (fill.shape [1], typecode = Float)
                   [nv, xyzv, col] = slice3mesh (x, y, z, color = fill,
                       smooth = 1)
                   if "s4" in opt :
                      pl4cont (nv, xyzv, col, contours = contours,
                         scale = scale, clear = 0, edges = edges,
                         cmin = zmin, cmax = zmax,
                         caxis_min = c_axis_min, caxis_max = c_axis_max)
                   else :
                      pl4cont (nv, xyzv, col, contours = contours,
                         scale = scale, clear = 0, edges = edges, color = "bg",
                         cmin = zmin, cmax = zmax,
                         caxis_min = c_axis_min, caxis_max = c_axis_max)
                       # Do not clear multiple curves
                else :
                   if not hasattr (s, "s3") :
                      if x is None :
                         [nv, xyzv, val] = slice3mesh (z, color = c)
                      else :
                         [nv, xyzv, val] = slice3mesh (x, y, z, color = c)
                      setattr (s, "s3", Slice (nv, xyzv, val))
                   if zmax is None :
                      zmax = max (ravel (s.s3.val))
                      zmin = min (ravel (s.s3.val))
                   else :
                      zmax = max (zmax, max (ravel (s.s3.val)))
                      zmin = min (zmin, min (ravel (s.s3.val)))
                   pl3tree (s.s3.nv, s.s3.xyzv, s.s3.val, cmin = zmin, cmax = zmax,
                      split = spl, edges = edges)
             elif s.type () == Mesh3dType :
                raise self._GistError, "Narcisse style mesh plot not implemented."
             elif s.type () == Slice3dType :
                opt = s.opt_3d
                if type (opt) != ListType : opt = [opt]
                if "wm" in opt :
                   edges = 1
                else :
                   edges = 0
                # fval is what we plot if we ask for "flat" mode
                if "f3" in opt :
                   fval = s.xyzv[:,2]
                elif "f4" in opt :
                   fval = s.val
                if "s3" in opt or "i3" in opt or "s4" in opt or "i4" in opt :
                   if s.contours is None:
                      contours = 10
                      ncol = contours + 1
                   else :
                      contours = s.contours
                      if type (contours) == IntType :
                         ncol = contours + 1
                      else :
                         ncol = len (contours)
                   if "s3" in opt :
                      plzcont (s.nv, s.xyzv, contours = contours,
                         scale = s.scale, clear = 0, edges = s.edges,
                         cmin = zmin, cmax = zmax, zaxis_max = self._zmax,
                         zaxis_min = self._zmin, split = spl)
                   elif "i3" in opt :
                      plzcont (s.nv, s.xyzv, contours = contours,
                         scale = s.scale, clear = 0, edges = s.edges,
                         color = "bg", cmin = zmin, cmax = zmax,
                         split = spl)
                   elif "s4" in opt :
                      pl4cont (s.nv, s.xyzv, s.val, contours = contours,
                         scale = s.scale, clear = 0, edges = s.edges,
                         cmin = zmin, cmax = zmax,
                         caxis_min = self._cmin, caxis_max = self._cmax,
                         split = spl)
                   else :
                      pl4cont (s.nv, s.xyzv, s.val, contours = contours,
                         scale = s.scale, clear = 0, edges = s.edges,
                         cmin = zmin, cmax = zmax,
                         color = "bg", caxis_min = self._cmin, caxis_max = self._cmax,
                         split = spl)
                elif ("f3" in opt or "f4" in opt) and fval is not None \
                   and type (fval) == ArrayType :
                   if graf._split != 0 and isosurfaces_present :
                      if s.plane is not None :
                         pl3tree (s.nv, s.xyzv,
                            split_bytscl (fval, upper = 0, cmin = zmin,
                               cmax = zmax).astype('b'), s.plane,
                               cmin = zmin, cmax = zmax, split = 0,
                               edges = edges)
                      else :
                         pl3tree (s.nv, s.xyzv,
                            split_bytscl (fval, upper = 1, cmin = zmin,
                               cmax = zmax).astype('b'), s.plane,
                               cmin = zmin, cmax = zmax, split = 0,
                               edges = edges)
                   else :
                      pl3tree (s.nv, s.xyzv, bytscl (fval, top = 199,
                         cmin = zmin, cmax = zmax).astype('b'), s.plane,
                               cmin = zmin, cmax = zmax, split = 0,
                         edges = edges)
                else :
                   pl3tree (s.nv, s.xyzv, s.val, s.plane,
                         cmin = zmin, cmax = zmax, split = graf._split,
                         edges = edges)
             else :
                raise self._GistError, \
                   "object for 3d plot must be a Surface, a Mesh3d, " + \
                   "or one or more Slices."
          if zmax is not None and graf._color_bar != 0 :
             self.set_style ("z_nobox.gs")
          else :
             self.set_style (graf._style)
          if draw == 0 :
             return
          self.lims = draw3 (1)
          if hasattr (self, "_hcp") and self._hcp != "" :
             hcp ()
          if self.lims is not None :
             limits (self.lims [0], self.lims [1],
                self.lims [2], self.lims [3])
          if zmax is not None and zmin != zmax and graf._color_bar != 0 :
             color_bar (zmin, zmax, spl, ncol)

       # text must be done last or it might be obscured.
       self.plot_text ( )
       self._plot_titles ( )

   def move_light (self, i) :
       if i >= self.nframes : return 0
       theta = pi / 4. + (i - 1) * self.angle
       light3 (sdir = array ( [cos (theta), .25, sin(theta)], Float))
       draw3 (not self.making_movie)
       if hasattr (self, "_hcp") and self._hcp != "" :
          hcp ()
       return 1

   def move_light_source (self, graf, angle, nframes) :
       self.angle = angle
       self.nframes = nframes
       # Have to call plot3d to do calculations.
       self.plot3d (graf, draw = 0)
       limits (self.lims [0], self.lims [1], self.lims [2], self.lims [3])
       self.making_movie = 1
       set_draw3_ (0)
       movie (self.move_light, lims = self.lims)
       self.making_movie = 0
       fma ()
       self.move_light (1)

   def rotate_graph (self, graf, axis, angle, nframes) :
       # Have to call plot3d to do calculations.
       self.plot3d (graf, draw = 0)
       [lims0, lims1, lims2, lims3] = draw3 (1)
       limits (lims0, lims1, lims2, lims3)
       if hasattr (self, "_hcp") and self._hcp != "" :
          hcp ()
       spin3 (axis = axis, angle = angle * nframes,
          nframes = nframes) # , lims = [lims0, lims1, lims2, lims3])
       [lims0, lims1, lims2, lims3] = draw3 (1)
       limits (lims0, lims1, lims2, lims3)
       if hasattr (self, "_hcp") and self._hcp != "" :
          hcp ()
