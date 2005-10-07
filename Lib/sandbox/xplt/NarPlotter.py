# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

import narcisse
from Numeric import *
from scipy_base.fastumath import *
# We need types to check args to some routines
from types import *
from graftypes import *
from shapetest import *
from arrayfns import *
from string import uppercase
import os

def minmax1 ( x ) :
   """minmax1 (x) where x is a one-dimensional array computes the minimum and
   maximum values in the array and returns them as a list [min, max].
   """
   max = x [0]
   min = x [0]
   for i in range (len (x)) :
      if x [i] > max : max = x[i]
      if x [i] < min : min = x [i]
   return [floor (min), ceil (max)]


def minmax2 ( x ) :
   """minmax2 (x) where x is a two-dimensional array computes the minimum and
   maximum values in the array and returns them as a list [min, max]. I use
   this routine because there are apparently some circumstances in which Gist
   fails to calculate default axis limits correctly.
   """
   max = x [0, 0]
   min = x [0, 0]
   for i in range (shape (x) [0]) :
      for j in range (shape (x) [1]) :
         if x [i, j] > max : max = x [i, j]
         if x [i, j] < min : min = x [i, j]
   return [floor (min), ceil (max)]

NarFloat = 'f'
NarInt = 'i'

# Define a simple title function: missing arguments become blanks.
class Plotter :

   def type (self) :
       return NarType

   def open ( self , filename = ' ' ) :
       """open ( string ) opens a connection to Narcisse (if it can)
       using filename 'string.'
       """
       if self._file_open :
          if self._file_name == filename :
             return # quietly
          else :
             raise self.ConnectException , \
                 "This instance already open with filename '" + \
                 self._file_name + "'."
       else :
          self._file_descr = narcisse.naropen ( filename )
          if self._file_descr >= 0 :
             self._file_open = 1
             self._file_name = filename
          else :
             raise self.ConnectException , \
                 "Unable to open graphics file '" + filename + "'."

   _cgm_warning = 0
   _ps_warning = 0

   def __init__ ( self , filename = ' ' , ** kw ) :
       self.NarError = "NarError"
       if filename == "none" :
          if not self._cgm_warning :
             print "Sorry, Narcisse does not write cgm files."
             print "...This will be your only warning."
             self._cgm_warning = 1
       elif len (filename) >= 3 and filename [-3:] == ".ps" :
          if not self._ps_warning :
             print "Sorry, Narcisse does not write postscript files"
             print "except from the graphical user interface."
             print "...This will be your only warning."
             self._ps_warning = 1
       self._file_open = 0
       self._frozen = 0
       self._freeze_each = 0
       self._mono = 0    #defaults to color
       self._file_descr = -1
       self.ConnectException = "ConnectException"
       self.open ( filename )
       self.freeze_graph ( )
       self.set_grid_type ( "axes" )
       self._xyequal = 0
       self.set_default_axes_limits () # let Narcisse determine limits
       self.set_axis_lin ("all")       # all axes linear scales
       narcisse.narsetar ( "curve_label_x_min", 0.2 )
       narcisse.narsetar ( "curve_label_x_max", 0.2 )
       narcisse.narsetar ( "curve_label_y_min", 0.2 )
       narcisse.narsetar ( "curve_label_y_max", 0.2 )
       narcisse.narsetvals ( self._file_descr )
       self._x_axis_min = 0.
       self._y_axis_min = 0.
       self._yr_axis_min = 0.
       self._z_axis_min = 0.
       self._c_axis_min = 0.
       self._x_axis_max = 0.
       self._y_axis_max = 0.
       self._yr_axis_max = 0.
       self._z_axis_max = 0.
       self._c_axis_max = 0.
       self.clear_text ( )
       self.set_text_color (2, 0) #black or nearly so
       self.set_axis_labels ()         # To English defaults
       self.set_titles ( )
       self.set_title_colors ( )
       self.plot_curve = self.plot_object
       self.add_curve = self.add_object
       self._graph_type = 0
       if kw.has_key ("style") :
          self._style = kw ["style"]
       else :
          self._style = " "
       self._next_letter = 0

   def close ( self ) :
       "close () closes the connection to Narcisse."
       if self._file_open :
          narcisse.narclose ( self._file_descr )
          self._file_descr = -1
          self._file_open = 0
          self._file_name = ""

   def __del__ ( self ) :
       self.close ( )

   def new_frame (self) :
       return

   def set_tosys (self, *x) :
       return

   def set_mono ( self ) :
       """set_mono () will set the 3d display mode permanently to
       monochrome mesh. This is the only meaningful display
       mode if you are only displaying 3d data on a monochrome
       monitor. Calls to set_3d_options will do nothing
       (silently). Call set_color () to allow color options
       again."""

       self.set_3d_options ( color_bar, color_bar_pos, "wm" )
       self._mono = 1

   def synchronize ( self ) :
       if self._file_open :
          narcisse.narsync ( self._file_descr )
       else :
          print "synchronize: sorry, nothing is open to synchronize with."

   def query ( self ) :
       if not self._file_open :
          return -1
       else :
          return narcisse.narquery ( self._file_name )
  
   def set_color ( self ) :
       """set_color ( ) will allow you to use the color 3d options
       which are disabled by set_mono ( )."""

       self._mono = 0

   # Everything on a 2d graph shares the same color card:
   # (This dictionary is used to convert Narcisse color card names
   # to the numbers required by the plotting routines. It will also
   # convert Gist names.)
   narcisse_color_card_dict = { "absolute" : 0 , "binary" : 1 ,
   "bluegreen" : 2 , "default" : 6 , "negative" : 4 , "positive" : 5 ,
   "rainbow" : 6 , "rainbowhls" : 7 , "random" : 8 , "redblue" : 9 ,
   "redgreen" : 10 , "shifted" : 11 ,"earth.gp" : 8 , "stern.gp" : 2 ,
   "rainbow.gp" : 7 , "heat.gp" : 10 , "gray.gp" : 0 , "yarg.gp" : 4 }

   def set_color_card ( self , h , now = 0) :
       """set_color_card ( n ) indicates a predefined color card
       for a plot. See the manual for the values of n and the
       color card selected (sec. 4.2.134, parametre_map)."""

       if self.narcisse_color_card_dict.has_key (h) :
          h = self.narcisse_color_card_dict [h]
       narcisse.narsetai ("parameter_map", h)
       narcisse.narsetvals (self._file_descr)
  
   def set_titles ( self , * vals ):
       """set_titles ('bottom', 'top', 'left', 'right')
         All arguments are optional. Missing ones default to ' '."""

       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if len (vals) == 0 :
          vals = []
       elif type (vals [0]) == StringType :
          vals = [vals [0]]
       else :
          vals = vals [0]
       if len (vals) == 0 :
          vals = [ " " , " " , " " , " " ]
       elif len (vals) == 1 :
          vals = vals + [ " " , " " , " " ]
       elif len (vals) == 2 :
          vals = vals + [ " " , " " ]
       elif len (vals) == 3 :
          vals = vals + [ " " ]
       elif len (vals) <> 4 :
          print "titles must be one string or a list of up to four strings!"
          return
       narcisse.narsetac ( "title_value_bottom" , vals [0] )
       narcisse.narsetac ( "title_value_top" , vals [1] )
       narcisse.narsetac ( "title_value_left" , vals [2] )
       narcisse.narsetac ( "title_value_right" , vals [3] )
       narcisse.narsetvals ( self._file_descr )

   # Translation table from color names to Narcisse (only works for rainbowhls)
   gist_to_narcisse_col = { "bg" : 0, "background" : 0, "fg" : 1,
                            "foreground" : 1, "blue" : 2, "green" : 3,
                            "yellow" : 4 , "orange" : 5 , "red" : 6,
                            "magenta" : 7, "purple" : 7, "black" : 8,
                            "white" : 9, "cyan" : 20 , "yellowgreen" : 39,
                            "gold" : 42 , "orangered" : 47, "redorange" : 48,
                            -1 : 0 , -2 : 1 , -3 : 8 , -4 : 9 , -5 : 6 ,
                            -6 : 3 , -7 : 2 , -8 : 20 , -9 : 7 , -10 : 4 }

   def _figure_color (self , col) :
       """_figure_color ( col ) does the best job it can to return
       a correct color. If the value is legal for Narcisse (even though
       it may mean something else in another system) then it is
       returned unchanged. If it is a Gist value, it is converted
       to Narcisse if possible. In all other cases, return 1.
       """
       if type (col) == IntType and 0 <= col <= 63 :
          return col
       if self.gist_to_narcisse_col.has_key (col) :
          return self.gist_to_narcisse_col [col]
       return 1

   def set_title_colors ( self , * vals ) :
       """set_title_colors (bottom_color, top_color, left_color, right_color)
         All arguments are optional, integers from 0 to 63 representing
         a color in some color map. Missing arguments default
         to foreground."""
 
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if len (vals) == 0 :
          vals = []
       elif type (vals [0]) == IntType :
          vals = [vals [0]]
       else :
          vals = vals [0]
       if len (vals) == 0 :
          vals = [ 1 , 1 , 1 , 1 ]
       elif len (vals) == 1 :
          vals = [vals [0]] + [ 1 , 1 , 1 ]
       elif len (vals) == 2 :
          vals = [vals [0]] + [vals [1]] + [ 1 , 1 ]
       elif len (vals) == 3 :
          vals = [vals [0]] + [vals [1]] + [vals [2]] + [ 1 ]
       elif len (vals) <> 4 :
          raise self.NarError ,\
             "Title color must be list of size 4 or less."
       else :
          vals = [vals [0]] + [vals [1]] + [vals [2]] + [vals [3]]
       for i in range (4) :
          vals[i] = self._figure_color (vals [i])
       narcisse.narsetai ( "title_color_bottom" , vals [0] )
       narcisse.narsetai ( "title_color_top" , vals [1] )
       narcisse.narsetai ( "title_color_left" , vals [2] )
       narcisse.narsetai ( "title_color_right" , vals [3] )
       narcisse.narsetvals ( self._file_descr )
 
   def set_grid_type ( self , * val ) :
       """set_grid_type ( string ) determines how intrusive the axes
       and grids are. The legal arguments are:
       'none'--no axes and grids are drawn.
       'axes'--axes with tick marks.
       'wide'--widely spaced grid in x and y (2d or 3d).
       'full'--narrowly spaced grid in x and y (2d or 3d).
       If no argument is specified, the default is 'axes'."""
 
       if len ( val ) > 1 :
          raise self.NarError , "Too many arguments to set_grid_type."
       if len ( val ) == 0 or val [0] == "axes" :
          narcisse.narsetai ( "grid_type" , 1 )
       elif val [0] == "none" :
          narcisse.narsetai ( "grid_type" , 0 )
       elif val [0] == "wide" :
          narcisse.narsetai ( "grid_type" , 2 )
       elif val [0] == "full" :
          narcisse.narsetai ( "grid_type" , 3 )
       else :
          raise self.NarError , val [0] + \
             " is an inappropriate argument for set_grid_type."
       narcisse.narsetvals ( self._file_descr )
 
   def set_3d_grid_type ( self , val ) :
       """set_3d_grid_type (gt) sets what the wire grid will look like
       in a 3d surface plot in one of the wire modes. The choices
       for gt are 'x' (x lines only), 'y' (y lines only) and 'xy'
       (both x and y lines)."""
 
       if val == "x" :
          narcisse.narsetai ( "option_3d_grid_type" , 0 )
       elif val == "y" :
          narcisse.narsetai ( "option_3d_grid_type" , 1 )
       else : # You'll get "xy" if you goof.
          narcisse.narsetai ( "option_3d_grid_type" , 2 )
       narcisse.narsetvals ( self._file_descr )
 
   def set_connect ( self , val ) :
       """set_connect (cn) tells whether to connect two or more
       surface plots, which presumably improves masking.
       cn=1 to connect, cn=0 to disconnect."""
 
       narcisse.narsetai ("option_3d_conv_mode" , val)
       narcisse.narsetvals ( self._file_descr )
 
   def set_link ( self , val ) :
       """set_link (ln) tells whether to link two or more surfaces
       plotted with different 3d options into one plot (otherwise
       all surfaces will have the same options). ln=1 to link,
       ln = 0 not to link. This needs to be set to 1 for all surfaces
       except the last. Connection must not be set (see set_connect ()).
       The axes must not be plotted for surfaces after the first."""
 
       narcisse.narsetai ("parameter_scene", val)
       narcisse.narsetvals ( self._file_descr )
 
   def set_z_c_switch ( self , val ) :
       """set_z_c_switch (sw) tells whether to switch the roles
       of the z and c variables in a 4d plot. sw=1 to do the
       switch, sw=0 not to do it."""
 
       narcisse.narsetai ("option_3d_z_or_c", val)
       narcisse.narsetvals ( self._file_descr )
 
   # routine to label the axes
   def set_axis_labels ( self , * vals ):
       """set_axis_labels ('x_label', 'y_label', 'z_label', 'yr_label')
         All arguments are optional. Default values (from right):
         ' ', 'Z axis', 'Y axis', 'X axis'."""
 
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if len (vals) == 1 and (type (vals [0]) == TupleType or
                               type (vals [0]) == ListType) :
          valsin = vals[0]
       else :
          valsin = vals
       vals = [ "X axis" , "Y axis" , "Z axis" , " " ]
       if len (valsin) >= 1 :
          vals [0] = valsin [0]
       if len (valsin) >= 2 :
          vals [1] = valsin [1]
       if len (valsin) >= 3 :
          vals [2] = valsin [2]
       narcisse.narsetac ( "x_axis_title" , vals [0] )
       narcisse.narsetac ( "y_axis_title" , vals [1] )
       narcisse.narsetac ( "z_axis_title" , vals [2] )
       narcisse.narsetac ( "yr_axis_title" , vals [3] )
       narcisse.narsetvals ( self._file_descr )
 
   # routines to set axis scales -- linear scales
   def set_axis_lin ( self , ax ) :
       """set_axis_lin (ax) where ax can be 'x', 'y', 'yr', 'z', 'c', or 'all'.
          The specified axis will have a linear scale."""
 
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if ( ax == "x" ) :
          narcisse.narsetai ( "x_axis_log" , 0 )
       elif ( ax == "y" ) :
          narcisse.narsetai ( "y_axis_log" , 0 )
       elif ( ax == "yr" ) :
          narcisse.narsetai ( "yr_axis_log" , 0 )
       elif ( ax == "z" ) :
          narcisse.narsetai ( "z_axis_log" , 0 )
       elif ( ax == "c" ) :
          narcisse.narsetai ( "c_axis_log" , 0 )
       elif ax == "all" :
          narcisse.narsetai ( "x_axis_log" , 0 )
          narcisse.narsetai ( "y_axis_log" , 0 )
          narcisse.narsetai ( "yr_axis_log" , 0 )
          narcisse.narsetai ( "z_axis_log" , 0 )
          narcisse.narsetai ( "c_axis_log" , 0 )
       else :
          raise self.NarError , "set_axis_lin: axis must be x, y, yr, z, or c."
       narcisse.narsetvals ( self._file_descr )
 
   # routines to set axis scales -- log scales
   def set_axis_log ( self , ax ) :
       """set_axis_log (ax) where ax can be 'x', 'y', 'yr', 'z', 'c', or 'all'.
          The specified axis will have a logarithmic scale."""
 
       if ( ax == "x" ) :
          narcisse.narsetai ( "x_axis_log" , 1 )
       elif ( ax == "y" ) :
          narcisse.narsetai ( "y_axis_log" , 1 )
       elif ( ax == "yr" ) :
          narcisse.narsetai ( "yr_axis_log" , 1 )
       elif ( ax == "z" ) :
          narcisse.narsetai ( "z_axis_log" , 1 )
       elif ( ax == "c" ) :
          narcisse.narsetai ( "c_axis_log" , 1 )
       elif ax == "all" :
          narcisse.narsetai ( "x_axis_log" , 1 )
          narcisse.narsetai ( "y_axis_log" , 1 )
          narcisse.narsetai ( "yr_axis_log" , 1 )
          narcisse.narsetai ( "z_axis_log" , 1 )
          narcisse.narsetai ( "c_axis_log" , 1 )
       else :
          raise self.NarError , "axis_log: axis must be x, y, yr, z, or c."
       narcisse.narsetvals ( self._file_descr )
 
   # special routines to set both x and y scales at once
   def set_linlin ( self ) :
       "set_linlin () sets both x and y axes to linear scale."
       self.set_axis_lin ( "x" )
       self.set_axis_lin ( "y" )
 
   def set_linlog ( self ) :
       'set_linlog () sets x axis to linear, y axis to logarithmic.'
       self.set_axis_lin ( "x" )
       self.set_axis_log ( "y" )
 
   def set_loglin ( self ) :
       'set_loglin () sets x axis to logarithmic, y axis to linear.'
       self.set_axis_log ( "x" )
       self.set_axis_lin ( "y" )
 
   def set_loglog ( self ) :
       'set_loglog () sets both x and y axes to logarithmic scale.'
       self.set_axis_log ( "x" )
       self.set_axis_log ( "y" )
 
   #determine which y axis to use for a curve
   def set_y_axis ( self , val1 , * val2 ) :
       """use set_y_axis ( 'left' , n ) or set_y_axis ( 'right' , n )
       to cause curve number n to be associated with the left or
       right y axis."""
 
       if len ( val2 ) == 2 :
          n = val2 [1]
       else :
          n = 0       # set for curve 0 if not specified
       if len ( val2 ) == 0 or val2 [0] == "left" or val2 [0] != "right" :
          narcisse.narsetaii ( "curve_y_axis" , 0 , n )
       else :
          narcisse.narsetaii ( "curve_y_axis" , 1 , n )
       narcisse.narsetvals ( self._file_descr )

   def set_bytscl ( self, cmin, cmax ) :
       return

   def add_text (self, str, x, y, size, color="fg", tosys = 1) :
       """add_text (str, x, y, size [, color]) adds a text to a graph."""
       return
 
   # set the maximum value of an axis
   def set_axis_max ( self , ax , * val1 ) :
       """set_axis_max (ax, val) where ax is 'x', 'y', 'z', 'yr', or 'c'.
          The maximum of the specified axis will be set to val.
          val should be a PyFloat object."""
 
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if len ( val1 ) == 0 :
          val = 0.0
       else :
          val = val1 [0]
       if ( ax == "x" ) :
          self._x_axis_max = val
          narcisse.narsetar ( "x_axis_max" , val )
       elif ( ax == "y" ) :
          self._y_axis_max = val
          narcisse.narsetar ( "y_axis_max" , val )
       elif ( ax == "yr" ) :
          self._yr_axis_max = val
          narcisse.narsetar ( "yr_axis_max" , val )
       elif ( ax == "z" ) :
          self._z_axis_max = val
          narcisse.narsetar ( "z_axis_max" , val )
       elif ( ax == "c" ) :
          self._c_axis_max = val
          narcisse.narsetar ( "c_axis_max" , val )
       else :
          raise self.NarError , "set_axis_max: axis must be x, y, yr, z, or c."
#      narcisse.narsetvals ( self._file_descr )
 
   # set the minimum value of an axis
   def set_axis_min ( self , ax , * val1 ) :
       '''set_axis_min (ax, val) where ax is "x", "y", "z", "yr", or "c".
          The minimum of the specified axis will be set to val.
          val should be a PyFloat object.'''
 
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if len ( val1 ) == 0 :
          val = 0.0
       else :
          val = val1 [0]
       if ( ax == "x" ) :
          self._x_axis_min = val
          narcisse.narsetar ( "x_axis_min" , val )
       elif ( ax == "y" ) :
          self._y_axis_min = val
          narcisse.narsetar ( "y_axis_min" , val )
       elif ( ax == "yr" ) :
          self._yr_axis_min = val
          narcisse.narsetar ( "yr_axis_min" , val )
       elif ( ax == "z" ) :
          self._z_axis_min = val
          narcisse.narsetar ( "z_axis_min" , val )
       elif ( ax == "c" ) :
          self._c_axis_min = val
          narcisse.narsetar ( "c_axis_min" , val )
       else :
          raise self.NarError , "set_axis_min: axis must be x, y, yr, z, or c."
#      narcisse.narsetvals ( self._file_descr )

   # Send axes limits at the last moment before a plot
   def _send_axes_limits ( self ) :
       narcisse.narsetar ( "x_axis_max" , self._x_axis_max )
       narcisse.narsetar ( "x_axis_min" , self._x_axis_min )
       narcisse.narsetar ( "y_axis_max" , self._y_axis_max )
       narcisse.narsetar ( "y_axis_min" , self._y_axis_min )
       narcisse.narsetar ( "yr_axis_max" , self._yr_axis_max )
       narcisse.narsetar ( "yr_axis_min" , self._yr_axis_min )
       narcisse.narsetar ( "z_axis_max" , self._z_axis_max )
       narcisse.narsetar ( "z_axis_min" , self._z_axis_min )
       narcisse.narsetar ( "c_axis_max" , self._c_axis_max )
       narcisse.narsetar ( "c_axis_min" , self._c_axis_min )
       # narsetvals will be done in send_graph
 
   # Allow Narcisse to calculate the axis limits
   def set_default_axes_limits ( self , * h ) :
       '''set_default_axes_limits () sets narcisse to compute the maximum
       and minimum of the axes depending on the data.'''
 
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if self._xyequal : # compute xy limits myself
          xdist = self._x_axis_max - self._x_axis_min
          ydist = self._y_axis_max - self._y_axis_min
          if xdist > ydist :
             self._y_axis_max = self._y_axis_max + xdist - ydist
          elif ydist > xdist :
             self._x_axis_max = self._x_axis_max + ydist - xdist
          narcisse.narsetar ( "x_axis_max" , self._x_axis_max )
          narcisse.narsetar ( "y_axis_max" , self._y_axis_max )
          narcisse.narsetar ( "x_axis_min" , self._x_axis_min )
          narcisse.narsetar ( "y_axis_min" , self._y_axis_min )
#         narcisse.narsetvals ( self._file_descr )
       else :
          narcisse.narsetar ( "x_axis_max" , 0.0 )
          narcisse.narsetar ( "y_axis_max" , 0.0 )
          narcisse.narsetar ( "yr_axis_max" , 0.0 )
          narcisse.narsetar ( "z_axis_max" , 0.0 )
          narcisse.narsetar ( "c_axis_max" , 0.0 )
          narcisse.narsetar ( "x_axis_min" , 0.0 )
          narcisse.narsetar ( "y_axis_min" , 0.0 )
          narcisse.narsetar ( "yr_axis_min" , 0.0 )
          narcisse.narsetar ( "z_axis_min" , 0.0 )
          narcisse.narsetar ( "c_axis_min" , 0.0 )
#         narcisse.narsetvals ( self._file_descr )
 
   # routines to set the limits on individual axes
   def set_x_axis_limits ( self , val1 , * val2i ) :
       '''set_x_axis_limits (min, max) sets the limits on the x axis to
       the specified (pyFloat) sizes.'''
 
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if len ( val2i ) == 0 :
          val2 = 0.0
       else :
          val2 = val2i [0]
       self._x_axis_min = val1
       self._x_axis_max = val2
       narcisse.narsetar ( "x_axis_max" , val2)
       narcisse.narsetar ( "x_axis_min" , val1 )
#      narcisse.narsetvals ( self._file_descr )
 
   def set_y_axis_limits ( self , val1 , * val2i ) :
       '''set_y_axis_limits (min, max) sets the limits on the y axis to
       the specified (pyFloat) sizes.'''
 
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if len ( val2i ) == 0 :
          val2 = 0.0
       else :
          val2 = val2i [0]
       self._y_axis_min = val1
       self._y_axis_max = val2
       narcisse.narsetar ( "y_axis_max" , val2)
       narcisse.narsetar ( "y_axis_min" , val1 )
#      narcisse.narsetvals ( self._file_descr )
 
   def set_yr_axis_limits ( self , val1 , * val2i ) :
       '''set_yr_axis_limits (min, max) sets the limits on the yr axis to
       the specified (pyFloat) sizes.'''
 
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if len ( val2i ) == 0 :
          val2 = 0.0
       else :
          val2 = val2i [0]
       self._yr_axis_min = val1
       self._yr_axis_max = val2
       narcisse.narsetar ( "yr_axis_max" , val2)
       narcisse.narsetar ( "yr_axis_min" , val1 )
#      narcisse.narsetvals ( self._file_descr )
 
   def set_z_axis_limits ( self , val1 , * val2i ) :
       '''set_z_axis_limits (min, max) sets the limits on the z axis to
       the specified (pyFloat) sizes.'''
 
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if len ( val2i ) == 0 :
          val2 = 0.0
       else :
          val2 = val2i [0]
       self._z_axis_min = val1
       self._z_axis_max = val2
       narcisse.narsetar ( "z_axis_max" , val2)
       narcisse.narsetar ( "z_axis_min" , val1 )
#      narcisse.narsetvals ( self._file_descr )
 
   def set_c_axis_limits ( self , val1 , * val2i ) :
       '''set_c_axis_limits (min, max) sets the limits on the c axis to
       the specified (pyFloat) sizes.'''
 
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if len ( val2i ) == 0 :
          val2 = 0.0
       else :
          val2 = val2i [0]
       self._c_axis_min = val1
       self._c_axis_max = val2
       narcisse.narsetar ( "c_axis_max" , val2)
       narcisse.narsetar ( "c_axis_min" , val1 )
#      narcisse.narsetvals ( self._file_descr )
 
   # stuff to help set 3d options
   # (1) These are the legal arguments and their values if wire shows
   legal_3d_options = { 'wm' : 0 , 'w3' : 1 , 'w4' : 3 , 'f3' : 8 , \
                        'f4' : 16 , 'i3' : 32 , 'i4' : 64 , 's3' : 128 , \
                        's4' : 256 , 'none' : 0}
   # (2) These are the values of the other arguments if there is no wire
   legal_3d_no_wire = { 'f3' : 7 , 'f4' : 15 ,  'i3' : 31 , 'i4' : 63 , \
                        's3' : 127 , 's4' : 255 }
   # (3) The following arguments can occur together; the values given
   #     are used if there is no wire showing. (If wire is present,
   #     the values in legal_3d_options are simply or'ed.
   legal_3d_double = { 'f3' : { 'i3' : 39 , 'i4' : 71 } ,
                       'f4' : { 'i3' : 47 , 'i4' : 79 } ,
                       'i3' : { 'f3' : 39 , 'f4' : 47 } ,
                       'i4' : { 'f3' : 71 , 'f4' : 79 } }
 
   def set_3d_options ( self , color_bar , color_bar_pos , * vals ) :
       """set_3d_options (args) may be called with no argument,
       a single string argument, or a sequence of up to three strings.
       If called with no arguments, the graph display is erased.
       A surface is colored by height in z if a 3d option is
       specified, and by the value of a given function if a 4d
       option is specified. With a wire grid option, the grid
       is colored; with a flat option, the quadrilaterals set
       off by grid lines are colored; with a smooth option,
       the surface itself is colored by height; and with an iso
       option, the contour lines are colored. flat and iso options
       may be used together in any combination. wire grid options
       are independent of the other options. Legal arguments for
       set_3d_options are:
       'wm'--monochrome wire grid; 'w3' and 'w4'--3d and 4d
             coloring of wire grid.
       'f3' and 'f4'--flat 3d and 4d coloring options.
       'i3' and 'i4'--3d and 4d isoline (contour line) options.
       's3' and 's4'--3d and 4d smooth coloring options."""
 
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if self._mono == 1 :
          return
       if len (vals) == 0 :
          vals = ["wm"]
       elif is_scalar (vals) :
          vals = [vals [0]]
       else :
          vals = vals [0]
          if is_scalar (vals) :
             vals = [vals]
       if len (vals) > 3 :
          raise self.NarError , "set_3d_options: too many arguments"
       wire_option = -1 # If this ever gets sent, the graph vanishes
       option = 0
       c_color_bar = 0
       z_color_bar = 0
       for i in range ( len (vals) ) :
          if vals [i] == "s4" or vals [i] == "i4" :
             c_color_bar = color_bar
          if vals [i] == "s3" or vals [i] == "i3" :
             z_color_bar = color_bar
          if not self.legal_3d_options.has_key ( vals [i] ) :
             raise self.NarError , "set_3d_options: "\
                   + vals [i] + " is an illegal option."
          if self.legal_3d_options [vals [i]] <= 3 :
             wire_option = self.legal_3d_options [vals [i]]
       if wire_option != -1 :
          for i in range ( len (vals) ) :
             option = option | self.legal_3d_options [vals [i]]
       elif len (vals) >= 1 :
          if len (vals) == 1 :
             option = self.legal_3d_no_wire [vals [0]]
          elif not self.legal_3d_double.has_key [vals [0]] or \
               not self.legal_3d_double [vals [0]].has_key (vals [1]) :
             print "set_3d_options: illegal combination of options: " \
                    + vals [0] + " and " + vals [1] + "."
             return
          else :
             option = self.legal_3d_double [vals [0]][vals [1]]
       else : # cause graph to commit suicide if no args given
          option = -1
       # at this point the arguments were legal and 'option' has been set.
       narcisse.narsetai ( "option_3d" , option )
       # check out whether a color bar is wanted :
       if c_color_bar :
          if color_bar_pos is not None :
             c_color_bar = 1
             narcisse.narsetai ( "height_c_x_min", color_bar_pos [0, 0])
             narcisse.narsetai ( "height_c_y_min", color_bar_pos [0, 1])
             narcisse.narsetai ( "height_c_x_max", color_bar_pos [1, 0])
             narcisse.narsetai ( "height_c_x_max", color_bar_pos [1, 1])
          else :
             c_color_bar = 2
          narcisse.narsetai ( "height_c_type", c_color_bar )
       elif z_color_bar :
          if color_bar_pos is not None :
             z_color_bar = 1
             narcisse.narsetai ( "height_z_x_min", color_bar_pos [0, 0])
             narcisse.narsetai ( "height_z_y_min", color_bar_pos [0, 1])
             narcisse.narsetai ( "height_z_x_max", color_bar_pos [1, 0])
             narcisse.narsetai ( "height_z_x_max", color_bar_pos [1, 1])
          else :
             z_color_bar = 2
          narcisse.narsetai ( "height_z_type", z_color_bar )
       else :
          narcisse.narsetai ( "height_c_type", 0)
          narcisse.narsetai ( "height_z_type", 0)
       narcisse.narsetvals ( self._file_descr )
 
   # Some other routines to set stuff relating to 3d options
   def set_z_contours ( self , val ) :
       """set_z_contours (arg) sets various properties when doing 3d contour
       (iso), smooth, or flat plots. It accepts one argument, as
       follows:
          if an integer n, sets the number of contours to n. This also
                       clears the contour levels array. Countour levels
                       will be computed automatically from the data.
          if a string: 'lin' plots the contours linearly spaced.
                       'log' plots the contours logarithmically spaced.
          if an array NarFloat: sets the contour levels to the values in the
                       array."""
 
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if type ( val ) == IntType :
          for i in range ( val ) :
             narcisse.narsetari ("height_z", 0.0 , i )
          narcisse.narsetvals ( self._file_descr )
          return
       elif type ( val ) == StringType :
          if val == "log" :
             narcisse.narsetai ("height_z_log", 1)
             narcisse.narsetvals ( self._file_descr )
             return
          elif val == "lin" :
             narcisse.narsetai ("height_z_log", 0)
             narcisse.narsetvals ( self._file_descr )
             return
       elif type ( val ) == ArrayType :
          val = val.astype (NarFloat)
          if len (val.shape) == 1 :
             # Note: when setting a Narcisse array you must do a narsetvals
             # after setting each element. If instead you send a whole list
             # of values all at once, then only the last takes effect and
             # all lower values in the table are cleared.
             for i in range (val.shape [0]) :
                narcisse.narsetari ("height_z", val [i] , i)
                narcisse.narsetvals ( self._file_descr )
             narcisse.narsetar ("height_z_h_min", val [0])
             narcisse.narsetar ("height_z_h_max", val [val.shape [0]-1])
             narcisse.narsetvals ( self._file_descr )
             return

       raise self.NarError , "Wrong type of argument to set_z_contours."
 
   def set_c_contours ( self , val ) :
       """set_c_contours (arg) sets various properties when doing 4d contour
       (iso), smooth, or flat plots. It accepts one argument, as
       follows:
          if an integer n, sets the number of contours to n. This also
                       clears the contour levels array. Countour levels
                       will be computed automatically from the data.
          if a string: 'lin' plots the contours linearly spaced.
                       'log' plots the contours logarithmically spaced.
          if an Array NarFloat: sets the contour levels to the values in the
                       array."""
 
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if type ( val ) == IntType :
          for i in range (val) :
             narcisse.narsetari ("height_c", 0.0 , i )
          narcisse.narsetvals ( self._file_descr )
          return
       elif type ( val ) == StringType :
          if val == "log" :
             narcisse.narsetai ("height_c_log", 1)
             narcisse.narsetvals ( self._file_descr )
             return
          elif val == "lin" :
             narcisse.narsetai ("height_c_log", 0)
             narcisse.narsetvals ( self._file_descr )
             return
       elif type ( val ) == ArrayType :
          val = val.astype (NarFloat)
          if len (val.shape) == 1 :
             # Note: when setting a Narcisse array you must do a narsetvals
             # after setting each element. If instead you send a whole list
             # of values all at once, then only the last takes effect and
             # all lower values in the table are cleared.
             for i in range (val.shape [0]) :
                narcisse.narsetari ("height_c", val [i] , i)
                narcisse.narsetvals ( self._file_descr )
             narcisse.narsetar ("height_c_h_min", val [0])
             narcisse.narsetar ("height_c_h_max", val [val.shape [0]-1])
             narcisse.narsetvals ( self._file_descr )
             return
 
       raise self.NarError , "Wrong type of argument to set_c_contours."
 
   # set the mask (hidden line remover) for 3d
   def set_mask ( self , * val ) :
       """set_mask (arg) determines whether hidden parts of the surface
       will be shown on the graph, and if not, what algorithm
       will be used to determine what is hidden. The allowed
       arguments and masking algorithm are as follows:
       'none'--no masking. in wire grid mode, all grid lines
               are visible.
       'min'--the surface is traced beginning in the corner
              closest to the observer.
       'max'--the surface is traced beginning in the corner
              farthest from the observer.
       'sort'--a cell sorting is carried out to determine the
               masking."""
 
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if len ( val ) == 0 or val [0] == "none" :
          narcisse.narsetai ( "option_3d_mask_type" , 0 ) # default: no mask
       elif val [0] == "min" :
          narcisse.narsetai ( "option_3d_mask_type" , 1 ) # minimum mask
       elif val [0] == "max" :
          narcisse.narsetai ( "option_3d_mask_type" , 2 ) # maximum mask
       elif val [0] == "sort" :
          narcisse.narsetai ( "option_3d_mask_type" , 3 ) # sorted mask
       else :
          raise self.NarError , val [0] + " is not a valid mask type."
       narcisse.narsetvals ( self._file_descr )
 
   # Set language
   def set_language ( self , * val ) :
       """set_language (arg) determines what language the Narcisse GUI will
       be displayed in. Called with no argument, it sets the language
       to English. Otherwise it may be called with 'English', 'French',
       'anglais', or 'francaise'. In a concession to the lazy among us,
       'english' and 'french' are also allowed."""
 
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if len ( val ) == 0 or val [0] == "English" or val [0] == "english" :
          narcisse.narsetac ( "language" , "anglais" )
       elif val [0] == "French" or val [0] == "french" :
          narcisse.narsetac ( "language" , "francais" )
       else : # let the user commit suicide
          narcisse.narsetac ( "language" , val [0] )
       narcisse.narsetvals ( self._file_descr )
 
   # commands to set the angle of view:
   def set_phi ( self , * val ) :
       """set_phi (arg) sets the angle of view, measured from the positive z
       axis. If called with no argument, phi is set to 45 degrees.
       Otherwise it should be called with an integer argument (the angle
       in degrees)."""
 
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if len ( val ) == 0 or val [0] is None :
          narcisse.narsetar ( "height" , 30.0 )
       else :
          narcisse.narsetar ( "height" , 90.0 - val [0] )
       narcisse.narsetvals ( self._file_descr )
 
   def set_theta ( self , * val ) :
       """set_theta (arg) sets the angle of view, measured from the positive x
       axis. If called with no argument, theta is set to 45 degrees.
       Otherwise it should be called with an integer argument (the angle
       in degrees)."""
 
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if len ( val ) == 0 or val [0] is None :
          narcisse.narsetar ( "theta" , -45.0 )
       else :
          narcisse.narsetar ( "theta" , val [0] )
       narcisse.narsetvals ( self._file_descr )
 
   def set_roll ( self , * val ) :
       """set_roll (arg) is the angle of rotation around the line determined
       by set_phi and set_theta. If called with no argument, roll
       is set to zero degrees. Otherwise it should be called with
       an integer argument (the angle in degrees)."""
 
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if len ( val ) == 0 or val [0] is None :
          narcisse.narsetar ( "roll" , 0.0 )
       else :
          narcisse.narsetar ( "roll" , val [0] )
       narcisse.narsetvals ( self._file_descr )

   def set_gnomon (self, val) :
       """set_gnomon (val) does nothing in Narcisse."""
       return
 
   # set the distance of view
   def set_distance ( self , * val ) :
       """set_distance (arg) sets the distance of the view point from a 3d
       graph. If called with no argument, or 0.0, this distance is
       effectively infinite. Otherwise it should be called with a
       real number."""
 
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if len ( val ) == 0 :
          narcisse.narsetar ( "distance" , 0. )
       else :
          narcisse.narsetar ( "distance" , val [0] )
       narcisse.narsetvals ( self._file_descr )
 
   # set whether a curve is drawn as a line, step, or one of a set
   # of symbols. val1 specifies the curve(s) and val2 the type(s).
   # If they're both scalars, set that one curve. If they are both
   # vectors, the shorter length will be used. If val1 is a vector
   # and val2 a scalar, then set all curves to the same type.
   # Note: narsetvals has to be called after each call to one of
   # the indexed routines, or else only the last one set is effective.
   # Bug or feature? I don't know.
   ###################NOTE:
   # Currently val2 is an integer value. Eventually I want to replace
   # it with a character designation.
   ###################
   def set_curve_type ( self , val1 , val2 ) :
       """set_curve_type (arg1, arg2) is used to determine how one or a family of
       curves is to be plotted. It must be called with two arguments.
       The first argument is an integer scalar or array Int giving
       the curve number(s) and the second is an integer scalar
       or array Int describing how the curve(s) should be graphed.
       Curves are numbered starting with 0. The allowed values for
       the second argument are: -1 (do not graph), 0 (normal graph),
       1 (graph as a step function), or else a number of options
       to draw the graph as a set of points denoted by symbols:
       2 (+), 3 (*), 4 (o) , 5 (x) , 6 (.)."""
 
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if type (val1) == IntType and type (val2) == IntType :
           narcisse.narsetaii ( "curve_type" , val2 , val1 )
           narcisse.narsetvals ( self._file_descr )
       elif type (val1) == ArrayType and val1.typecode () == Int \
                                   and type (val2) == IntType :
           for i in range (len (val2)) :
               narcisse.narsetaii ( "curve_type" , val2 , val1 [i] )
               narcisse.narsetvals ( self._file_descr )
       elif not is_scalar (val1) and not is_scalar (val2) :
           # both must be > 1 in length
           r = range (len (val1))
           if len (val2) < len (val1) :
               r = range (len (val2))
           for i in r :
               narcisse.narsetaii ( "curve_type" , val2 [i] , val1 [i] )
               narcisse.narsetvals ( self._file_descr )
       else :
           raise self.NarError, "bad arguments to curve_type."
 
   # set the curve color(s) for one or a set of curves.
   # val1 specifies the curve(s) and val2 the color(s).
   # If they're both scalars, set that one curve. If they are both
   # vectors, the shorter length will be used. If val1 is a vector
   # and val2 a scalar, then set all curves to the same color.
   ###################NOTE:
   # Currently val2 is an integer value. Eventually I want to replace
   # it with a character designation.
   ###################
   def set_curve_color ( self , val1 , val2 ) :
       """set_curve_color (arg1, arg2) is used to determine how one or a family of
       curves is to be colored. It must be called with two arguments.
       The first argument is an integer scalar or array Int giving
       the curve number(s) and the second is an integer scalar
       or array Int describing how the curve(s) should be colored.
       Curves are numbered starting with 0. The allowed values for
       the second argument are 0 to 63, denoting the index into
       the current palette."""
 
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if is_scalar (val1) and is_scalar (val2) :
           val2 = self._figure_color (val2)
           narcisse.narsetaii ( "curve_color" , val2 , val1 )
           narcisse.narsetvals ( self._file_descr )
       elif not is_scalar (val1) and len (val1) > 1 and is_scalar (val2) :
           val2 = self._figure_color (val2)
           for i in range (len (val2)) :
               narcisse.narsetaii ( "curve_color" , val2 , val1 [i] )
               narcisse.narsetvals ( self._file_descr )
       elif not is_scalar (val1) and not is_scalar (val2) :
           # both must be > 1 in length
           r = range (len (val1))
           if len (val2) < len (val1) :
               r = range (len (val2))
           for i in r :
               val2 [i] = self._figure_color (val2 [i])
               narcisse.narsetaii ( "curve_color" , val2 [i] , val1 [i] )
               narcisse.narsetvals ( self._file_descr )
       else :
           raise self.NarError , "bad parameters to set_curve_color."

   # set the label type for the curves. "end" and "box".
   def set_label_type ( self , val ) :
       """set_label_type (arg) determines whether curve labels will be attached
       to the ends of curves, or enclosed in a box. The allowed
       arguments are thus 'end' and 'box'."""
 
       if ( val == "end" ) :
          narcisse.narsetai ( "curve_label_type" , 0 )
       elif val == "box" :
          narcisse.narsetai ( "curve_label_type" , 1 )
       else :
          raise self.NarError ,\
                   "set_label_type: 'end' and 'box' are the allowed options."
 
   # set the curve label(s) for one or a set of curves.
   # val1 specifies the curve(s) and val2 the label(s).
   # If they're both scalars, set that one curve. If they are both
   # vectors, choose the shorter of the two lengths.
   def set_curve_label ( self , val1 , val2 ) :
       """set_curve_label (arg1, arg2) is used to label one or a set of curves.
       It requires two arguments. The first is an integer scalar or
       array specifying the curve numbers (starting with 1). The
       second is a scalar string or list of strings specifying
       the label(s) of the curve(s)."""
 
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if type (val1) == IntType and type (val2) == StringType :
           narcisse.narsetaci ( "curve_label" , val2 , val1 )
       elif type (val1) == ArrayType and val1.typecode () == Int and \
            type (val2) == ListType and type (val2 [0]) == StringType :
           r = range (len (val1))
           if len (val2) < len (val1) :
               r = range (len (val2))
           for i in r :
               narcisse.narsetaci ( "curve_label" , val2 [i] , val1 [i] )
               narcisse.narsetvals ( self._file_descr )
       else :
           print "Val1: " , `val1`
           print "Val2: " , `val2`
           raise self.NarError ,\
                "set_curve_label: arguments have inconsistent types or sizes."
 
   def set_xyequal (self) :
      """set_xyequal () sets a parameter that makes the axes equal scale."""
      self._xyequal = 1

   def reset_xyequal (self) :
      """set_xyequal () resets a parameter that makes the axes equal scale."""
      self._xyequal = 0

   narcisse_marks = { "+" : 2 , "*" : 3 , "o" : 4 , "x" : 5 , "." : 6 }
   narcisse_types = { "none" : -1 , "hide" : -1 , "line" : 0 ,
                      "normal" : 0 , "step" : 1 }

   def _figure_type ( self , crv ) :
      """_figure_type (crv) makes sure to return a valid type for a
      Narcisse curve.
      """
      if crv.hide :
         return -1
      if crv.marks and crv.marker is None :
         if (type (crv.line_type) == IntType and \
            crv.line_type == 0 or \
            type (crv.line_type) == StringType and \
            (crv.line_type == "line" or crv.line_type == "normal" or \
            crv.line_type == "solid")) :
            if crv.label == " " :
               crv.label = uppercase [self._next_letter]
               self._next_letter = (self._next_letter + 1) % 26
            return 0
      if crv.marks and crv.marker is not None :
         # if a marker is specified but a curve is desired, set the
         # curve's label to the marker
         if type (crv.marker) == StringType and \
            (type (crv.line_type) == IntType and \
            crv.line_type == 0 or \
            type (crv.line_type) == StringType and \
            (crv.line_type == "line" or crv.line_type == "normal" or \
            crv.line_type == "solid")) :
            if crv.label == " " :
               crv.label = crv.marker
            return 0
         if type (crv.marker) == IntType and 2 <= crv.marker <= 6 :
            return crv.marker
         if type (crv.marker) == StringType and \
            self.narcisse_marks.has_key (crv.marker) :
            return self.narcisse_marks [crv.marker]
      if type (crv.line_type) == IntType and -1 <= crv.line_type <= 6 :
         return crv.line_type
      if type (crv.line_type) == StringType and \
         self.narcisse_types.has_key (crv.line_type) :
         return self.narcisse_types [crv.line_type]
      return 0 # incomprehensible so draw a line
         
   def plot_object ( self , crv ) :

       """plot_object (crv) is a general purpose plotting routine. It should
       be called with one argument, a curve (all that Narcisse currently
       accepts). In the case of multiple objects on one graph, the
       first call only should be to this routine, subsequent calls to
       add_object. For Narcisse, plot_object and add_object accumulate
       information about the various curves, then send all the freight
       when send_graph is called.
       """
       try :
          dum = crv.type ()
       except :
          raise self.NarError , \
             "Unknown object has been sent to Narcisse."
       if dum != CurveType :
          raise self.NarError , \
             "Narcisse does not know how to graph a " + dum + "."
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       self._graph_type = 2
       # We compute new axis limits if user wants equal scales
       if self._xyequal :
          new_x_limits = minmax (crv.x)
          new_y_limits = minmax (crv.y)
          self._x_axis_min = new_x_limits [0]
          self._x_axis_max = new_x_limits [1]
          self._y_axis_min = new_y_limits [0]
          self._y_axis_max = new_y_limits [1]
       # start a list of attributes for curves
       self._types = [self._figure_type (crv)]
       self._labels = [crv.label]
       self._colors = [crv.color]
       self._axispref = [crv.axis]
       self._ylist = [crv.y]
       self._xlist = [crv.x]

   # add a curve to an existing plot
   def add_object ( self , crv ) :
       """add_object (crv) will add a curve to an existing graph.
       The curve's attributes are saved up; nothing is sent to Narcisse
       until send_graph is called.
       """
       try :
          dum = crv.type ()
       except :
          raise self.NarError , \
             "Unknown object has been sent to Narcisse."
       if dum != CurveType :
          raise self.NarError , \
             "Narcisse does not know how to graph a " + dum + "."
       if self._graph_type != 2 :
          raise self.NarError , \
             "plot_object must be called for the first curve on a graph."
       if self._xyequal :
          new_x_limits = minmax (crv.x)
          new_y_limits = minmax (crv.y)
          if self._x_axis_min > new_x_limits [0] :
             self._x_axis_min = new_x_limits [0]
          if self._x_axis_max < new_x_limits [1] :
             self._x_axis_max = new_x_limits [1]
          if self._y_axis_min > new_y_limits [0] :
             self._y_axis_min = new_y_limits [0]
          if self._y_axis_max < new_y_limits [1] :
             self._y_axis_max = new_y_limits [1]
       self._types.append (self._figure_type (crv))
       self._labels.append (crv.label)
       self._colors.append (crv.color)
       self._ylist.append (crv.y)
       self._xlist.append (crv.x)
       self._axispref.append (crv.axis)

   def _send_2d_info ( self ) :
       """_send_2d_info ( ) sends the accumulated curve information
       out to Narcisse.
       """
       n = len (self._ylist)
       if n <= 0 :
          raise self.NarError, \
             "There is nothing to graph!"
       elif n == 1 :
          arg1 = 0
          arg2t = self._types [0]
          arg2c = self._colors [0]
          arg2l = self._labels [0]
       else :
          arg1 = arange (n, typecode = Int)
          arg2t = self._types
          arg2c = self._colors
          arg2l = self._labels
       y = self._ylist [0].astype (NarFloat)
       x = self._xlist [0].astype (NarFloat)
       narcisse.nar1curve (self._file_descr, y, x)
       self.set_y_axis (0, self._axispref [0])
       for i in range (1, n) :
          y = self._ylist [i].astype (NarFloat)
          x = self._xlist [i].astype (NarFloat)
          narcisse.narsetai ( "option_2d_concatenate" , 1 )
          narcisse.narsetvals ( self._file_descr )
          narcisse.nar1curve (self._file_descr, y, x)
          narcisse.narsetai ( "option_2d_concatenate" , 0 )
          narcisse.narsetvals ( self._file_descr )
          self.set_y_axis (i, self._axispref [i])
       self.set_curve_type (arg1, arg2t)
       self.set_curve_color (arg1, arg2c)
       self.set_curve_label (arg1, arg2l)

   def plot_text ( self ) :
       "plot_text does nothing in Narcisse."
       return

   def set_text ( self , txt , n ) :
       "set_text (str, ix) sets the ix'th text to str."
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if txt == " " : txt = ""
       narcisse.narsetaci ( "text_value" , txt , n )
       narcisse.narsetvals ( self._file_descr )

   def clear_text (self) :
       "clear_text ( ) sets the number of texts to 0."
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       narcisse.narsetaci ( "text_value" , "" , 0 )
       narcisse.narsetai ( "text_number" , 0 )
       narcisse.narsetvals ( self._file_descr )

   def set_text_color ( self , txt , n ) :
       """set_text_color (col, ix) sets the ix'th text color to col,
       which is a number between 0 and 63 associated with a color table."""

       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       txt = self._figure_color (txt)
       narcisse.narsetaii ( "text_color" , txt , n )
       narcisse.narsetvals ( self._file_descr )

   def set_text_size ( self , txt , n ) :
       """set_text_size (sz, ix) sets the ix'th text size to sz.
       sz represents essentially the number of characters that
       will fill the width of the graphics screen, so the larger
       the number, the smaller the text."""

       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       narcisse.narsetaii ( "text_size" , txt , n )
       narcisse.narsetvals ( self._file_descr )

   def set_text_pos ( self , x , y , ix ) :
       """set_text_pos (x, y, ix) positions the ix'th text at (x, y),
       which are real numbers between 0 and 1 giving relative
       position in the graphics window."""

       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       narcisse.narsetari ("text_pos_x", x, ix)
       narcisse.narsetari ("text_pos_y", y, ix)
       narcisse.narsetvals ( self._file_descr )


   # Here's the grandaddy of them all, a perfectly general surface
   # plotting routine. Note that it passes a lot of information
   # to narcissemodule for error checking.
   ###############################################################
   # Eventually these routines should probably all be rewritten
   # to accept numerical sequences of any kind as inputs,
   # convert them to array types as appropriate, check for
   # appropriate dimensions, etc. The problem is that anybody
   # can call the low level routines directly from Python, so they
   # need to do error checking anyhow, just in case.
   ###############################################################
   def plot_surface ( self , arg1 , * args2 ) :
       """plot_surface (args) is a general-purpose 3d/4d plotting routine.
       The type of plot depends on the numbers and types of the
       arguments (which all must be of type array NarFloat except for the
       cell information for unstructured grids). Here we go:
          1. single argument, two dimensional array z: Plot z as a
             surface versus equally spaced x and y coordinates.
          2. three arguments, two vectors x and y and a two dimensional
             matrix z: plot z as a surface versus the given x and y.
          3. three arguments, matrices x, y, and z (whose dimensions
             must match): plot z as a surface versus the given x and y.
          4. four arguments, two vectors x and y and two two dimensional
             matrices z and c: plot z as a surface versus the given x
             and y; use the variable c to color the graph.
          5. four arguments, matrices x, y, z, and c (whose dimensions
             must match): plot z as a surface versus the given x and y;
             use the variable c to color the graph.
          6. four arguments, three vectors x, y, and z specifying a
             structured grid, and a three-dimensional array c defined
             at each grid point: draw the grid and color according to
             the variable c.
          7. six arguments, vectors x, y, and z of the same size
             specifying a nonstructured grid, and c of the same size
             specifying a value at each point; cd, an integer vector
             specifying connectivity (see the Narcisse manual for
             details), and nc an integer specifying the number of
             cells in the grid, draw the nonstructured grid and color
             according to the variable c."""

       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if len ( args2 ) == 0 :
          narcisse.narsurf ( self._file_descr , arg1.astype (NarFloat) )
          return
       if len ( args2 ) == 2 or len ( args2 ) > 2 and args2 [2] is None :
          x = arg1.astype (NarFloat)
          y = args2 [0].astype (NarFloat)
          z = args2 [1].astype (NarFloat)
          if len (x.shape) == 1 and len (y.shape) == 1 :
             narcisse.nar3drect ( self._file_descr , x , y , z )
             return
          else :
             narcisse.nar3dtetra ( self._file_descr , x , y , z )
             return
       if len ( args2 ) == 3 :
          x = arg1.astype (NarFloat)
          y = args2 [0].astype (NarFloat)
          z = args2 [1].astype (NarFloat)
          c = args2 [2].astype (NarFloat)
          if ( len ( z.shape ) == 2 ) :
             if len ( x.shape ) == 1 :
                narcisse.nar4drect ( self._file_descr , x , y , z , c )
                return
             else :
                narcisse.nar4dtetra ( self._file_descr , x , y , z , c )
                return
          else :
             narcisse.narstructmesh ( self._file_descr , x , y , z , c )
             return
       if len ( args2 ) == 5 :
          x = arg1.astype (NarFloat)
          y = args2 [0].astype (NarFloat)
          z = args2 [1].astype (NarFloat)
          c = args2 [2].astype (NarFloat)
          cd = args2 [3].astype (Int)
          nc = args2 [4]
          narcisse.narnonstructmesh ( self._file_descr , x , y , z , c , cd , nc)
          return
       else :
          n=1+len (args2)
          raise self.NarError ,\
                   "plot_surface: inappropriate number of arguments: " + `n`

   def set_palette (self, col) :
       """set_palette (col) sets the color palette to col. The first
       entry in col tells how long the rest of the array is; then
       there are col [0] / 3 entries for red, followed by the same
       number of greens, followed by the same number of blues.
       """
       for i in range (col [0]) :
           narcisse.narsetaii ("parameter_map_pal", col [i + 1], i)
       narcisse.narsetai ("parameter_map", -1)
       narcisse.narsetvals ( self._file_descr )

   def set_no_concat ( self ) : #called by a graphics object initially
       "set_no_concat () turns off the 2d and 3d concatenation mode."
       narcisse.narsetai ("option_3d_concatenate", 0)
       narcisse.narsetai ("option_2d_concatenate", 0)
       narcisse.narsetvals ( self._file_descr )

   # add a surface to an existing plot
   def add_surface ( self , arg1 , * args2 ) :
       """add_surface (args) will add one or more surfaces to an existing graph.
       Its arguments are the same form as the arguments of
       plot_surface. See plot_surface documentation for details."""

       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       narcisse.narsetai ("option_3d_concatenate" , 1 )
       narcisse.narsetvals ( self._file_descr )
       if len (args2) == 0 :
          self.plot_surface ( arg1 )
       elif len (args2) == 2 :
          self.plot_surface ( arg1 , args2 [0] , args2 [1] )
       elif len (args2) == 3 :
          self.plot_surface ( arg1 , args2 [0] , args2 [1] , args2 [2] )
       elif len (args2) == 5 :
          self.plot_surface ( arg1 , args2 [0] , args2 [1] , args2 [2] , args2 [3] , args2 [4] )
       else :
          raise self.NarError , "add_surface: inappropriate number of arguments ("\
                            + `len (args2)` + ")."
       narcisse.narsetai ("option_3d_concatenate" , 0 )
       narcisse.narsetvals ( self._file_descr )

   # routine to freeze the graph
   # i. e., arguments and graphs sent will not be plotted until
   # send_graph is called.
   def freeze_graph ( self ) :
       """freeze_graph () keeps a graph from being plotted until
       send_graph () is called."""

       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if not self._frozen :
          narcisse.narsetai ( "plot_now" , 0 )
          narcisse.narsetvals ( self._file_descr )
          self._frozen = 1


   # routine to release the graph
   # The current graph will be plotted and any arguments will
   # be sent. if _freeze_each has been set, then the next graph
   # will be _frozen too.
   def send_graph ( self, graf ) :
       """send_graph () causes a plot that has been accumulated
       after freeze_graph () was called, to be plotted."""
 
       if not self._file_open : raise self.ConnectException , \
          "You are not connected to Narcisse."
       if self._graph_type == 0 :
          raise self.NarError, \
             "There is nothing to graph!"
       if graf.type () == Graph3dType :
          self._dims = 3
       else :
          self.dims = 2
       self._send_axes_limits ( )
       if self._graph_type == 2 :
          self._send_2d_info ( )
       if self._frozen :
          self._frozen = 0
          narcisse.narsetai ( "plot_now" , 1 )
       narcisse.narsetvals ( self._file_descr )
       if self._freeze_each :
          self.freeze_graph ( )

   def set_freeze_each ( self , val ) :
       """set_freeze_each ( fe ) tells whether or not to re-freeze the
       graphics after each send_graph call. 1 to re-freeze, 0 not to."""
       self._freeze_each = val

   def send_generics (self, graf) :
       """send_generics ( graf ) sets up all the things that are generic to any
       graph. It does not actually do any plotting yet.
       """
       self.set_titles ( graf._titles )
       self.set_title_colors ( graf._title_colors )
       # The following line is redundant for Gist. Not sure about Narcisse.
       # self._plot_titles ( )
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

   # The following is the equivalent of the Gist split palette.
   # The lower half is the rainbow, the upper half is greyscale.
   split_palette = array ([162, # The last 54 colors of the palette
                           #27 reds, 27 greys:
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           19, 67, 115, 163, 211,
                           255, 255, 255, 255, 255, 255, 
                           255, 255, 255, 255, 255,
                           0, 9, 19,  29,  39,  49,  58,  68,
                           78,  88,  98, 107, 117, 127, 137,
                           147, 156, 166, 176, 186, 196, 205,
                           215, 225, 235, 245, 255,
                           #27, greens, 27 greys:
                           24, 72, 120, 168, 216, 
                           255, 255, 255, 255, 255, 255,
                           255, 255, 255, 255, 255,
                           226, 178, 129, 81, 33, 0, 0, 0, 0, 0, 0,
                           0, 9, 19,  29,  39,  49,  58,  68,
                           78,  88,  98, 107, 117, 127, 137,
                           147, 156, 166, 176, 186, 196, 205,
                           215, 225, 235, 245, 255,
                           #27 blues, 27 greys:
                           255, 255, 255, 255, 255,
                           245, 197, 149, 101, 52, 4, 
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           14, 62, 110, 158, 206, 255,
                           0, 9, 19,  29,  39,  49,  58,  68,
                           78,  88,  98, 107, 117, 127, 137,
                           147, 156, 166, 176, 186, 196, 205,
                           215, 225, 235, 245, 255
                          ], Int)

   def do_generic (self, graf) :
       self.set_freeze_each (1)
       self.freeze_graph ( ) #freeze everything until entire graph is sent
       self.set_no_concat ( )
       self.send_generics ( graf )
       self.set_axis_labels ( graf._axis_labels )
       self.set_x_axis_limits (graf._axis_limits [0][0],
                               graf._axis_limits [0][1])
       self.set_y_axis_limits (graf._axis_limits [1][0],
                               graf._axis_limits [1][1])
       if self._dims == 2 :
          self.set_yr_axis_limits (graf._axis_limits [2][0],
                                  graf._axis_limits [2][1])
       elif self._dims == 3 :
          self.set_z_axis_limits (graf._axis_limits [2][0],
                                  graf._axis_limits [2][1])
          self.set_c_axis_limits (graf._axis_limits [3][0],
                                  graf._axis_limits [3][1])
          self.set_yr_axis_limits (graf._axis_limits [4][0],
                                  graf._axis_limits [4][1])
       if self._dims == 2:
          for i in range (graf._no_of_axes) :
             if graf._axis_scales [i] == "lin" :
                self.set_axis_lin (graf._axes [i])
             elif graf._axis_scales [i] == "log" :
                self.set_axis_log (graf._axes [i])
             else :
                raise graf._AxisSpecError , \
                   graf._axis_scales [i] + " is not a valid axis scale."
       elif self._dims == 3:
           sc = graf._axis_scales
           if is_scalar (sc) :
              sc = [sc] + ["lin", "lin", "lin", "lin"]
           else :
              for i in range (5 - no_of_dims (sc)) :
                 sc = sc + ["lin"]
           for i in range (5) :
              if sc [i] == "log" :
                 self.set_axis_log (graf._axes [i])
              else : # anything else will be lin
                 self.set_axis_lin (graf._axes [i])
       try:
           no_color = os.environ["NO_COLOR"]
       except KeyError:
           no_color = 0
       if no_color == 0 or no_color == "no" or no_color == "n" :
           self.set_color_card (graf._color_card , 1)
       if self._dims == 3:
           self.set_phi ( graf._phi )
           self.set_theta ( graf._theta )
           self.set_roll ( graf._roll )
       return

   def quick_plot (self, graf) :
       "quick_plot (graf) plots without recomputing."
       if graf.type () == Graph2dType :
          self._dims = 2
       else :
          self._dims = 3
       self.do_generic (graf)
       if hasattr (graf, "n") and self._dims == 3 :
          if graf.opt_3d_change :
             self.set_3d_options ( graf._color_bar,
                                   graf._color_bar_pos,
                                   graf._s [graf.n - 1].opt_3d )
          if graf.mask_change :
             self.set_mask ( graf._s [graf.n - 1].mask )
          if graf.mesh_type_change :
             self.set_3d_grid_type ( graf._s [graf.n - 1].mesh_type )
       if hasattr (graf, "n") and self._dims == 2 :
          if graf.type_change :
             self.set_curve_type ( graf.n - 1 , graf._c[graf.n - 1].line_type )
          if graf.color_change :
             self.set_curve_color ( graf.n - 1 , graf._c[graf.n - 1].color )
          if graf.label_change :
             self.set_curve_label ( graf.n - 1 , graf._c[graf.n - 1].label )
       if graf._label_type != " " :
          self.set_label_type ( graf._label_type )
       if (graf._sync) :
          self.synchronize ( )
       self.send_graph (graf)

   def plot2d (self, graf) :
       """A Graph2d object calls plot2d with itself as argument.
       plot2d sorts out everything for the graph and then does the plot.
       The bulk of this work used to be done in Graph and Graph2d,
       but I decided it was too graphics-dependent.
       """
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
          else :
             self.add_object ( graf._c [i] )
       self.plot_text ( )
       # Finally do the graph
       if (graf._sync) :
          self.synchronize ( )
       self.send_graph (graf)

   def split_bytscl (self, val, top) :
       """
       split_bytscl (val, top) scales the values in val to the top
       half of the palette (values 27 to 53) if top = 1, and to
       the bottom half (values 0 to 26) if top = 0.
       """
       retval = ( (val - min (val)).astype(Float) /
                  max( (val - min (val)).astype(Float)*26. +
                  0.5)).astype (Int) + top * 26
       

   def plot3d (self, graf) :
       """plot3d (graf) plots a 3d graph object.
       """
       self._dims = 3
       # (1) Do graph-generic stuff first
       self.do_generic (graf)

       self.set_phi ( graf._phi )
       self.set_theta ( graf._theta )
       self.set_roll ( graf._roll )
       self.set_distance ( graf._distance )
       n = graf._s_ln
       if n > 1 and graf._connect :
          self.set_connect ( 1)
       else :
          self.set_connect ( 0)
       if graf._s [0].type () == Slice3dType :
          # This is a graph of one or more isosurface and/or plane slices.
          # Basically, we just need to put the vertices and cell
          # information into the form recognized by SpxNonStruct4d.
          # For now, Narcisse will not allow slices to be combined
          # with other surfaces.
          # send out surface characteristics, then each surface
          self.set_link ( 0 )
          self.set_mask ( graf.mask )
          self.set_3d_options ( graf._color_bar,
                                graf._color_bar_pos,
                                graf._s [0].opt_3d )
          self.set_3d_grid_type ( graf._s[0].mesh_type )
          self.set_z_c_switch ( graf._s[0].z_c_switch )
          self.set_z_contours ( graf._s[0].z_contours_scale )
          self.set_c_contours ( graf._s[0].c_contours_scale )
          if graf._s[0].z_contours_array is None :
             if graf._s[0].number_of_z_contours is None :
                self.set_z_contours (20)
             else :
                self.set_z_contours (graf._s[0].number_of_z_contours)
          else :
             self.set_z_contours ( graf._s[0].z_contours_array )
          if graf._s[0].c_contours_array is None :
             if graf._s[0].number_of_c_contours is None :
                self.set_c_contours (20)
             else :
                self.set_c_contours (graf._s[0].number_of_c_contours)
          else :
             self.set_c_contours ( graf._s[0].c_contours_array )

          isosurfaces_present = 0
          self._graph_type = 4
          for i in range (graf._s_ln) :
             if graf._s [i].type () == Slice3dType :
                if graf._s [i].plane is None and graf._s [i].iso is not None :
                   isosurfaces_present = 1
             else :
                raise self.NarError, \
                   "If one component is a Slice, all must be."
          for i in range (graf._s_ln) :
             s = graf._s [i]
             opt_3d = s.opt_3d
             if type (opt_3d) != ListType :
                opt_3d = [opt_3d]
             if i == 0 :
                nv = s.nv
                x = s.xyzv [:, 0]
                y = s.xyzv [:, 1]
                z = s.xyzv [:, 2]
                if (max (abs (x)) < 10.e-30) :
                   x [0: len(x)] = 0.
                if (max (abs (y)) < 10.e-30) :
                   y [0: len(y)] = 0.
                if (max (abs (z)) < 10.e-30) :
                   z [0: len(z)] = 0.
                if not isosurfaces_present or s.iso is None and \
                   s.plane is None :
                   if "i3" in opt_3d or "s3" in opt_3d or \
                      "w3" in opt_3d or "f3" in opt_3d :
                      val = z
                   else :
                      val = s.val
                elif s.plane is not None :
                   if len(s.val) == len (s.nv) :
                      val = to_corners (s.val, s.nv, sum (s.nv))
                   else :
                      val = s.val
                else :
                   val = ones (sum (s.nv), Float) * s.iso
             else :
                nv = concatenate ( (nv, s.nv))
                x = concatenate ( (x, s.xyzv [:, 0]))
                y = concatenate ( (y, s.xyzv [:, 1]))
                z = concatenate ( (z, s.xyzv [:, 2]))
                if not isosurfaces_present or s.iso is None and \
                   s.plane is None :
                   val = concatenate ( (val, s.val))
                elif s.plane is not None :
                   if len(s.val) == len (s.nv) :
                      val = concatenate ( (val,
                         to_corners (s.val, s.nv, sum (s.nv))))
                   else :
                      val = concatenate ( (val, s.val))
                else :
                   val = concatenate ( (val, ones (sum (s.nv), Float) * s.iso))
          nc = len (nv)
          nv = concatenate ( (cumsum (nv), arange (len (x))))
##        if isosurfaces_present :
##           self.set_palette (self.split_palette)
          self.set_color_card (graf._color_card)
          self.plot_surface (x, y, z, val, nv, nc)
          if (graf._sync) :
             self.synchronize ( )
          self.plot_text ( )
          self.send_graph (graf)
       elif graf._link :
          # got to send out one surface and its characteristics at a time
          self.set_link ( 1 )
          for i in range ( n ) :
             # Do not replot axes for subsequent components
             if i > 0 : self.set_grid_type ("none")
             self.set_mask ( graf._s[i].mask )
             self.set_3d_options ( graf._color_bar,
                                   graf._color_bar_pos,
                                   graf._s[i].opt_3d )
             self.set_3d_grid_type ( graf._s[i].mesh_type )
             if graf._s[i].z_c_switch :
                self.set_z_c_switch ( 1 )
             else :
                self.set_z_c_switch ( 0 )
             self.set_z_contours ( graf._s[i].z_contours_scale )
             self.set_c_contours ( graf._s[i].c_contours_scale )
             if graf._s[i].z_contours_array is None :
                if graf._s[i].number_of_z_contours is None :
                   self.set_z_contours (20)
                else :
                   self.set_z_contours (graf._s[i].number_of_z_contours)
             else :
                self.set_z_contours ( graf._s[i].z_contours_array )
             if graf._s[i].c_contours_array is None :
                if graf._s[i].number_of_c_contours is None :
                   self.set_c_contours (20)
                else :
                   self.set_c_contours (graf._s[i].number_of_c_contours)
             else :
                self.set_c_contours ( graf._s[i].c_contours_array )
             # always send coordinates of linked surfaces
             if not hasattr (graf._s[i], "x") or \
                graf._s[i].x is None : # just graphing z
                self._graph_type = 3
                self.plot_surface ( array ( graf._s[i].z, Float))
             elif graf._s[i].c is None : #surface alone
                self._graph_type = 3
                self.plot_surface ( array ( graf._s[i].x, Float),
                                 array ( graf._s[i].y, Float),
                                 array ( graf._s[i].z, Float))
             else : # 4d plot (surface or structured mesh plot)
                self._graph_type = 4
                if graf._s[i].type () == SurfaceType or \
                   graf._s[i].structured :
                   # (surface or structured mesh plot)
                   self.plot_surface ( array ( graf._s[i].x, Float),
                                    array ( graf._s[i].y, Float),
                                    array ( graf._s[i].z, Float),
                                    array ( graf._s[i].c, Float))
                else : # Nonstructured mesh
                   graf._s[i].create_Narcisse_format ()
                   self.plot_surface ( array ( graf._s[i].x, Float),
                       array ( graf._s[i].y, Float),
                       array ( graf._s[i].z, Float),
                       array ( graf._s[i].c, Float),
                       array ( graf._s[i].cell_descr, Int ),
                       graf._s[i].number_of_cells)
             if (graf._sync) :
                self.synchronize ( )
             self.send_graph (graf)
             if i == 0 :
                self.plot_text ( )
                self.set_link (0)
       else : # not graf._link
          # send out surface characteristics, then each surface
          self.set_link ( 0 )
          self.set_mask ( graf._s[n - 1].mask )
          self.set_3d_options ( graf._color_bar,
                                graf._color_bar_pos,
                                graf._s[n - 1].opt_3d )
          self.set_3d_grid_type ( graf._s[n - 1].mesh_type )
          self.set_z_c_switch ( graf._s[n - 1].z_c_switch )
          self.set_z_contours ( graf._s[n - 1].z_contours_scale )
          self.set_c_contours ( graf._s[n - 1].c_contours_scale )
          if graf._s[n - 1].z_contours_array is None :
             if graf._s[n - 1].number_of_z_contours is None :
                self.set_z_contours (20)
             else :
                self.set_z_contours (graf._s[n - 1].number_of_z_contours)
          else :
             self.set_z_contours ( graf._s[n - 1].z_contours_array )
          if graf._s[n - 1].c_contours_array is None :
             if graf._s[n - 1].number_of_c_contours is None :
                self.set_c_contours (20)
             else :
                self.set_c_contours (graf._s[n - 1].number_of_c_contours)
          else :
             self.set_c_contours ( graf._s[n - 1].c_contours_array )
          if graf._send_coordinates :
             for i in range ( n ) : # now send out surfaces
                if not hasattr (graf._s[i], "x") or \
                   graf._s[i].x is None : # just graphing z
                   self._graph_type = 3
                   if i == 0 :
                      self.plot_surface (array (graf._s[i].z, Float))
                   else :
                      self.add_surface (array (graf._s[i].z, Float))
                elif graf._s[i].c is None : # 3d plot
                   self._graph_type = 3
                   if i == 0 :
                      self.plot_surface ( array ( graf._s[i].x, Float),
                                       array ( graf._s[i].y, Float),
                                       array ( graf._s[i].z, Float))
                   else :
                      self.add_surface ( array ( graf._s[i].x, Float),
                                       array ( graf._s[i].y, Float),
                                       array ( graf._s[i].z, Float))
                else : # 4d plot (surface or structured mesh plot)
                   self._graph_type = 4
                   if graf._s[i].type () == SurfaceType or \
                      graf._s[i].structured :
                      if i == 0 :
                         self.plot_surface ( array ( graf._s[i].x, Float),
                                          array ( graf._s[i].y, Float),
                                          array ( graf._s[i].z, Float),
                                          array ( graf._s[i].c, Float))
                      else :
                         self.add_surface ( array ( graf._s[i].x, Float),
                                          array ( graf._s[i].y, Float),
                                          array ( graf._s[i].z, Float),
                                          array ( graf._s[i].c, Float))
                   else : # Nonstructured mesh plot
                      graf._s[i].create_Narcisse_format ()
                      if i == 0 :
                         self.plot_surface ( array ( graf._s[i].x, Float),
                             array ( graf._s[i].y, Float),
                             array ( graf._s[i].z, Float),
                             array ( graf._s[i].c, Float),
                             array ( graf._s[i].cell_descr, Int),
                             graf._s[i].number_of_cells )
                      else :
                         self.add_surface ( array ( graf._s[i].x, Float),
                             array ( graf._s[i].y, Float),
                             array ( graf._s[i].z, Float),
                             array ( graf._s[i].c, Float),
                             array ( graf._s[i].cell_descr, Int),
                             graf._s[i].number_of_cells )
          if (graf._sync) :
             self.synchronize ( )
          self.plot_text ( )
          self.send_graph (graf)

   def move_light_source (self, graf, angle, nframes) :
      raise self.NarError, \
         "Sorry, Narcisse does not yet support a moving light source."

   def rotate_graph (self, axis, angle, nframes) :
       # In Narcisse, only the angle counts.
       narcisse.narsetai ("animation_number", nframes)
       angle = angle * 180. / pi
       narcisse.narsetai ("animation_azimuth", angle)
       narcisse.narsetai ("animation_elevation", angle)
       if self._frozen :
          self._frozen = 0
          narcisse.narsetai ( "plot_now" , 1 )
       narcisse.narsetvals ( self._file_descr )
       if self._freeze_each :
          self.freeze_graph ( )
