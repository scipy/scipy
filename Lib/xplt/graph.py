# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

GraphicsError = "GraphicsError"
from types import *
from graftypes import *
import os
try:
   graphics = os.environ["PYGRAPH"]
except KeyError:
   print "env variable PYGRAPH is not set, so assuming Gist graphics."
   print "...this can be overruled by the 'graphics' keyword argument"
   print "   to Graph2d or Graph3d."
   graphics = "Gist"

 
if graphics [0:3] == "Nar" :
   import NarPlotter
   Graphics = NarPlotter
elif graphics == "Gist" :
   import GistPlotter
   Graphics = GistPlotter
else :
   raise GraphicsError , \
      graphics + " is an unknown graphics package. Check PYGRAPH " + \
         "environment variable."
 
# The following is so I know about arrays:
from Numeric import *
from scipy_base.fastumath import *
from shapetest import *

class Graph :
 
   """
   g = Graph ( <keyword arguments> ) abstracts all the information
   that is common to 2d, 3d, and 4d graphs, and will be inherited from
   by all of those types. It is never intended to be instantiated
   except via a derived class. The keyword arguments are:
 
       plotter = <Plotter object> or a sequence of <Plotter object>s
              if you have a plotter object or a sequence of them that
              you want the graph to use when plotting itself. I
              recommend against this; a graph can create its own
              plotter if you don't give it one.
       filename = <string> or a sequence of <string>s if you want to
              connect to a Plotter with the named file(s). (The default
              is " ".) Only one of the keywords 'plotter' or 'filename'
              is allowed, and both are optional.
              NOTE: the possibility of a sequence of file names or
              plotters allows one to display Plotter on one's
              own machine as well as one or more remote machines.
              In Gist, the filename argument is the same as the
              display argument (below).
       hardcopy = <string> Allows the user to specify a hardcopy file
              for this Plotter object. If string ends in ".ps", this will
              be postscript; it ".cgm", then it will be cgm.
              This is valid only for Gist; ignored for Narcisse.
       display = <string> or a sequence of <string>s if you want to
              display on the named hosts. The form of <string> is
              the usual "hostname:server.screen". The purpose of
              this argument is to allow you to continue without
              exiting python if you have to open a Gist window
              without the DISPLAY environment variable having been
              set, or if you want to open Gist windows on more
              than one host.
       graphics = <string> or a sequence of <string>s if you want to
              specify which graphics the particular plotter connects to.
              Currently the values allowed are "Nar" and "Gist".
              This argument is meaningless if you supply a list
              of plotters. If it is a scalar and you supply a list
              of filenames, all plotters opened will be that type.
              If it is a vector, then it must match the list of
              filenames in size and correspond to the filename.
       titles = <value> where <value> is a string or a sequence
                of up to four strings, giving the titles in the
                order bottom, top, left, right.
       title_colors = <value> where value is an integer or string
                or a sequence of up to four integers or strings giving
                the colors of the titles.
       axis_labels = <value> where <value> is a string or
                for a 2d Graph: a sequence of up to three strings
                representing the labels of the x axis, the y axis,
                and the right y axis. for a 3d Graph: a sequence of
                up to four strings representing the labels of the
                x axis, the y axis, the z axis, and the right y axis.
       axis_limits = <value> where <value> is a pair [xmin, xmax] or
                 for a 2d Graph: a sequence of up to three pairs,
                 where the second would be the y limits and the third
                 the right y limits. for a 3d or 4d graph:
                 a sequence of up to five pairs, where the second
                 would be the y limits, the third the z limits,
                 the fourth the c limits, and the fifth the right
                 y limits. axis_limits = "defaults" will let the
                 graphics calculate the limits.
          x_axis_limits, y_axis_limits, yr_axis_limits, z_axis_limits,
                 and c_axis_limits are used to change the limits on
                 individual axes. They may be set to a pair or to
                 "defaults".
       axis_scales = "linlin", "linlog", "loglin", or "loglog" for
                2d graphs with no separate right y axis; if there is
                a separate yr axis, and it is not to default to
                linear, then the scales must be expressed as a
                sequence of three values "log" or "lin" with the
                third one (for yr) "log". For 3d and 4d plots,
                one or a sequence of up to five values "log" or "lin"
                which are specified in the order x, y, z, c, and yr.
          x_axis_scale, y_axis_scale, yr_axis_scale, z_axis_scale,
                 and c_axis_scale are used to change the scale on
                 individual axes. They may be set to "log" or "lin".
       text = <value> where value is one or a sequence of strings
                representing texts to be placed on the plot.
       text_color = <value> where <value> is one or a sequence
                of color numbers or names giving colors for the
                texts.
       text_size = <value> where <value> is one or a sequence of
                integers giving (roughly) the number of characters
                in a line on the graph.
       text_pos = <value> where <value> is a pair or a sequence or
                reals between 0. and 1.0 giving the relative
                position of the lower left corner of a text
                in the graphics window.
       xyequal = 0/1 If 1, the axis limits will be adjusted so that
                both axes are to the same scale.
       sync = 0 or 1 (1 to synchronize before plotting)
                defaults to 1, otherwise plots may get garbled.
       color_bar = 0 or 1 (1 enables plotting of a color bar on
                any graphs for which it is meaningful (colored contour
                plots, filled contour plots, cell arrays, filled
                meshes and polygons).
       color_bar_pos (ignored unless a color bar is actually plotted)
                is a 2d array [ [xmin, ymin], [xmax, ymax]] specifying
                where (in window coordinates) the diagonally opposite
                corners of the color bar are to be placed.
   """

   def _initialize_generics ( self , keywords ) :
      """_initialize_generics ( keywords ) is called from __init__ to
      set all those things which are generic to a graph: sync, 
      axis_labels, axis_scales, axis_limits, etc., and
      titles, title_colors, text, text_color, text_size, and text_pos.
      Gives no default if keyword is not present, so can be called
      by change.
      Formerly, this routine used to do color conversions between
      graphics versions. Not any more; that has been relegated to
      the Plotter level.
      This routine has also been changed to delete generic items from
      the keywords so that other items can be changed as well. It was
      not intuitive why calling 'change' would not change any item
      you pleased.
      """
      if keywords.has_key ("sync") :
         self._sync = keywords ["sync"]
         del keywords ["sync"]
      if keywords.has_key ("titles") :
         self._titles = keywords ["titles"]
         if type (self._titles) == StringType :
            self._titles = [self._titles]
         del keywords ["titles"]
      if keywords.has_key ("title_colors") :
         self._title_colors = keywords ["title_colors"]
         del keywords ["title_colors"]
      if keywords.has_key ("text") :
         self._text = keywords ["text"]
         del keywords ["text"]
      if keywords.has_key ("text_color") :
         self._text_color = keywords ["text_color"]
         del keywords ["text_color"]
      if keywords.has_key ("text_size") :
         self._text_size = keywords ["text_size"]
         del keywords ["text_size"]
      if keywords.has_key ("text_pos") :
         self._text_pos = keywords ["text_pos"]
         del keywords ["text_pos"]
      if keywords.has_key ("tosys") :
         self._tosys = keywords ["tosys"]
         del keywords ["tosys"]
      self._init_axis_labels (keywords)
      # handle keywords which take care of all axes
      if keywords.has_key ("axis_scales") :
         axs = keywords ["axis_scales"]
         if is_scalar (axs) and self._scale_dict_2d.has_key (axs) :
            axs = self._scale_dict_2d [axs]
         if is_scalar (axs) :
            self._axis_scales [0] = axs
         else :
            for i in range (len (axs)) :
               self._axis_scales [i] = axs [i]
         del keywords ["axis_scales"]
      if keywords.has_key ("axis_limits") :
         axl = keywords ["axis_limits"]
         if axl == "defaults" :
            for i in range (self._no_of_axes) :
                self._axis_limits [i] = [0., 0.]
         elif is_scalar (axl) :
            raise self._AxisSpecError , \
               "Axis limits must be a point."
         elif type (axl) == ListType or type (axl) == ArrayType :
            if type (axl [0]) != ListType and type (axl [0]) != ArrayType :
               self._axis_limits [0] = axl
            else :
               for i in range (len (axl)) :
                  self._axis_limits [i] = axl [i]
         else :
            raise   self._AxisSpecError , \
               "Axis limits are incomprehensible."
         del keywords ["axis_limits"]
      if keywords.has_key ("color_bar") :
         self._color_bar = keywords ["color_bar"]
      else :
         self._color_bar = 0
      if keywords.has_key ("color_bar_pos") :
         self._color_bar_pos = keywords ["color_bar_pos"]
      else :
         self._color_bar_pos = None
      # handle individual axis keywords
      for i in range (self._no_of_axes) :
         if keywords.has_key (self._scale_keywords [i]) :
            self._axis_scales [i] = keywords [self._scale_keywords [i]]
            del keywords [self._scale_keywords [i]]
         if keywords.has_key (self._limits_keywords [i]) :
            self._axis_limits [i] = keywords [self._limits_keywords [i]]
            del keywords [self._limits_keywords [i]]
      if keywords.has_key ( "xyequal") :
         self._xyequal = keywords ["xyequal"]
         del keywords ["xyequal"]

   # common to all Graph classes
   _GraphInitError = "GraphInitError"
   _GraphSpecError = "GraphSpecError"
   _AxisSpecError = "AxisSpecError"

   _scale_dict_2d = { "linlin" : ["lin", "lin"] ,
                      "linlog" : ["lin", "log"] ,
                      "loglin" : ["log", "lin"] ,
                      "loglog" : ["log", "log"] }

   def __init__ ( self , * kwds , ** keywords ) :
      global Graphics
      self._send_coordinates = 1
      if len (kwds) == 1 :
         keywords = kwds[0]
      if not hasattr (self, "_filename_list") :
         self._filename_list = []
      if not hasattr (self, "_plotter_list") :
         self._plotter_list = []
      if ( keywords.has_key ("filename") or keywords.has_key ("display")) \
         and keywords.has_key ("plotter") :
         raise self._GraphInitError , \
            "Only one keyword, 'filename (display)' or 'plotter' allowed."
      if keywords.has_key ("hardcopy") :
         self._hardcopy = keywords ["hardcopy"]
      else :
         self._hardcopy = ""
      self._graphics_list = [Graphics]
      if keywords.has_key ("filename") or keywords.has_key ("display") :
         if keywords.has_key ("graphics") :
            if keywords ["graphics"] == "" or keywords ["graphics"] is None :
                self._graphics_list = [Graphics]
            elif keywords ["graphics"] == "Nar" or \
                 keywords ["graphics"] == "Narcisse" :
                self._graphics_list = [NarPlotter]
            elif keywords ["graphics"] == "Gist" :
                self._graphics_list = [GistPlotter]
            else :
                self._graphics_list = keywords ["graphics"]
            if is_scalar (self._graphics_list) :
               self._graphics_list = [self._graphics_list]
            for i in range (len (self._graphics_list)) :
               if self._graphics_list [i] == "Gist" :
                  self._graphics_list [i] = GistPlotter
               elif self._graphics_list [i] == "Nar" or \
                    self._graphics_list [i] == "Narcisse" :
                  self._graphics_list [i] = NarPlotter
               else :
                  raise self._GraphInitError , "Illegal graphics option: " + \
                      self._graphics_list [i]
         if keywords.has_key ("filename") :
            fn = keywords ["filename"]
         else :
            fn = keywords ["display"]
         if is_scalar (fn) and not fn in self._filename_list :
            self._filename_list.append (fn)
            self._plotter_list.append (-1)
         else :
            for i in range (len (fn)) :
               if not fn [i] in self._filename_list :
                  self._filename_list.append (fn [i])
                  self._plotter_list.append (-1)
         self._plotters_started = 0
      elif keywords.has_key ("plotter") :
         pl = keywords ["plotter"]
         if is_scalar (pl) and not pl in self._plotter_list :
            self._plotter_list.append (pl)
            self._filename_list.append (" ")
         else :
            for i in range (len (pl)) :
               if not pl [i] in self._plotter_list :
                  self._plotter_list.append (pl [i])
                  self._filename_list.append (" ")
         self._plotters_started = 1
      self._titles = [" ", " ", " ", " "]
      self._title_colors = "fg"
      self._text = " "
      self._text_color = 1
      self._text_size = 0
      self._text_pos = [0.,0.]
      self._tosys = 0
      self._sync = 1
      self._xyequal = 0
      self._initialize_generics ( keywords )

   def _init_axis_labels ( self , keywords ) :
      """_init_axis_labels ( kw ) sets the axis labels if keywords are present.
      This routine can be called by quick_plot in a dervied class.
      """
      if keywords.has_key ("axis_labels") :
         axl = keywords ["axis_labels"]
         if is_scalar (axl) :
            self._axis_labels [0] = axl
         else :
            for i in range (len (axl)) :
               self._axis_labels [i] = axl [i]
         del keywords ["axis_labels"]
      if keywords.has_key ("gnomon") :
         self._gnomon = keywords ["gnomon"]
         del keywords ["gnomon"]
      else :
         self._gnomon = 0
      for i in range (self._no_of_axes) :
         if keywords.has_key (self._label_keywords [i]) :
            self._axis_labels [i] = keywords [self._label_keywords [i]]
            del keywords [self._label_keywords [i]]
 
   def add_file ( self, fn ) :
       """add_file ( "filename" ) allows the user to add a plotter
       contacted via "filename" to the list of plotters being used
       by this object.
       """
       self._filename_list.append (fn)
       self._plotter_list.append (-1)
       self._plotters_started = 0

   add_display = add_file
 
   def delete_file ( self, fn ) :
       """delete_file ( "filename" ) allows the user to delete a plotter
       contacted via "filename" from the list of plotters being used
       by this object.
       """
       if fn in self._filename_list :
          self._filename_list.remove (fn)
       return

   delete_display = delete_file
 
   def add_plotter ( self, pl ) :
       """add_plotter ( pl ) allows the user to add the specified
       plotter to the list of plotters being used by this object.
       """
       self._plotter_list.append (pl)
       self._filename_list.append (" ")
 
   def delete_plotter ( self, pl ) :
       """delete_plotter ( pl ) allows the user to delete the specified
       plotter from the list of plotters being used by this object.
       """
       if pl in self._plotter_list :
          self._plotter_list.remove (pl)
       return
 
   def change ( self, * kwds , ** keywords ) :
       """change ( <keyword arguments> ) allows some of the graph's
       generic characteristics to be changed.
       Actually it allows any keyword values to be changed. It should
       be used with care, as it does no error checking.
       """
       if len (kwds) == 1 :
          keywords = kwds[0]
       self._initialize_generics ( keywords )
       for k in keywords.keys ( ) :
          setattr ( self, "_" + k , keywords [k] )

   def _init_plotter ( self ) :
       """ _init_plotter (): initializer routine which will start up
       plotter(s) if necessary. This routine is called from the plot
       routine of a derived class.
       """
       try :
          Graphics = self._graphics_list [0]
       except :
          raise self._GraphInitError, "No graphics specified??"
       if self.type () == Graph2dType :
          dpi = 100
       else :
          dpi = 75
       if len (self._plotter_list) == 0 :
          # If user has not specified a plotter, create a default
          self._filename_list.append (" ")
          self._plotter_list.append (Graphics.Plotter (" " ,
             style = self._style, dpi = dpi, hcp = self._hardcopy))
       else :
          for i in range (len (self._plotter_list)) :
              # start up plotters for any unspecified ones.
              if self._plotter_list [i] == -1 :
                 if len (self._plotter_list) > 1 :
                    Graphics = self._graphics_list [i]
                 self._plotter_list [i] = Graphics.Plotter (
                    self._filename_list [i] , style = self._style, dpi = dpi,
                    hcp = self._hardcopy)
       self._plotters_started = 1
 
