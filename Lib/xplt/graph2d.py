# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

from graph import *
from graftypes import *
 
class Graph2d ( Graph ) :
 
   """
   g = Graph2d ( <curve list>, ...keyword arguments...) will create
   a two-dimensional graphics object consisting of several curves
   plus a global environment for the object. It will accept one or
   a list of Plotter objects or plotter_identifiers, or will try
   to complete a generic connection of its own if asked to plot
   without such a plotter specification.
   <curve list> is one or a sequence of Curve, Lines, QuadMesh, 
   CellArray, or PolyMap objects.
   The keyword arguments for Graph2d are:
 
       plotter = <Plotter object> or a sequence of <Plotter object>s
              if you have a plotter object or a sequence of them that
              you want the curve to use when plotting itself. I
              recommend against this; a curve can create its own
              plotter if you don't give it one.
       filename = <string> or a sequence of <string>s if you want to
              connect to Narcisse with the named file(s). (The default
              is " ".) Only one of the keywords 'plotter' or 'filename'
              is allowed, and both are optional.
              NOTE: the possibility of a sequence of file names or
              plotters allows one to display Narcisse graphs on one's
              own machine as well as one or more remote machines.
              In Gist, the filename argument is the same as the
              display argument (below).
       display = <string> or a sequence of <string>s if you want to
              display on the named hosts. The form of <string> is
              the usual "hostname:server.screen". The purpose of
              this argument is to allow you to continue without
              exiting python if you have to open a Gist window
              without the DISPLAY environment variable having been
              set, or if you want to open Gist windows on more
              than one host.
       style = one of "vg.gs", "boxed.gs", "vgbox.gs", "nobox.gs",
              "work.gs". For Gist only, and only if a plotter has
              not been specified, each plotter opened will be opened
              with the specified style. The default is "work.gs".
       grid_type = <string> where "none" means no axis grid;
                "axes" means a pair of axes with tick marks;
                "wide" means a widely spaced 2d grid; and
                "full" means a closely spaced 2d grid.
       label_type = "end" (to label the curve at its end)
                    "box" (to put the labels in a box)
       titles = <value> where <value> is a string or a sequence
                of up to four strings, giving the titles in the
                order bottom, top, left, right.
       title_colors = <value> where value is an integer or string
                or a sequence of up to four integers or strings giving
                the colors of the titles.
       grid_type = <string> where "none" means no axis grid;
                "axes" means a pair of axes with tick marks;
                "wide" means a widely spaced 2d grid; and
                "full" means a closely spaced 2d grid.
       axis_labels = <value> where <value> is a string or
                sequence of up to three strings representing
                the labels of the x axis, the y axis, and the
                right y axis.
       x_axis_label, y_axis_label, and yr_axis_label may be used
                to label individual axes.
       axis_limits = <value> where <value> is a pair [xmin, xmax] or
                a sequence of up to three pairs, where the second
                would be the y limits, and the third the yr limits.
       x_axis_limits, y_axis_limits, and yr_axis_limits may be used
                to specify limits on individual axes.
       axis_scales = "linlin", "linlog", "loglin", or "loglog"
                or, if all three axes are to be specified, a
                triple of the values "lin" and "log".
       x_axis_scale, y_axis_scale, and yr_axis_scale may be used
                to specify individual axis scales.
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
       color_card = <value> specifies which color card you wish
                to use, e. g., "rainbowhls" (the default), "random",
                etc. Note that for curves, color_card is a Graph2d
                keyword, since it is not possible to specify
                different color cards on the same 2d graph,
                whereas linked 3d and 4d graphs can have
                different color cards; therefore for these,
                it is a Surface keyword.
       sync = 0 or 1 (1 to synchronize before sending a plot)
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
 
   _axes = ["x", "y", "yr"]
 
   _NotA2dObject = "NotA2dObject"

   # The following are actually used in Graph
   _limits_keywords = ["x_axis_limits", "y_axis_limits", "yr_axis_limits"]
   _scale_keywords = ["x_axis_scale", "y_axis_scale", "yr_axis_scale"]
   _label_keywords = ["x_axis_label", "y_axis_label", "yr_axis_label"]
   _no_of_axes = 3

   def type (self) :
      return Graph2dType
 
   def __init__ ( self , curve_list , *kwds , ** keywords ) :
      if len ( kwds ) == 1 :
         keywords = kwds[0]
      if is_scalar ( curve_list ) :
         self.check_curve (curve_list)
         self._c = [curve_list]
      else :
         for i in range ( len ( curve_list ) ) :
            self.check_curve (curve_list [i])
         self._c = curve_list
      self._c_ln = len (self._c)
      if keywords.has_key ("label") :
         self._label = keywords ["label"]
      else :
         self._label = " "
      if keywords.has_key ("label_type") :
         self._label_type = keywords ["label_type"]
      else :
         self._label_type = " "
      if keywords.has_key ( "color_card" ) :
         self._color_card = keywords ["color_card"]
      else :
         self._color_card = "default"
            # Everything else is Graph generic:
      self._axis_limits = [[0., 0.], [0., 0.], [0., 0.]]
      self._axis_scales = ["lin", "lin", "lin"]
      self._axis_labels = ["X axis", "Y axis", "YR axis"]
      Graph.__init__ ( self , keywords )
      if keywords.has_key ("grid_type") :
         self._grid_type = keywords ["grid_type"]
      else :
         self._grid_type = "axes"
      if keywords.has_key ("style") :
         self._style = keywords ["style"]
      else :
         self._style = "work.gs"

   def new ( self , curve_list , ** keywords ) :
      """new ( curves, <keyword arguments> ) cleans out a Graph2d
      and reinitializes it. This has the same argument list as
      Graph2d. Do not change the plotter list or filename list (to
      avoid the tedious on-and-off flickering of windows) unless
      this is actually requested.
      """
      del self._c
      pl = self._plotter_list
      fl = self._filename_list
      self._plotter_list = []
      self._filename_list = []
      self.__init__ ( curve_list , keywords )
      if self._plotter_list == [] :
         self._plotter_list = pl
         self._filename_list = fl
 
   def set ( self , ** keywords ) :
       """ set (...keyword arguments...) allows you to set individual
       Graph2d characteristics. No error checking is done.
       It will only change the plotter if specifically asked to do so.
       """
       self._filename = ""
       self._display = ""
       self._plotter = None
       for k in keywords.keys ():
          setattr (self, "_" + k, keywords [k])
       if is_scalar (self._filename) and self._filename != "" :
          self._filename_list = [self._filename]
          self._plotter_list = [-1]
       elif is_scalar (self._display) and self._display != "" :
          self._filename_list = [self._display]
          self._plotter_list = [-1]
       elif self._filename != "" :
          self._filename_list = self._filename
          self._plotter_list = [-1] * len (self._filename_list)
       elif self._display != "" :
          self._filename_list = self._display
          self._plotter_list = [-1] * len (self._filename_list)
       elif self._plotter is not None :
          if is_scalar (self._plotter) :
             self._plotter_list = [self._plotter]
             self._filename_list = [" "]
          else :
             self._plotter_list = self._plotter
             self._filename_list = [" "] * len (self._plotter_list)
          
 

   def check_curve (self, crv) :
      """check_curve (crv) raises an exception if crv is not a legal
      2d object.
      """
      try :
         dum = crv.type ()
      except:
         raise self._NotA2dObject , "Illegal object sent to Graph2d."
      if crv.type () != CurveType and \
         crv.type () != LinesType and \
         crv.type () != PolyMapType and \
         crv.type () != CellArrayType and \
         crv.type () != QuadMeshType and \
         crv.type () != Animation2dType :
         raise self._NotA2dObject , "Illegal object (type " + \
            `crv.type ()` + ") sent to Graph2d."
      return
 
   def add ( self , curve ) :
      """ add ( curve ) adds a curve with its characteristics to the
      existing plot. Curves are numbered in the order that they are
      added, beginning with 1.
      """
      self.check_curve (curve)
      self._c.append ( curve )
      self._c_ln = len (self._c)
 
   def delete ( self , n ) :
      """delete ( n ) n integer: deletes the nth curve from the Graph. Note
      that curves are numbered beginning with 1.
      n a curve (etc.) object: deletes that object, if it is present.
      """
      DeleteError = "DeleteError"
      if type (n) == IntType and 1 <= n <= self._c_ln :
         self._c [n-1:n] = []
         self._c_ln = len (self._c)
      elif type (n) == IntType :
         raise DeleteError , "There is no curve numbered " + `n` + \
            " in the current graph, which has only " + `self._c_ln` + \
            " curves."
      else :
         for i in range (self._c_ln) :
            if n == self._c [i] :
               del self._c [i]
               self._c_ln = self._c_ln - 1
               return
 
   def replace ( self , n , curve ) :
      """replace ( n , curve ) n integer : replaces the nth curve in the Graph
      with the specified curve. Note that curves are numbered beginning with 1.
      n a curve (etc.) object: replaces that object, if it is present.
      """
      self.check_curve (curve)

      ReplaceError = "ReplaceError"
      if IntType == type (n) and 1 <= n <= self._c_ln :
         self._c [n-1] = curve
      elif IntType == type (n) :
         raise ReplaceError , "There is no curve numbered " + `n` + \
            " in the current graph, which has only " + `self._c_ln` + \
            " curves."
      else :
         for i in range (self._c_ln) :
            if n == self._c [i] :
               self._c [i] = n
               return
         else :
            raise ReplaceError , "Attempt to remove an object not on the list."
 
   _CurveChangeError = "CurveChangeError"

   def change_curves (self, curve_or_list) :
      """ change_curves ( <curve or curve list> ) is designed to replace
      the entire curve list.
      """
      if is_scalar ( curve_or_list ) :
         self.check_curve (curve_or_list)
         self._c = [curve_or_list]
      else :
         for i in range ( len ( curve_or_list ) ) :
            self.check_curve (curve_or_list [i])
         self._c = curve_or_list
      self._c_ln = len (self._c)
      return
 
   def change_plot ( self , ** keywords ) :
      """change_plot ( <keyword arguments> ) is used to change any Graph2d
      characteristics except the curves being graphed. Use the add,
      delete, and/or replace commands to do that. change_plot will
      draw the graph without sending object surface coordinates, unless
      keyword send is 1. Generally, change_plot should be used when
      the graph needs to be recomputed, and quick_plot when it does not.
      change_plot does no error checking and does not conver user-friendly
      names of colors and such into numbers.
      change_curves ( <curve or curve list> ) may be used to replace
      thr rntire curve list.
      """
      for k in keywords.keys ():
         if k == "curve" :
            raise self._CurveChangeError, \
               "Use add, delete, or replace to change curves in a graph."
         setattr (self, "_" + k, keywords [k])
      if "send" in keywords.keys ():
         send = keywords ["send"]
         self._send_coordinates = send
      self.plot ( )
      self._send_coordinates = 1
 
   def quick_plot ( self , ** keywords ) :
      """quick_plot ( <keyword arguments> ) is used to change some Graph2d
      characteristics which do not demand that the graph be recomputed.
      You can change the characteristics of a curve in the graph by
      specifying its number (curve = n) and any combination of the
      traits type, color, and label. Or you can change such overall
      graph characteristics as label_type, titles, title_colors,
      text, text_color, text_size, text_pos, color_card, grid_type,
      sync, and axis_labels. The changes will be effected and the graph
      redrawn.
      Things that you cannot change include axis limits and scales,
      and the coordinates of a curve. Use change_plot if axis limits
      and scales are among the things you want to change, and use add,
      delete, or replace followed by a call to plot, if you wish to
      change a curve.
      quick_plot will tell you if you try to change something illegal.
      """
      ChangeError = "ChangeError"
      PlottersNotStarted = "PlottersNotStarted"
      if not self._plotters_started :
         raise PlottersNotStarted , \
            "quick_plot requires that all plotters have already been started."
      if keywords.has_key ( "curve" ) :
         n = keywords ["curve"]
         del keywords ["curve"]
         self.type_change = 0
         self.color_change = 0
         self.label_change = 0
         if 1 <= n <= self._c_ln :
            if keywords.has_key ("type") :
               self._c[n - 1].line_type = keywords ["type"]
               del keywords ["type"]
               self.type_change = 1
            if keywords.has_key ("color") :
               self._c[n - 1].color = keywords ["color"]
               del keywords ["color"]
               self.color_change = 1
            if keywords.has_key ("label") :
               self._c[n - 1].label = keywords ["label"]
               self.label_change = 1
               del keywords ["label"]
         else :
            raise ChangeError , "There is no curve numbered " + `n` + \
            " in the current graph, which has only " + `self._c_ln` + \
            " curves."
      if keywords.has_key ("label_type") :
         self._label_type = keywords ["label_type"]
         del keywords ["label_type"]
      if keywords.has_key ( "color_card" ) :
         self._color_card = keywords ["color_card"]
         del keywords ["color_card"]
      if keywords.has_key ("grid_type") :
         self._grid_type = keywords ["grid_type"]
         del keywords ["grid_type"]
      Graph.change ( self , keywords ) # Change generic traits
      if len (keywords.keys ()) > 0 :
         print "Note: quick_plot will ignore keywords" , keywords.keys () , "."
         print "-- Use change_plot instead."
      for ipl in range (len (self._plotter_list)) :
         pl = self._plotter_list [ipl]
         try :
            dum = Nar
         except :
            self.plot (pl)
         else :
            if self._graphics_list [ipl] == Nar :
               pl.quick_plot (self)
            else :
               self.plot (pl)

   def plot ( self , plotter = None) :
      """plot ( ) plots a 2d graph object. If the user has not by
      now specified plotter(s) or filename(s) then a generic plotter
      object will be created, if it is possible to find a local
      Graphics routine.
      """
      self._init_plotter ( )
      if plotter is not None :
         plotter.plot2d (self)
         return
      for ipl in range (len (self._plotter_list)) :
         pl = self._plotter_list [ipl]
         pl.plot2d (self)
