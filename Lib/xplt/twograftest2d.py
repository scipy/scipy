# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

from Numeric import *
from scipy_base.fastumath import *
GraphicsError = "GraphicsError"

import os
try:
   graphics = os.environ["PYGRAPH"]
except KeyError:
   graphics = "Gist"

if graphics [0:3] == "Nar" :
   import NarPlotter
   gr = NarPlotter
   print "This is a test of the Python interface to the Limeil Lab graphics"
   print "package, Narcisse. You  need Narcisse to be running. Fire it up by"
   print "doing setenv PORT_SERVEUR 0, typing  /dist/basis/Narcisse/bin/Narcisse &,"
   print "and then do another senetv PORT_SERVEUR to the port number which"
   print "appears in the Narcisse GUI."
elif graphics == "Gist" :
   import GistPlotter
   gr = GistPlotter
   print "This is a test of the Python interface to the Gist graphics"
   print "package."
else :
   raise GraphicsError , \
      graphics + " is an unknown graphics package. Check PYGRAPH " + \
         "environment variable."

print "Test of 2D plotting routine. Type demo () to see the test."
paws ()

from curve import Curve
from graph2d import *
def paws ( ) :
    i = raw_input ("Type in any string to continue; ^C to return to prompt. ")
    return
 
def demo () :
    
   vsf = 0.
   c = 1
   s = 1000.
   kmax = 25
   lmax = 35
    
   pl1 = gr.Plotter (" ")
   c1 = Curve ( y = arange (1, kmax+1, typecode = Float) , color = "yellow" )
   g2 = Graph2d ( c1 , plotter = [pl1], titles = ["Bottom of the Barrel", "Top Dog",
                                 "Leftist", "Reaganist"] )
   g2.plot ( )
   paws ( )
   c2 = Curve ( y = sqrt (arange (1, kmax+1, typecode = Float)**3) , color = "blue")
   g2.add ( c2 )
   g2.plot ( )
   paws ( )
   c1.set (marks = 1, marker = "A")
   c2.set (marks = 1, marker = "B")
   g2.plot ( )
   paws ( )
   c1.set (marks = 0)
   c2.set (marks = 0)
   c1.set (width = 3.0)
   c2.set (width = 4.0)
   g2.plot ( )
   paws ( )
   c1.set (type = "dash")
   c2.set (type = "dashdot")
   g2.plot ( )
   paws ( )
   g2.delete ( 1 )
   g2.plot ( )
   paws ( )
   mycolor= ["red", "blue", "yellow", "green", "orange", "purple"]
   markers= ["+", "*", "o", "x", "."]
   c5 = [ ]
   for i in range (len (markers)) :
      c5.append ( Curve (
         y = exp (array (multiply ((1.0+(i+1)/10.), log (range (1, kmax))))) ,
                  color = mycolor [i] , marks = 1, marker = markers [i] ,
                  type = "none" )
                )
   g2.new ( c5 , text="Line style and grid test" , text_color="black" ,
             plotter = [pl1],
             text_pos = array ( [.2, .8]) , text_size = 60 ,
             titles = "Figure 1")
   g2.plot ( )
   paws ( )
   g2.quick_plot ( grid_type = "wide" )
   paws ( )
   part1 = abs (array ( [-5, -4, -3, -2, -1, .0000000001] ) )
   c1.new (  y = part1, x = arange (-5, 1, typecode = Float) ,
            type="line" , color = "blue" )
   part2 = array ( [.0000000001, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10] )
   c2.new ( y = part2 , x = arange (0, 11, typecode = Float) ,
            type = "line" , color = "blue" )
   g2.new ( [c1, c2] , text = " " , axis_scales = "linlog" , titles = " " ,
             old_plotter = 1,
             axis_limits = array ( [[-6., 10.], [1.e-6, 10.]]))
   g2.plot ( )
   paws ( )
   c1.new ( y = array (exp (multiply (log (10.0), range (-50, 1)))) ,
            x = arange (1, 52, typecode = Float) )
   g2.new ( c1 , axis_scales = "loglog" ,
             old_plotter = 1,
                 axis_limits = array ( [[1., 100.], [1.e-8, 1.]]))
   g2.plot ( )
   paws ( )
   c1.new ( y = array ( [1.e-48, 1.e-21, 0., 0., 0.]) )
   g2.new ( c1 , axis_scales = "linlin" ,
             old_plotter = 1,
                 axis_limits = array  ( [[0., 0.], [0., 0.]]))
   g2.plot ( )
