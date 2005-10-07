# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

from Numeric import *
from scipy_base.fastumath import *
from mesh3d import Mesh3d
from graph3d import Graph3d
 
def paws ( ) :
    i = raw_input ("Type in any string to continue; ^C to return to prompt. ")
    return
 
print "This is a test of the Python interface to the Limeil Lab graphics"
print "package, Narcisse. You  need Narcisse to be running. Fire it up by"
print "doing setenv PORT_SERVEUR 0, typing  /dist/basis/Narcisse/bin/Narcisse &,"
print "and then do another senetv PORT_SERVEUR to the port number which"
print "appears in the Narcisse GUI."
 
ans = raw_input ("The next graph takes a *looong* time. Want to graph it? (y/n)")

def demo () :
   from PR import PR
   f = PR ('./berts_plot')
   x = f.x
   y = f.y
   z = f.z
   c = f.c
   c_cont = f.c_cont
   m1 = Mesh3d (x = x, y = y, z = z, c = c, opt_3d = "i4",
                  c_contours_array = c_cont, mask = "none")
   g1 = Graph3d (m1, titles = "Energy Filaments, cutoff = 7.5",
                  phi = 60., theta = 45.)
   g1.plot ( )
   f.close ( )

if ( ans == "y" or ans == "Y" or ans == "yes" or
     ans == "Yes" or ans == "YES" ) :
   demo ()
