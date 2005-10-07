# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

from lines import *
from polymap import *
from cellarray import *
from graph2d import *

def paws ( ) :
    i = raw_input ("Type in any string to continue; ^C to return to prompt. ")
    return

def demo () :
   x0 = arange(50, typecode = Float)
   y0 = zeros(50, Float)
   x1 = 49 * ones(50, Float)
   y1 = arange(50, typecode = Float)
   ly = Lines (x0 = x0,y0 = y0,x1 = x1,y1 = y1)
   g0 = Graph2d ( ly , titles = "Just Lines")
   g0.plot ()
   paws ()
   ly.set (color = "red")
   g0.change (titles = "Lines colored red")
   g0.plot ()
   paws ()
   ly.set (width=4.0,type="dashdotdot")
   g0.change(titles = "Wide lines, dashdotdot style")
   g0.plot ()
   paws ()
   n = [4,3,3,3,3,4,4,3,3,3,3]
   x = [0.,1.,1.,0.,0.,1.5,1.,1.5,3.,0.,1.5,3.,2.,1.5,2.,1.,2.,3.,3.,2.,
        1.,2.,2.,1.,2.,3.,1.5,1.,2.,1.5,1.,1.5,0.,0.,3.,1.5]
   y = [1.,2.,7.,8.,1.,1.,2.,0.,1.,1.,1.,1.,2.,1.,2.,2.,2.,1.,8.,7.,
        2.,2.,7.,7.,7.,8.,8.,7.,7.,8.,7.,8.,8.,8.,8.,9.]
   z = [2.5,1.2,1.5,1.2,.5,2.5,2.,1.2,.5,1.2,1.5]
   p1 = Polymap(x=x,y=y,z=z,n=n)
   g0.replace(1,p1)
   g0.change(titles = "Stained glass window")
   g0.plot ()
   paws ()
   p1.set(z = array("AMZMmAzMmMZ",'b') )
   g0.plot ()
   paws ()
   nx= 10
   ny= 19
   # ndx = zeros ( (nx, ny ), Int)
   # for jj in range (ny) :
   #    for ii in range (nx) :
   #        ndx [ii, jj] = ii + (jj - 1) * nx
   ndx = reshape (arange (nx * ny, typecode = 'b'), (nx, ny))
   cla = CellArray ( z = ndx )
   gca = Graph2d ( cla , titles = "Cell Array", axis_scales = "linlin" )
   gca.plot ( )
   paws ( )
   # -------------------- generate sombrero function
   nz = 20
   x = arange (-nz, nz+1, typecode = Float )
   y = x
   z = zeros ((2*nz + 2, 2*nz + 2), Float)
   for i in range ( len (x) ) :
      for j in range ( len (y)) :
         r = sqrt ( x [i] * x[i] + y [j] * y [j] ) + 1.e-12
         z [i, j] = sin (r) / r
   # cell array plot
   cla.set (z = z)
   gca.change (titles = "Sombrero Function",axis_limits="defaults")
   gca.plot ( )
   paws ( )

print "Type demo() to proceed."
