# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

# ed williams' colorbar stuff
from Numeric import *
from arrayfns import *
from gist import *
from slice3 import *

def nice_levels (z, n = 8) :
   """nice_levels(z, n = 8) finds approximately n "nice values"
   between min(z) and max(z) for axis labels. n defaults to eight.
   """
   zmax = max (ravel (z))
   zmin = min (ravel (z))
   finest = abs (zmax - zmin) / float (n)
   # blows up on zmin=zmax
   unit = 10. ** floor (log10 (finest))
   finest = finest / unit
   if finest > 5.0 :
      finest = 10.
   elif finest > 2. :
      finest = 5.
   elif finest > 1. :
      finest = 2.
   unit = unit * finest
   cmin = unit * ceil (zmin / unit)
   if (abs (cmin - zmin) < 0.01 * unit) :
      cmin = cmin + unit
   cmax = unit * floor (zmax / unit)
   if (abs (cmax - zmax) < 0.01 * unit) :
      cmax = cmax - unit
   n = int ( ( (cmax - cmin) / unit + 0.5) + 1)
   levs = span (cmin, cmax, n)
   list = nonzero (less (abs (levs), 0.1 * unit))
   if len (list) > 0 :
      array_set (levs, list, 0.0)
   return levs

def color_bar (minz, maxz, split = 0, ncol = None) :
   """
   color_bar (minz, maxz) plots a color bar to the right of the plot square
   labelled by the z values from minz to maxz.

   plf (z, y, x)
   color_bar (z (min, min), z (max, max))

   or
   plf (z, y, x, cmin = MINZ, cmax = MAXZ)
   color_bar (MINZ, MAXZ)

   are typical usage
   """
   if ncol == None:
      ncol = 100 + (1 - split) * 100
   plsys (0)
   ymax = 0.85
   ymin = 0.44
   if type (minz) == type (maxz) == type (1) : # Do not change!!!
      plotval = reshape (arange (minz, maxz + 1, typecode = 'b'),
         (maxz + 1 - minz, 1))
      pli (plotval, .60, ymin, .62, ymax) # draw bar
   elif not split :
      pli (reshape (span (0, 1, ncol), (ncol, 1)),
         .60, ymin, .62, ymax) # draw bar
   else :
      pli (reshape (split_bytscl (span (0, 1, ncol), 0).astype ('b'), (ncol, 1)),
         .60, ymin, .62, ymax) # draw bar
   pldj (array ( [.60, .60]), array ( [ymin, ymax]), array ( [.62, .62]),
         array ( [ymin, ymax]))
   plsys (1)
   levs = nice_levels (array ( [minz, maxz]))
   scales = []
   for i in range (len (levs)) :
      scales.append ( "% .5g" % levs [i])
   ys = ymin + (ymax - ymin) * (levs - minz) / (maxz - minz)
   llev = len (levs)
   rllev = range (llev)
   for i in rllev :
      plt (scales [i], .625, ys [i])   # labels
   xmin = zeros (llev, Float)
   xmax = zeros (llev, Float)
   xmin [0:llev] = .60
   xmax [0:llev] = .625
   plsys (0)
   pldj (xmin, ys, xmax, ys) # ticks
   plsys (1)
   # Write the max and min on bar
   if max (ys) > 0.84 :
      plt ("% .5g" % maxz, .59, ymax + 0.020)
   else :
      plt ("% .5g" % maxz, .59, ymax + 0.005)

   plt ("% .5g" % minz, .59, ymin - 0.005, justify = "LT")
