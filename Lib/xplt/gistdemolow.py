# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.
from gist import *
import time

def grtest():
  print "\n\n         Pygist Comprehensive (well, sort of) Graphics Test\n"
  print "Each frame will be described at the terminal."
  print "Compare what you see with the description, then"
  print "hit <RETURN> to see the next test, or q <RETURN> to quit.\n"
  if quitnow(): return
  pldefault(marks=1, width=0, type=1, style="work.gs", dpi=100)
  winkill(0)
  print "Test 1: Commands: window (0, wait=1, dpi=75); plg([0, 1])"
  window (0, wait=1, dpi=75)
  plg([0, 1])
  print "A small (75 dpi) window with line marked A from (1,0) to (2,1)."
  if quitnow(): return
  winkill(0)
  ###
  print "Test 2: Commands: plg([0,1])"
  plg([0,1])
  print "A large (100 dpi) window with line marked A from (1,0) to (2,1)."
  if quitnow(): return
  unzoom()
  ###
  print "Test 3:  Commands: plg([1, 0])"
  plg([1, 0])
  print "Added line marked B from (1,1) to (2,0) to previous plot."
  if quitnow(): return
  unzoom()
  ###
  print "Test 4: Commands: logxy(1, 0)"
  logxy(1, 0)
  print "X axis now a log scale."
  if quitnow(): return
  unzoom()
  ###
  print "Test 5:  Commands: logxy(0, 0)"
  logxy(0, 0)
  print "X axis back to linear scale."
  if quitnow(): return
  unzoom()
  ###
  print "Test 6:  Commands: limits(1.2, 1.8, 0.2, 0.8)"
  limits(1.2, 1.8, 0.2, 0.8)
  print "Limits changed to 1.2<x<1.8, 0.2<y<0.8."
  if quitnow(): return
  unzoom()
  ###
  ylimits(0.4, 0.6)
  print "Test 7: Commands: ylimits(0.4, 0.6)"
  ylimits(0.4, 0.6)
  print "Limits changed to 1.2<x<1.8, 0.4<y<0.6."
  if quitnow(): return
  unzoom()
  ###
  limits()
  print "Test 8: Commands: limits"
  print "Limits back to extreme values (1,0) to (2,1)."
  if quitnow(): return
  unzoom()
  ###
  fma()
  x = 10*pi*arange(200, typecode = Float)/199.0
  plg(sin(x), x)
  print "Test 9:  Commands: fma(); plg(sin(x), x)"
  print "where x = 10*pi*arange(200, typecode = Float)/199.0"
  print "Five cycles of a sine wave on a new frame."
  print "Before you continue, try clicking with the mouse buttons:"
  print "Left button zooms in, right button zooms out, middle no zoom"
  print "In each case, the point where you press the mouse button will"
  print "be translated to the point where you release the mouse button."
  print "To zoom one axis only, click over the tick marks on that axis."
  sts = quitnow()
  unzoom()
  if (sts): return
  ###
  pledit(marks=0, width=6, type="dash")
  print "Test 10     Commands: pledit(marks=0, width=6, type=\"dash\")"
  print "Marker A on sine curve disappears, curve becomes bold dashed."
  if quitnow(): return
  unzoom()
  ###
  fma()
  x = 2*pi*arange(200, typecode = Float)/199.0
  for i in range(1,7):
    r = 0.5*i - (5-0.5*i)*cos(x)
    s = `i`
    plg( r*sin(x), r*cos(x), marks=0, color=-4-i, legend = s)
  print "Test 11: Commands: plg(r*sin(x), r*cos(x), color=-4-i)"
  print "A set of nested cardioids in the primary and secondary colors."
  if quitnow(): return
  unzoom()
  ###
  plt ("Colored Nested Cardioids", -2.0, -3.0, orient = 1, opaque = 1, tosys = 1)
  print "Test 12: Commands: pltitle(\"Colored nested cardioids\"); plq()"
  plq()
  print "Adds the title above the upper tick marks."
  print "Also prints legends for the six curves at the terminal (plq)."
  if quitnow(): return
  unzoom()
  ###
  print "Test 13     Commands: pledit(color=\"fg\", type=0, marker=i)"
  for i in range(1,6): pledit(i, color="fg", type=0, marker=i)
  pledit(i, color="fg", type=0, marker='A')
  print "Changes the colors to foreground, types to no lines."
  print "Markers are point, plus, asterisk, O, X, A."
  if quitnow(): return
  unzoom()
  ###
  print "Test 14: Commands: pledit(marks=0, type=i)"
  for i in range(1,6): pledit(i, marks=0, type=i)
  pledit(i, color="fg", type=1, width=4)
  print "Changes line types to solid, dash, dot, dashdot, dashdotdot."
  print "Outermost cardioid becomes a thick, solid line."
  if quitnow(): return
  unzoom()
  ###
  fma()
  limits()
  x = a3(-1, 1, 26)
  y = transpose (x)
  z = x+1j*y
  z = 5.*z/(5.+z*z)
  xx = z.real
  yy = z.imaginary
  print "Test 15     Commands: plm(y, x)"
  print "Quadrilateral mesh -- round with bites out of its sides."
  plm (yy, xx)
  if quitnow(): return
  fma()
  unzoom()
  ###
  plmesh(yy, xx)
  print "Test 16     Commands: plv(v, u, y, x)"
  plv(x+.5, y-.5)
  print "Velocity vectors.  Try zooming and panning with mouse."
  if quitnow(): return
  fma()
  unzoom()
  ####
  print "Test 17     Commands: plc(z, y,x); plm(y,x, boundary=1, type=2)"
  ireg = ones (xx.shape, Int)
  ireg [0, :] = 0
  ireg [:, 0] = 0
  plfc(mag(x+.5,y-.5),yy,xx,ireg,contours=8)
  plc (mag(x+.5,y-.5), marks=1)
  plm (boundary=1, type=2)
  print "Contours A-H, with mesh boundary dashed."
  if quitnow(): return
  fma()
  unzoom()
  ####
  #print "Test 18     Commands: plf(zncen(z), y,x);  plc(z, y,x)"
  print "Test 18     Commands: plf(z, y,x);  plc(z, y,x)"
  z = mag(x+.5,y-.5)
  #plf(zncen(z))
  plf(z)
  fma()
  unzoom()
  if quitnow(): return
  plfc(z,yy,xx,ireg,contours=array([0.5, 1.0, 1.5]))
  plc(z, marks=0, type=2, color="bg", levs=[0.5, 1.0, 1.5])
  print "Filled mesh (color only) with three dashed contours overlayed."
  if quitnow(): return
  fma()
  unzoom()
  ####
  print "Test 19     Commands: palette( \"<various>.gp\")"
  print "After each new palette is installed, hit <RETURN> for the next,"
  print "or q<RETURN> to begin test 20.  There are 6 palettes altogether."
  z = mag(x+.5,y-.5)
  #plf(zncen(z))
  plf(z)
  plc(z, marks=0, type=2, color="bg", levs=[0.5, 1.0, 1.5])
  if quitnow(): return
  fma ()
  plfc(z,yy,xx,ireg,contours=array([0.5, 1.0, 1.5]))
  plc(z, marks=0, type=2, color="bg", levs=[0.5, 1.0, 1.5])
  pal = ["heat.gp", "stern.gp", "rainbow.gp", "gray.gp", "yarg.gp", "earth.gp"]
  for i in range(6):
    palette (pal[i])
    print "Palette name: ", pal[i]
    if quitnow(): return
  ###
  palette("earth.gp")
  fma()
  unzoom()
  ###
  print "Test 20     Commands: window, style=\"<various>.gs\""
  print "After each new style is installed, hit <RETURN> for the next,"
  print "or q<RETURN> to begin test 21.  There are 5 styles altogether."
  pal= ["vg.gs", "boxed.gs", "vgbox.gs", "nobox.gs", "work.gs"]
  for i in range(5):
    window (style=pal[i])
    plc (mag(x+.5,y-.5), marks=1)
    plm (boundary=1, type=2)
    print "Style name: ", pal[i]
    if quitnow(): return
  ###
  window(style="work.gs")
  fma()
  unzoom()
  ####
  print "Test 21     Commands: pli(image)"
  x = a3 (-6,6,200)
  y = transpose (x)
  r = mag(y,x)
  theta = arctan2 (y, x)
  funky = cos(r)**2 * cos(3*theta)
  pli(funky)
  print "Cell array image (color only).  Three cycles in theta, r."
  if quitnow(): return
  fma()
  unzoom()
  ####
  print "Test 22     Commands: pldj(x0, y0, x1, y1)"
  theta = a2(0, 2*pi, 18)
  x = cos(theta)
  y = sin(theta)
  pldj(x, y, transpose (x), transpose (y))
  pltitle("Seventeen Pointed Stars")
  limits(square = 1)
  print "All 17 pointed stars."
  if quitnow(): return
  limits(square = 0)
  fma()
  unzoom()
  ###
  print "Test 23: Lissajous animation"
  print "First run without animation mode"
  print "Second run with animation mode"
  print "Press RETURN to continue"
  if quitnow(): return
  lissajous(0)
  print "Press RETURN to continue"
  if quitnow(): return
  lissajous(1)
  if quitnow(): return
  ####
###

def quitnow(prompt=""):
  s = raw_input(prompt)
  if s == 'q': return 1
  return 0

def lissajous(animation):
#/* DOCUMENT lissajous
#     runs the Yorick equivalent of an old graphics performance test
#     used to compare PLAN, ALMA, and Basis with LTSS TMDS graphics.
#   SEE ALSO: testg, grtest
# */

#  /* Two figures with (x,y)= (cx,cy) + size*(cos(na*t), sin(nb*t+phase))
#     -- the centers describe semi-circular arcs of radius rc.  */
  import os
  t= 2*pi*arange(400, typecode = Float)/399.0
  na1= 1
  nb1= 5
  na2= 2
  nb2= 7
  rc1= 40.
  rc2= 160.
  size= 40.
  phase= theta= 0.

  n = 200;   #/* number of frames in animation */
  dtheta= pi/(n-1)
  dphase= 2*pi/(n-1)

  if animation == 1:
    animate(1)
  else:
    animate(0)

  fma()

  t0 = os.times()

  for i in range(n):
    cost= cos(theta)
    sint= sin(theta)
    x= rc1*cost+size*cos(na1*t)
    y= rc1*sint+size*sin(nb1*t+phase);
    plg( y, x)
    x= rc2*cost+size*cos(na2*t)
    y= rc2*sint+size*sin(nb2*t+phase);
    plg( y, x)
    fma()
    theta= theta + dtheta;
    phase= phase + dphase;
  
  t1 = os.times()
  print "Timing summary for Lissajous"
  print "Wall clock:", t1[4] - t0[4]
  print "User cpu:", t1[0] - t0[0]
  print "Sys cpu:", t1[1] - t0[1]

  if animation == 1:
    #/* turn off animation and pop up final frame again */
    animate(0)
    x= -rc1+size*cos(na1*t)
    y= size*sin(nb1*t);
    plg(y, x)
    x= -rc2+size*cos(na2*t)
    y= size*sin(nb2*t);
    plg(y, x)

# (nr X nc) array, each row an integer in the sequence [0, 1, ... ]
def a1(nr,nc):
  return reshape (array(concat1(nr,nc), Float), (nr,nc))

# (n-1) square array, with each row == spanz(lb,ub,n)
def a2(lb,ub,n):
  return reshape (array((n-1)*spanz(lb,ub,n), Float), (n-1,n-1))
  
# (n) square array, with each row == span(lb,ub,n)
def a3(lb,ub,n):
  return reshape (array(n*span(lb,ub,n), Float), (n,n))
  
# Half-hearted attempt at span(), which returns N equally-spaced
# values in the sequence [lb, ..., ub]
def span(lb,ub,n):
  if n < 2: raise ValueError, '3rd arg must be at least 2'
  b = lb
  a = (ub - lb)/(n - 1.0)
  #return a*arange(n) + b # if you want an array
  return map(lambda x,A=a,B=b: A*x + B, range(n)) # if you want a list

# Half-hearted attempt at span()(zcen), which returns N-1 "zone-centered"
# values in sequence (lb, ..., ub)
def spanz(lb,ub,n):
  if n < 3: raise ValueError, '3rd arg must be at least 3'
  c = 0.5*(ub - lb)/(n - 1.0)
  b = lb + c
  a = (ub - c - b)/(n - 2.0)
  return map(lambda x,A=a,B=b: A*x + B, range(n-1))

# magnitude or distance function
def mag(*args):
  r = 0
  for i in range(len(args)):
    r = r + args[i]*args[i]
  return sqrt(r)

def run():
  grtest()

def demo():
  grtest()
#################################
#if __name__ == '__main__': run()
