# Gist axes.gs drawing style
# $Id$

# A single coordinate system on a portrait page with central axes
# Legends: two columns below viewport, contours in single column to right

# See work.gs for description of meanings

landscape= 0

default = {
  viewport= { 0.19, 0.60, 0.44, 0.85 },

  ticks= {

    horiz= {
      nMajor= 7.5,  nMinor= 50.0,  logAdjMajor= 1.2,  logAdjMinor= 1.2,
      nDigits= 3,  gridLevel= 1,  flags= 0x13c,
      tickOff= 0.0007,  labelOff= 0.009,
      tickLen= { 0.0143, 0.0091, 0.0052, 0.0026, 0.0013 },
      tickStyle= { color= -2,  type= 1,  width= 1.0 },
      gridStyle= { color= -2,  type= 1,  width= 1.0 },
      textStyle= { color= -2,  font= 0x08,  height= 0.0182,
        orient= 0,  alignH= 0,  alignV= 0,  opaque= 0 },
      xOver= 0.395,  yOver= 0.4024 },

    vert= {
      nMajor= 7.5,  nMinor= 50.0,  logAdjMajor= 1.2,  logAdjMinor= 1.2,
      nDigits= 4,  gridLevel= 1,  flags= 0x13c,
      tickOff= 0.0007,  labelOff= 0.009,
      tickLen= { 0.0143, 0.0091, 0.0052, 0.0026, 0.0013 },
      tickStyle= { color= -2,  type= 1,  width= 1.0 },
      gridStyle= { color= -2,  type= 1,  width= 1.0 },
      textStyle= { color= -2,  font= 0x08,  height= 0.0182,
        orient= 0,  alignH= 0,  alignV= 0,  opaque= 0 },
      xOver= 0.150,  yOver= 0.4024 },

    frame= 0,
    frameStyle= { color= -2,  type= 1,  width= 1.0 }}}

# The one coordinate system matches the default template exactly
system= { legend= "System 0" }

legends= {
  x= 0.04698,  y= 0.360,  dx= 0.3758,  dy= 0.0,
  textStyle= { color= -2,  font= 0x00,  height= 0.0156,
    orient= 0,  alignH= 1,  alignV= 1,  opaque= 0 },
  nchars= 36,  nlines= 20,  nwrap= 2 }

clegends= {
  x= 0.6182,  y= 0.8643,  dx= 0.0,  dy= 0.0,
  textStyle= { color= -2,  font= 0x00,  height= 0.0156,
    orient= 0,  alignH= 1,  alignV= 1,  opaque= 0 },
  nchars= 14,  nlines= 28,  nwrap= 1 }
