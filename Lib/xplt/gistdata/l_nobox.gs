# Gist nobox.gs drawing style
# $Id$

# A single coordinate system on a landscape page with no ticks or labels
# Legends: none

# See work.gs for description of meanings

landscape= 1

default = {
  legend= 0,

  viewport= { 0.244, 0.789, 0.127, 0.672 },

  ticks= {

    horiz= {
      nMajor= 7.5,  nMinor= 50.0,  logAdjMajor= 1.2,  logAdjMinor= 1.2,
      nDigits= 3,  gridLevel= 1,  flags= 0x000,
      tickOff= 0.0007,  labelOff= 0.0182,
      tickLen= { 0.0143, 0.0091, 0.0052, 0.0026, 0.0013 },
      tickStyle= { color= -2,  type= 1,  width= 1.0 },
      gridStyle= { color= -2,  type= 3,  width= 1.0 },
      textStyle= { color= -2,  font= 0x08,  height= 0.0182,
        orient= 0,  alignH= 0,  alignV= 0,  opaque= 0 },
      xOver= 0.395,  yOver= 0.350 },

    vert= {
      nMajor= 7.5,  nMinor= 50.0,  logAdjMajor= 1.2,  logAdjMinor= 1.2,
      nDigits= 4,  gridLevel= 1,  flags= 0x000,
      tickOff= 0.0007,  labelOff= 0.0182,
      tickLen= { 0.0143, 0.0091, 0.0052, 0.0026, 0.0013 },
      tickStyle= { color= -2,  type= 1,  width= 1.0 },
      gridStyle= { color= -2,  type= 3,  width= 1.0 },
      textStyle= { color= -2,  font= 0x08,  height= 0.0182,
        orient= 0,  alignH= 0,  alignV= 0,  opaque= 0 },
      xOver= 0.100,  yOver= 0.350 },

    frame= 0,
    frameStyle= { color= -2,  type= 1,  width= 1.0 }}}

# The one coordinate system matches the default template exactly
system= { legend= "System 0" }

legends= { nlines= 0 }

clegends= { nlines= 0 }
