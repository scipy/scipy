# Gist vg.gs drawing style
# $Id$

# A single coordinate system on a portrait page with big labels
# This is more suitable for scale reduction or projection than work.gs.
# Legends: None

# See work.gs for description of meanings

landscape= 0

default = {
  legend= 0,

  viewport= { 0.22, 0.62, 0.45, 0.85 },

  ticks= {

    horiz= {
      nMajor= 4.5,  nMinor= 40.0,  logAdjMajor= 1.2,  logAdjMinor= 1.2,
      nDigits= 3,  gridLevel= 1,  flags= 0x033,
      tickOff= 0.0007,  labelOff= 0.024,
      tickLen= { 0.0143, 0.0091, 0.0052, 0.0026, 0.0013 },
      tickStyle= { color= -2,  type= 1,  width= 4.0 },
      gridStyle= { color= -2,  type= 3,  width= 2.0 },
      textStyle= { color= -2,  font= 0x08,  height= 0.0312,
        orient= 0,  alignH= 0,  alignV= 0,  opaque= 0 },
      xOver= 0.405,  yOver= 0.360 },

    vert= {
      nMajor= 4.5,  nMinor= 40.0,  logAdjMajor= 1.2,  logAdjMinor= 1.2,
      nDigits= 4,  gridLevel= 1,  flags= 0x033,
      tickOff= 0.0007,  labelOff= 0.024,
      tickLen= { 0.0143, 0.0091, 0.0052, 0.0026, 0.0013 },
      tickStyle= { color= -2,  type= 1,  width= 4.0 },
      gridStyle= { color= -2,  type= 3,  width= 2.0 },
      textStyle= { color= -2,  font= 0x08,  height= 0.0312,
        orient= 0,  alignH= 0,  alignV= 0,  opaque= 0 },
      xOver= 0.160,  yOver= 0.360 },

    frame= 0,
    frameStyle= { color= -2,  type= 1,  width= 1.0 }}}

# The one coordinate system matches the default template exactly
system= { legend= "System 0" }

legends= { nlines= 0 }

clegends= { nlines= 0 }
