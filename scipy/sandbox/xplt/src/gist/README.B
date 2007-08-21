
      The BENCH program benchmarks and tests Gist

      $Id$

The execute line for the program is

     bench [useps [host:server.screen]]

where useps is 0 (the default) to use a CGM file as the hardcopy
device, and 1 to use a PostScript file.  If the optional X display
argument is not given, the Gist default is to use the DISPLAY
environment variable, like other X programs.  (If you built bench
with NO_XLIB, both arguments are ignored-- the "display" device is
a CGM file, and the hardcopy device is a PostScript file.)  The
output files, if any, will be "bench.ps" for PostScript output,
and "bench.cgm" for CGM output.  (Gist automatically families large
CGM files to about 1 Meg per file, so check for "benci.cgm",
"bencj.cgm", etc.  The gist CGM browser automatically jumps
between files in a family.)

The program waits for you to type a command, then exercises the
specified feature of Gist.  Some of the commands merely set attributes
rather than producing any graphics output.  Some commands add a
graphical element to the current drawing; others produce movies
to get an idea of Gist performance during animation.  Most commands
test the `D' and `H' level routines, but a few `A' and `P' level
benchmarks allow you to see how much performance penalty is incurred
by the higher level routines.  These low level tests also provide
fair comparisons of Gist with other graphics packages.

The 'help' command produces the following output, describing the
various commands:

	miggle[50] bench
	Initializing 51x51 contour function, STANDBY...  DONE
	bench> help

	   This is the GIST benchmark/test program, commands are:

	    Movies to test low level Gist functions--
	0   - raw performance test GaMesh (direct)
	1   - raw performance test GaMesh (animated)
	2   - raw performance test GaFillMesh (direct)
	3   - raw performance test GaFillMesh (animated)
	clr - toggle interframe clear for tests 0-3
	c   - clear screen
	im  - raw performance test GpCells
	clip  - toggle floating point clip
	zoom  - toggle zoom for 0-3, im, pli

	    Property toggles for high level tests--
	bnd  - mesh boundary          odd0  - mesh region to 0
	odd2  - mesh region to 2      reg0  - plot region to 0
	reg1  - plot region to 1      reg2  - plot region to 2
	nocopy - toggle noCopy mode for mesh-based tests
	cmode - toggle dumping of colormap in hcp file
	earth - use earthtones palette (default).  Others are:
	     stern   rainbow   heat   gray   yarg
	work  - use work stylesheet (default).  Others are:
	     boxed   axes
	wide      - toggle wide lines in all line drawings
	linetype  - cycle through 5 line types in all line drawings
	closed    - toggle closed lines in plg tests
	rays      - toggle ray arrows in plg tests
	smooth    - cycle through 3 hcp smoothnesses in plg tests
	limits, xmin, xmax, ymin, ymax, nice=1/0, square=1/0,
		restrict=1/0   - set limits (default limits,e,e,e,e)
	logxy, 1/0, 1/0   - set/reset log scaling for x,y axes
	animate  - toggle animation mode
	marks    - toggle occasional curve markers for line drawings
	legends  - toggle dumping of legends into hcp file

	    Movies to test high level Gist functions--
	plg  - test GdLines
	plm  - test GdMesh
	plf  - test GdFillMesh
	plv  - test GdVectors
	plc  - test GdContours

	    Single frame tests of high level Gist functions--
	plg1  - test GdLines
	plm1  - test GdMesh
	plf1  - test GdFillMesh
	plc1  - test GdContours
	pli   - test GdCells
	plt   - test GdText
	pls   - test GdLines in polymarker mode
	pldj  - test GdDisjoint
	txin  - test GdText     txout  - 2nd test of GdText
	calib - text calibration frame
	cfont - toggle font for calibration frame

	    Tests of high level Gist control functions--
	fma    - frame advance
	hcp    - send current frame to hardcopy file
	hcpon  - send every frame to hardcopy file
	hcpoff - stop sending every frame to hardcopy file
	redraw - redraw the X window

	quit  - exit bench program (just q works too)

	bench>


Here are 3 sample sessions with the bench program.  Input lines that you
type are marked with bench>.  Most of the features of Gist are
demonstrated here.  The three sessions are:
 (1) A performance test, both for "fastest possible" rendering using
     A-level and P-level Gist routines, and for the D-level display
     list routines that could reasonably be used by an interpreter.
     Gist's "animation" mode is exhibited (smooth animation requires
     a picture to be built up offscreen, then "blitted" onscreen).
     The polyline and filled area primitives are exercised.
 (2) Generate a relatively small CGM file which contains nearly
     everything that Gist will ever throw at a CGM interpreter.
     This file is useful for exercising the gist CGM browser, as well
     as for seeing how badly Gist CGM file break your second favorite
     CGM interpreter.
 (3) Exercise other high level features of Gist, especially the
     tick drawing and labeling routine.  Also exhibits D-level style
     and palette files.

For the record, "miggle" is a Sun-4/330 running the X11R4 sample server
from MIT.  Gist was compiled by gcc with no optimization (optimization
speeds up the low level benchmarks by about 30%).  The display is to
an HP 730 workstation, which has a very fast X server.


(1)---------------------------------------------------------------------


IMPORTANT-- If you care about the timings in the first test, be sure
to clear all other jobs off of the machine which is running your
X server, AND DON'T MOVE THE MOUSE OR TYPE DURING A TEST!

	miggle[51] bench
	Initializing 51x51 contour function, STANDBY...  DONE
	bench> 0
	   Begin benchmark 0 -- GaMesh (direct).
	elapsed seconds: 3.370000 user, 0.710000 system, 13.409998 wall
	Plots per wall second= 3.728561

	bench> 1
	   Begin benchmark 1 -- GaMesh (animated).
	elapsed seconds: 3.210000 user, 0.820000 system, 15.399998 wall
	Plots per wall second= 3.246754

	bench> 2
	   Begin benchmark 2 -- GaFillMesh (direct).
	elapsed seconds: 9.770000 user, 2.260000 system, 36.759997 wall
	Plots per wall second= 1.360174

	bench> 3
	   Begin benchmark 3 -- GaFillMesh (animated).
	elapsed seconds: 9.790000 user, 2.380000 system, 22.969997 wall
	Plots per wall second= 2.176753

	bench> clip
	   Toggle floating point clip (now on)

	bench> 0
	   Begin benchmark 0 -- GaMesh (direct).
	elapsed seconds: 4.870000 user, 0.780000 system, 15.829997 wall
	Plots per wall second= 3.158560

	bench> 1
	   Begin benchmark 1 -- GaMesh (animated).
	elapsed seconds: 5.010000 user, 0.630000 system, 15.849998 wall
	Plots per wall second= 3.154575

	bench> 2
	   Begin benchmark 2 -- GaFillMesh (direct).
	elapsed seconds: 13.390000 user, 2.410000 system, 40.499997 wall
	Plots per wall second= 1.234568

	bench> 3
	   Begin benchmark 3 -- GaFillMesh (animated).
	elapsed seconds: 13.570000 user, 2.330000 system, 28.999997 wall
	Plots per wall second= 1.724138

	bench> fma
	   Frame advance

	bench> plm
	   Begin mesh benchmark.
	elapsed seconds: 6.940000 user, 0.820000 system, 15.439997 wall
	Plots per wall second= 3.238343

	bench> animate
	   Toggle animation mode (now on)

	bench> plm
	   Begin mesh benchmark.
	elapsed seconds: 7.000000 user, 0.710000 system, 17.179998 wall
	Plots per wall second= 2.910361

	bench> animate
	   Toggle animation mode (now off)

	bench> plf
	   Begin filled mesh benchmark.
	elapsed seconds: 19.500000 user, 2.180000 system, 37.429997 wall
	Plots per wall second= 1.335827

	bench> animate
	   Toggle animation mode (now on)

	bench> plf
	   Begin filled mesh benchmark.
	elapsed seconds: 19.230000 user, 2.320000 system, 37.099997 wall
	Plots per wall second= 1.347709

	bench> quit
	120.5u 18.8s 7:42 30% 0+568k 1+1io 1pf+0w

Several features of these benchmark results are worth contemplating:
Most important, for the X Engine (at least on my Sun-4/330), the time
is dominated by the rendering time of the server, NOT by the CPU time
required by Gist.  Secondly, of the CPU time used by Gist, the time to
do the floating point clip explains much of the additional time used
by the Gist D-level routines (there is no very practical way to avoid
the clipping operation in the D-level routines).  Thirdly, animation
mode costs little in performance.


(2)---------------------------------------------------------------------

Generate a CGM test file for exercising the gist browser:

	miggle[52] bench
	Initializing 51x51 contour function, STANDBY...  DONE
	bench> plg1
	   Plot small graph.

	bench> marks
	   Toggle markers (now on)

	bench> plg1
	   Plot small graph.

	bench> marks
	   Toggle markers (now off)

	bench> rays
	   Toggle rays (now on)

	bench> plg1
	   Plot small graph.

	bench> rays
	   Toggle rays (now off)

	bench> wide
	   Toggle wide lines (now on)

	bench> plg1
	   Plot small graph.

	bench> linetype
	   Toggle line type (now dashed)

	bench> plg1
	   Plot small graph.

	bench> linetype
	   Toggle line type (now dotted)

	bench> plg1
	   Plot small graph.

	bench> linetype
	   Toggle line type (now dash-dot)

	bench> plg1
	   Plot small graph.

	bench> linetype
	   Toggle line type (now dash-dot-dot)

	bench> plg1
	   Plot small graph.

	bench> linetype
	   Toggle line type (now solid)

	bench> wide
	   Toggle wide lines (now off)

	bench> plg1
	   Plot small graph.

	bench> linetype
	   Toggle line type (now dashed)

	bench> plg1
	   Plot small graph.

	bench> linetype
	   Toggle line type (now dotted)

	bench> plg1
	   Plot small graph.

	bench> linetype
	   Toggle line type (now dash-dot)

	bench> plg1
	   Plot small graph.

	bench> linetype
	   Toggle line type (now dash-dot-dot)

	bench> plg1
	   Plot small graph.

	bench> linetype
	   Toggle line type (now solid)

	bench> txout
	   Exterior text

	bench> hcp
	   Sent frame to hcp

	bench> fma
	   Frame advance

	bench> plt
	   Text test

	bench> hcp
	   Sent frame to hcp

	bench> fma
	   Frame advance

	bench> pli
	   Plot image.

	bench> hcp
	   Sent frame to hcp

	bench> fma
	   Frame advance

	bench> plf1
	   GdFillMesh test

	bench> hcp
	   Sent frame to hcp

	bench> fma
	   Frame advance

	bench> cmode
	   Toggle hcp color mode (now dump)

	bench> plg1
	   Plot small graph.

	bench> pls
	   Plot small scatter.

	bench> pls
	   Plot small scatter.

	bench> pls
	   Plot small scatter.

	bench> pls
	   Plot small scatter.

	bench> pls
	   Plot small scatter.

	bench> pls
	   Plot small scatter.

	bench> pls
	   Plot small scatter.

	bench> plg1
	   Plot small graph.

	bench> txout
	   Exterior text

	bench> hcp
	   Sent frame to hcp

	bench> fma
	   Frame advance

	bench> pli
	   Plot image.

	bench> hcp
	   Sent frame to hcp

	bench> stern
	OK

	bench> fma
	   Frame advance

	bench> plf1
	   GdFillMesh test

	bench> hcp
	   Sent frame to hcp

	bench> quit
	miggle[53] mv bench.cgm check.cgm



Repeat this test with the command line "bench 1" to generate a PostScript
file which you can send to your printer to be sure that Gist doesn't
produce anything which confuses it.  Alternatively, you can build the
gist previewer and run check.cgm through it to convert to PostScript:

	miggle[54] gist/gist check.cgm -b -ps check.ps

Or view the CGM file on your X server:

	miggle[55] gist/gist check.cgm
	check.cgm metafile description:
	Fri Oct 30 13:53:07 1992;  For: munro
	gist> f
	gist> f
	gist> f
	gist> f
	gist> f
	gist> f
	gist> draw 1
	gist> quit

You can also compare the relative sizes of PostScript and CGM output.
On a more extensive test, I found that:

  (1) A Gist PostScript file is roughly twice the size of the Gist CGM
      file produced by the same commands (two hex characters represent
      one byte).
  (2) The CPU time required by Gist to produce a PostScript file is two
      to three times the CPU time required by Gist to produce a CGM file
      from the same commands.
  (3) The compress utility will squeeze the PostScript file by about 60%,
      to a size just slightly smaller than the CGM file.  The CGM file
      compresses only 10-15% (not worth doing, but slightly smaller than
      the compressed PostScript).

The moral is- use Gist CGM output and don't bother to compress it.  The
gist CGM browser will translate Gist CGM files into exactly the same
PostScript file Gist would have produced in the first place.


(3)---------------------------------------------------------------------

The plg, plm, plc, plv, plf, plt, and pli commands to bench test
features of the proposed interactive interface routines of the
same names described in README.  The fma, redraw, limits, and
logxy commands are also similar to the proposed interface, except
that the limits command in bench takes the logarithms of the
limits, while the interactive routine would work in the intuitive
way.

Whenever the screen looks interesting, type "hcp" to send it to the
CGM file "bench.ps", which you can send to your printer.  (Note-
this will fail silently for the frame displayed at the end of any
of the movies, or for the "stale" frame displayed after fma.)


	miggle[67] bench
	Initializing 51x51 contour function, STANDBY...  DONE
	bench> plm1
	   GdMesh test

	bench> txin
	   Interior text

	bench> txout
	   Exterior text

	bench> fma
	   Frame advance

	bench> plt
	   Text test

	bench> plm1
	   GdMesh test

	bench> redraw
	   Redraw

	bench> fma
	   Frame advance

	bench> plm1
	   GdMesh test

	bench> limits, .02, 1.48, .31, .59
	   Setting limits

	bench> limits, -.2, 1.2, e, e
	   Setting limits

	bench> limits, -.2, e, e, e, square=1
	   Setting limits

	bench> limits, -.2, 1.5, e, e
	   Setting limits

	bench> limits, square=0
	   Setting limits

	bench> logxy, 1, 1
	   Setting logxy

	bench> limits, 105, 190, -.001, 1
	   Setting limits

	bench> limits, 17, 39.5, .001, .1
	   Setting limits

	bench> limits, 17, 32, 1, 100
	   Setting limits

	bench> limits, 16, 19, 1.e4, 1.e6
	   Setting limits

	bench> logxy, 0, 0
	   Setting logxy

	bench> limits, e, e, e, e
	   Setting limits

	bench> plg
	   Lissajous test
	elapsed seconds: 2.410000 user, 0.180000 system, 2.839997 wall
	Plots per wall second= 17.605653

	bench> animate
	   Toggle animation mode (now on)

	bench> plg
	   Lissajous test
	elapsed seconds: 2.320000 user, 0.150000 system, 4.079997 wall
	Plots per wall second= 12.254911

	bench> marks
	   Toggle markers (now on)

	bench> plg
	   Lissajous test
	elapsed seconds: 3.170000 user, 0.280000 system, 4.389997 wall
	Plots per wall second= 11.389529

	bench> wide
	   Toggle wide lines (now on)

	bench> limits, square=1
	   Setting limits

	bench> marks
	   Toggle markers (now off)

	bench> plg
	   Lissajous test
	elapsed seconds: 2.260000 user, 0.160000 system, 5.779997 wall
	Plots per wall second= 8.650523

	bench> wide
	   Toggle wide lines (now off)

	bench> rays
	   Toggle rays (now on)

	bench> plg
	   Lissajous test
	elapsed seconds: 3.090000 user, 0.270000 system, 3.409997 wall
	Plots per wall second= 14.662770

	bench> rays
	   Toggle rays (now off)

	bench> limits, square=0
	   Setting limits

	bench> plm
	   Begin mesh benchmark.
	elapsed seconds: 6.930000 user, 0.830000 system, 16.279997 wall
	Plots per wall second= 3.071254

	bench> marks
	   Toggle markers (now on)

	bench> plc
	   Begin contour benchmark.
	elapsed seconds: 15.810000 user, 0.310000 system, 16.179997 wall
	Plots per wall second= 3.090235

	bench> limits, .25, .75, .25, .75
	   Setting limits

	bench> plv
	   Begin vector benchmark.
	elapsed seconds: 12.650000 user, 0.860000 system, 14.929997 wall
	Plots per wall second= 3.348963

	bench> animate
	   Toggle animation mode (now off)

	bench> limits, e, e, e, e
	   Setting limits

	bench> pli
	   Plot image.

	bench> rainbow
	OK

	bench> redraw
	   Redraw

	bench> cmode
	   Toggle hcp color mode (now dump)

	bench> stern
	OK

	bench> redraw
	   Redraw

	bench> heat
	OK

	bench> gray
	OK

	bench> yarg
	OK

	bench> earth
	OK

	bench> fma
	   Frame advance

	bench> plg1
	   Plot small graph.

	bench> pls
	   Plot small scatter.

	bench> pls
	   Plot small scatter.

	bench> pls
	   Plot small scatter.

	bench> pls
	   Plot small scatter.

	bench> pls
	   Plot small scatter.

	bench> pls
	   Plot small scatter.

	bench> pls
	   Plot small scatter.

	bench> plg1
	   Plot small graph.

	bench> boxed
	OK

	bench> plg1
	   Plot small graph.

	bench> axes
	OK

	bench> plg1
	   Plot small graph.

	bench> work
	OK

	bench> plg1
	   Plot small graph.

	bench> quit
	55.9u 4.1s 9:32 10% 0+736k 6+1io 6pf+0w

Note that in non-dump cmode, a redraw is required after a palette change
to change the pixel numbers of the image.
