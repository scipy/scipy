Notes on test2d
---------------

You can set the macros at the top to 0 to omit any calls to the
corresponding feature.  This is supposed to allow you to begin
debugging before some functions are implemented -- maybe it would be
better to assume that there are at least stubs for every play
function.

------------------------------------------------------------------------

If TRY_PSTDINIT is 1, then the stdio stream/window should become
functional.  There will be a prompt (> or idle>, the latter if
TRY_ALARMS is 1) and test2d should respond to commands typed at the
prompt.  The commonds vary depending on which TRY_* switches are set;
here are all of the possibilities:

help
  prints all available commands
quit [value]
  exits the program with return value (default 0)
idle [ncalls]
  prints "on_idle called..." ncalls times (default 1)
set alarm_secs
  sets an alarm to ring timeout_secs from now
  - a message should appear on p_stdout when the alarm rings
clr [alarm_number]
  clears the n-th (default earliest) alarm
raise [signal_number]
  raises signal_number (see /usr/include/signal.h)
  - if signal number omitted, divides by 0 to generate SIGFPE
  - sending SIGINT to test2d (by ^c) should also by fielded

ls path
cd path
pwd
mv old_path new_path
rm path
mkdir path
rmdir path
  imitate 1-path versions of UNIX directory management utilities
  - these must all accept UNIX style / directory separators and
    .. and . relative pathname syntax
  - ls must correctly distinguish directories from ordinary files
    and give correct size in bytes (when the file is readable and
    its size is meaningful)

feep
  rings the bell (for the TRY_GRAPHICS window)
resize width height
  resizes the TRY_GRAPHICS window to widthXheight pixels
    (silently ignored for sizes <50 or >1000 pixels)
redraw
  redraws the test2d window and the test2d advanced window
  - should happen automatically on iconify/deiconify or exposure
private
  destroy test2d advanced window and recreate it with private
  colormap -- subsequent calls alternate private/shared
rgb
  destroy test2d advanced window and recreate it with 5-9-5 rgb
  colormap -- subsequent calls alternate 5-9-5/indexed
  - this is a noop except on pseudocolor (8-bit) displays
palette
  cycle from grayscale (1) to rainbow (2) to no palette (0) in
  test2d advanced window
  - no palette does not do a redraw, so the palette swatch may
  look strange
rgbcell
  draws a fading palette swatch under the indexed palette swatch
  using p_rgb_cells().  On pseudocolor (8-bit) displays, this will
  force the 5-9-5 rgb colormap to be installed

------------------------------------------------------------------------

If TRY_GRAPHICS is 1, then a window pops up on the default screen when
test2d starts.  The window title is "test2d window".  If all other
features are on, it looks like this:

[]                                          []
  Corner boxes show four sides?         FFFF
  String aligned within boxes?          FFFF
  Fills aligned (UR corner)?
  Dots aligned (.+.)?                   DDD
  Segments aligned (-+-)?               DDD
  Circs/rects OK (1 pixel gaps)?
                                        SSS

    _  	     _        _        ____     CCC
   /  |	    /  |     /  |      |X |
   | _/	    | _/     | _/      | X|     RRR
                               ----
  -------------=================XXXXXXXXXXX
    _  	     _        _        ___  
   /  |	    /  |     /  |      | / /|
   | _/	    | _/     | _/      |/ /_|


  on_focus  on_key  on_click  on_motion
[]                                          []

You need to examine this picture carefully for the following features
(some of which are suggested by the text):

1. The four boxes in the window corners should show all four sides,
each one pixel thick, with no gaps between the edge of the window and
the edges of the box.

2. "String aligned within boxes?" should be surrounded by a bounding
box that just touches the extreme characters on all sides.  The words
"aligned within" should be additionally enclosed within a box whose
bottom is the baseline of the text (not including descenders like the
outer box).

3. The FFF object in the upper right corner tests alignment of p_fill.
The pattern should be symmetric about both its horizontal and vertical
centerlines, and also symmetric about both 45 degree lines through the
center.  (Fills aligned?)

As best I can, here is a verbal description of what the object is
supposed to look like: Begin with a large black (foreground) square.
Divide each side in fourths, and remove the quarter size triangle at
each corner, leaving an octagon.  Next, remove the central four of the
sixteen smaller squares, leaving an octagon with a square hole.
Replace the four central triangles, to get an octagon with a square
hole that has a diamond-shaped plug.  Now remove a one pixel
horizontal and vertical centerline, dividing the pattern into four
pieces.  Finally, add one pixel dark (foreground) crosshairs aligned
with the one pixel gaps protruding beyond the edges of the octagon.
(Under X11, the downward facing back triangles have one extra row of
black pixels.  This is not easily visible, but it does break the
precise up-down symmetry of the figure.)

4. The DDD object (below FFF) is a crosshair with a one pixel gap,
then a dot at the tip of each of the four hair ends.  The dots should
be exactly one pixel, the gaps should be one pixel, and the dots
should be colinear with the crosshairs.  (Dots aligned?)

5. The SSS object (below DDD) is a crosshair with one pixel gaps in
each of the four lines.  (Segments aligned?)

6. The CCC object (below SSS) is two one pixel square boxes.  The left
box contains a filled circle; the right box contains a one pixel wide
empty circle.  There should be a one pixel gap between box and circle
on all four sides.  (Circs OK?)

6.5 The RRR object (below CCC) is a rectangular box.  The outer edge
of the box is 1 pixel wide, then a 1 pixel gap, then a 4 pixel wide
box, then a 1 pixel gap, then a 3 pixel wide box, a 1 pixel gap, a 1
pixel wide box, a 1 pixel gap, and a solid rectangular center.  The
corners of the thick boxes should be sharp.  (rects OK?)

7. The upper row of parenthesis shapes are one pixel wide curves
testing the line styles.  From left to right, the six styles are:
solid, dash, dot, dash-dot, dash-dot-dot, and solid.

8. To the right of 7 is the clip test.  This should look like a 2x2
checkerboard with 1 pixel border line that leaves a 1 pixel gap around
the dark squares of the checkeboard.

9. The next test is line width; the line increases in one pixel width
increments from one pixel wide to nine pixels wide.  The line should
remain centered, so that the nine pixel line has four pixels above and
four pixels below the row containing the original one pixel line at
the left.

10. The lower row of parenthesis shapes repeats the line style test
for five pixel wide lines.  The dash patterns should scale to keep the
patterns good looking.

11. To the right of the wide line style test are two triangles which
test the join style (P_SQUARE).  The upper left triangle should have
sharp corners; the lower right triangle should have rounded corners.
(The default for the previous thick lines in 9 and 10 is rounded
ends.)

12. on_focus: When you give the window keyboard focus, the word "IN"
appears above on_focus; when the window loses focus, the word "OUT"
appears.

13. on_key: When you type characters in the window, they are echoed
above on_key.  The key is indicated to the right; control characters
will have a ^ before the letter (^@ is \0, Ctl-SPACE on most
keyboards, and ^? is delete).  Function keys are prefixed by the
letter F, so that F1 thru F9 look like themselves, F10, F11,
etc. become Fa, Fb, etc.  Other function keys are: FH - home, FE -
end, FI - insert, FU - page up, FD - page down, F> - right arrow, F< -
left arrow, F^ - up arrow, FV - down arrow.  To the left of the key
name is a row of symbols indicating the shift keys being held when the
key was pressed.  If multiple shift keys are in effect, all are
indicated.  They are: S - shift, C - control, M - meta, A - alt.
Additionally, a # appears if the key is on the numeric keypad.
Finally, the key symbol is replaced by "->" if a mouse button is
clicked, in order to show the shift keys in effect at click time.

14. on_click: Indicates two masks: the mask for the buttons held down
at the time the current event arrived is on the left; the mask for the
button which caused the current event is on the right.  Left, center,
and right button masks are 010, 020, 040 respectively.  When you press
and hold the left button, you should see "000 010".  When you release
the left button, it should change to "010 010".  Chording need not
work perfectly (it doesn't on my X11 server).  However, the shift keys
(indicated over on_key) should indicate properly.

Double clicking on the word "on_click" should highlight "on_click" and
make it the primary selection.  Making a selection in another
application window should dehighlight "on_click".  When it deselects,
the new selection value should be printed using p_stdout to the
test2d stdout stream/window.  A second double click on on_click will
also dehighlight and deselect it.

15. on_motion: Moving the mouse around in the window indicates the
current x,y coordinates in pixels (0,0 is at upper left corner).  When
you hold down a mouse button, the odometer tracking should continue
outside the window (this is also permissible when no button is pressed
but the window has keyboard focus).

16. If you use the window manager to resize the test2d window, the
small squares should remain in the new corners when the window is
redrawn.

17. You should be able to use the window manager to destroy the
test2d window.  This will print a message on p_stdout.

18. If TRY_OFFSCREEN is 1, then a one pixel thick square appears below
the upper left corner square.  If TRY_ALARMS and TRY_GUI are also 1,
then clicking at the left edge of the window causes the square to
bounce up and down between the two left edge corner squares.  A second
click stops the animation.

19. If TRY_CURSOR is 1, then the area below the on_focus line is
divided into 13 regions from left to right; the cursor should change
shape as the most moves from left to right across the bottom of the
window.  In order, the shapes are: select arrow, crosshair, text tool,
up, down, right, left, up-down, and right-left scrollbar arrows,
four-way arrow, rotation indicator, death head, sliding hand, and
cursor off.

20. If TRY_MENU is 1, then double clicking in the upper part of the
window will pop up a menu below and right of the cursor.  If you
release the mouse button before moving into the menu, it should stay
popped up, popping down when you click again.  If you move the mouse
into the menu before releasing the button, the menu will pop down on
release.  When you move the mouse into the center menu entry, a
submenu should pop up to the right of the menu; releasing the mouse
button or clicking on the center menu entry of the primary menu should
leave the submenu up.  Click (or release) on a submenu entry to pop
both down.  Moving the mouse onto the primary menu to an entry other
than the center one should pop down the submenu but leave the primary
menu up.  The message printed on p_stdout should agree with the menu
slection ("item 0" means nothing was selected -- nothing was
highlighted when you clicked or released the mouse button.)

------------------------------------------------------------------------

A second window titled "test2d advanced" should also pop up.  This
one contains:

1. A more comprehensive font test, which looks like this:

    Courier 14 Bold Italic Bold-Italic
    Helvetica 10 Bold Italic Bold-Italic
    Times 12 Bold Italic Bold-Italic
    Newcentury 18 Bold Italic Bold-Italic
    Symbol 14

all of these are supposed to appear in the typeface mentioned if that
font is available.

2. A square consisting of the string "Square" repeated in all four
text orientations, with a one pixel bounding box drawn around each.
All four should be pixel-for-pixel identical except for orientation.

3. A white background containing a sample of each of the standard
colors.  The P_BG color is obtained as a P_XOR of the P_FG sample.

4. A set of six colored checkerboards of various sizes and shapes
testing the p_ndx_cells() function.  These should all be 5x5
checkerboards with the same set of rows and columns in the same order:

      bg   fg  black  white   red
    green blue cyan  magenta yellow
      bg   fg  black  white   red
    green blue cyan  magenta yellow
      bg   fg  black  white   red

Each cell array has a one pixel fg border drawn around it.

5. The upper left 20x20 pixels of the text orientation test (2.) is
reproduced to the right of the checkerboards.  This was read from the
window using p_rgb_read() then drawn with p_ndx_cells().

6. A palette testing swatch.  This shows the current palette, which
can be cycled from grayscale (1) to rainbow (2) to no palette (0) by
means of the p_stdin palette command.  During all these palette
changes, the standard colors should remain invariant, and minimal
screen flashing should occur, even when a private colormap is being
used.  With shared colors, you will want to run other clients which
steal large numbers of colors to observe the interactions.

7. The p_rgb_cells() test swatch, which is turned on only after the
rgbcell command is issued.
