# $Id$
# Copyright (c) 1996, 1997, The Regents of the University of California.
# All rights reserved.  See Legal.htm for full text and disclaimer.

from gist import *
from shapetest import *
from yorick import *

#
#  MOVIE.PY
#  Support functions for making animated sequences.
#
#  $Id$
#

#    Copyright (c) 1997.  The Regents of the University of California.
#                  All rights reserved.

def movie (draw_frame, time_limit = 120., min_interframe = 0.0,
   bracket_time = array ([2., 2.], Float ), lims = None, timing = 0) :

   """
   movie (draw_frame, time_limit = 120., min_interframe = 0.0,
      bracket_time = array ([2., 2.], Float ), lims = None, timing = 0)
     runs a movie based on the given DRAW_FRAME function.  The movie
     stops after a total elapsed time of TIME_LIMIT seconds, which
     defaults to 60 (one minute), or when the DRAW_FRAME function
     returns zero.

     note: All but the first argument are keyword arguments, with
     defaults as shown.

     def draw_frame(1) :
       # Input argument i is the frame number.
       # draw_frame should return non-zero if there are more
       # frames in this movie.  A zero return will stop the
       # movie.
       # draw_frame must NOT include any fma command if the
       # making_movie variable is set (movie sets this variable
       # before calling draw_frame)

     If MIN_INTERFRAME is specified, a pause will be added as
     necessary to slow down the movie.  MIN_INTERFRAME is a time
     necessary to slow down the movie.  MIN_INTERFRAME is a time
     in seconds (default 0).

     The keyword bracket_time= (again a time in seconds) can be
     used to adjust the duration of the pauses after the first
     and last frames.  It may also be a two element array [beg, end].
     If the pause at the end is greater than five seconds, you will
     be prompted to explain that hitting <RETURN> will abort the final
     pause. Well, the Python version does not have this capability.

     timing = 1 enables a timing printout for your movie.

     If every frame of your movie has the same limits, use the
     limits command to fix the limits before you call movie.

   BUG:  If you hit <RETURN> to start a movie early, it will not
         pause at the end of the movie at all.  You probably should
         not use long initial pauses.
   """

   global movie_timing, making_movie
   if is_scalar (bracket_time) :
      bracket_time = array ( [bracket_time, bracket_time], Float )
   elif len (bracket_time) == 1 :
      bracket_time = array ( [bracket_time [0], bracket_time [0]], Float )

   elapsed = zeros (3)
   this_frame = zeros (3)

   window (wait = 1, style = "nobox.gs") # make sure window is ready to draw

   fma ( )
   animate (1)
   making_movie = 1

   i = 0
   if time == 1 :
      elapsed = timer_ (elapsed)
      elapsed0 = array(elapsed, copy = 1)
   i = i + 1
   more = draw_frame (i)
   if lims != None:
      limits (lims [0], lims [1], lims [2], lims [3])
   else:
      limits (square = 1)
   fma ( )
   if time == 1 :
      [elapsed, this_frame] = timer_ (elapsed, this_frame)
      wait = bracket_time [0] - this_frame [2]
      waited = waited0 = 0.0
      if wait > 0 :
         if wait > 5 :
            print "Movie starts in", wait, "secs, or when you hit <RETURN>"
         time.sleep (wait)
         waited0 = waited0 + wait
      else :
         wait = min_interframe - this_frame [2]
         if wait > 0 :
            time.sleep (wait)
            waited0 = waited0 + wait

   while (more) :
      this_frame = zeros (3)
      i = i + 1
      more = draw_frame (i)
      if lims != None:
         limits (lims [0], lims [1], lims [2], lims [3])
      else :
         limits (square = 1)
      fma ( )
      if not more :
         break
      if timing == 1 :
         [elapsed, this_frame] = timer_ (elapsed, this_frame)
         if (not more or (elapsed [2] - elapsed0 [2] > time_limit)) :
            break
         wait = min_interframe - this_frame [2]
         if wait > 0 :
            time.sleep (wait)
            waited0 = waited0 + wait

   if timing == 1 :
      wait = bracket_time [1] - this_frame [2]
      if wait > 0 :
         if wait > 5 :
            print "Holding last frame for", wait, "secs, or hit <RETURN>"
         time.sleep (5)      # wait)
         waited0 = waited0 + wait

      elapsed = timer_ (elapsed)

      e = elapsed - elapsed0
      movie_timing = [e [0], e [1], e [2], i, waited, waited0]

      movie_stats ()

   animate (0)

def movie_stats ( *timing ) :

   """
   movie_stats ( ) or movie_stats ( timing )
     prints statistics from the last movie command, or from the
     command which produced TIMING.  TIMING is the contents of the
     movie_timing external variable after the movie command completes.
   """

   if len (timing) > 0 :
      timing = timing [0]
   else :
      timing = movie_timing
   cpu = timing [0] + timing [1]
   print "user", timing [0], "system", timing [1]
   wall = timing [2]
   nframes = timing [3]
   waited = timing [4] + timing [5]
   wait = timing [5]

   print "  Wall(sec)  Wait(sec)  CPU(sec)"
   print "  %9.3f  %9.3f  %8.3f   %ld frames\n" % (wall, waited, cpu, nframes)
   print "  %9.3f  %9.3f  %8.3f   per frame\n" % \
      ( (wall - waited) / nframes, wait / ( (nframes > 1) * (nframes - 1) + \
      (nframes <= 1) * 1), cpu / nframes)
