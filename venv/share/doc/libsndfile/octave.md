---
layout: page
---

# libsndfile and GNU Octave

[GNU Octave](http://www.octave.org/) is a high-level interactive language for
numerical computations. There are currently two development streams, a stable
2.0.X series and a development 2.1.X series. Octave reads and writes data in
binary formats that were originally developed for
[MATLAB](http://www.mathworks.com/). Version 2.0.X of Octave uses binary data
files compatible with MATLAB version 4.2 while Octave 2.1.X uses binary data
files compatible with MATLAB version 5.0 as well as being able to read the older
MATLAB 4.2 format.

From version 1.0.1 of libsndfile onwards, libsndfile has the ability of reading
and writing a small subset of the binary data files used by both versions of GNU
Octave. This gives people using GNU Octave for audio based work an easy method
of moving audio data between GNU Octave and other programs which use libsndfile.

For instance it is now possible to do the following:

* Load a WAV file into a sound file editor such as
  [Sweep](http://www.metadecks.org/software/sweep/).
* Save it as a MAT4 file.
* Load the data into Octave for manipulation.
* Save the modified data.
* Reload it in Sweep.

Another example would be using the MAT4 or MAT5 file formats as a format which
can be easily loaded into Octave for viewing/analyzing as well as a format which
can be played with command line players such as the one included with
libsndfile.

## Details

Octave, like most programming languages, uses variables to store data, and
Octave variables can contain both arrays and matrices. It is also able to store
one or more of these variables in a file. When reading Octave files, libsndfile
expects a file to contain two variables and their associated data. The first
variable should contain a variable holding the file sample rate while the second
variable contains the audio data.

For example, to generate a sine wave and store it as a binary file which is
compatible with libsndfile, do the following:

    octave:1 > samplerate = 44100 ;
    octave:2 > wavedata = sin ((0:1023)*2*pi/1024) ;
    octave:3 > save sine.mat samplerate wavedata

The process of reading and writing files compatible with libsndfile can be made
easier by use of two Octave script files:

    octave:4 > [data fs] = sndfile_load ("sine.mat") ;
    octave:5 > sndfile_save ("sine2.mat", data, fs) ;

In addition, libsndfile contains a command line program which which is able to
play the correct types of Octave files. Using this command line player
**sndfile-play** and a third Octave script file allows Octave data to be played
from within Octave on any of the platforms which **sndfile-play** supports (at
the moment: Linux, MacOS X, Solaris and Win32).

    octave:6 > sndfile_play (data, fs) ;

These three Octave scripts are installed automatically in Octave's site script
directory when libsndfile is installed (except on Win32) ie when libsndfile is
being installed into /usr/local, the Octave scripts will be installed in
/usr/local/share/octave/site/m/.

There are some other Octave scripts for audio to be found
[here](http://octave.sourceforge.net/audio/index.html).
