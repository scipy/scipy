Name

      qhull, rbox         2003.1           2003/12/30
  
Convex hull, Delaunay triangulation, Voronoi diagrams, Halfspace intersection
 
      Documentation:
        html/index.htm

      Available from:
        <http://www.qhull.org>
	<http://savannah.nongnu.org/projects/qhull>
 
     Version 1 (simplicial only):
        <http://www.qhull.org/download/qhull-1.0.tar.gz>
        <http://www.qhull.org/download/qhull.sit.hqx>
       
      News and a paper:
        <http://www.qhull.org/news>
        <http://citeseer.nj.nec.com/83502.html>

Purpose

  Qhull is a general dimension convex hull program that reads a set 
  of points from stdin, and outputs the smallest convex set that contains 
  the points to stdout.  It also generates Delaunay triangulations, Voronoi 
  diagrams, furthest-site Voronoi diagrams, and halfspace intersections
  about a point.  

  Rbox is a useful tool in generating input for Qhull; it generates 
  hypercubes, diamonds, cones, circles, simplices, spirals, 
  lattices, and random points.
  
  Qhull produces graphical output for Geomview.  This helps with
  understanding the output. <http://www.geomview.org>

    
Environment requirements

  Qhull and rbox should run on all 32-bit and 64-bit computers.  Use
  an ANSI C or C++ compiler to compile the program.  The software is 
  self-contained.  
  
  Qhull is copyrighted software.  Please read COPYING.txt and REGISTER.txt
  before using or distributing Qhull.

To cite Qhull, please use

  Barber, C.B., Dobkin, D.P., and Huhdanpaa, H.T., "The Quickhull 
  algorithm for convex hulls," ACM Trans. on Mathematical Software,
  22(4):469-483, Dec 1996, http://www.qhull.org.

To contribute to Qhull

  Qhull is on Savannah at http://savannah.nongnu.org/projects/qhull/

Qhull on Windows 95, 98, ME, NT, 2000, XP

  The zip file contains rbox.exe, qhull.exe, qconvex.exe, qdelaunay.exe, 
  qhalf.exe, qvoronoi.exe, documentation files, and source files.
  
  To install Qhull:
  - Unzip the files into a directory.  You may use WinZip32 <www.hotfiles.com>
  - Click on QHULL-GO
    
  - In Windows 95, the DOS window needs improvement.
      - Increase the size of the screen font to 8x12.
      - If the text is too dim, fix the screen colors with shareware (e.g., crt.exe)
  - If you use qhull a lot, consider using the Cygwin Unix shell (www.cygwin.com),

  To learn about Qhull:
  - Execute 'qconvex' for a synopsis and examples.
  - Execute 'rbox 10 | qconvex' to compute the convex hull of 10 random points.
  - Execute 'rbox 10 | qconvex i TO file' to write results to 'file'. 
  - If an error occurs, Windows 95 sends the error to stdout instead of stderr 
      - use 'TO xxx' to send normal output to xxx and error output to stdout
  - Browse the documentation: qhull\html\index.htm

Compiling with cygwin on Windows NT, 2000, XP
  - install cygwin [www.cygwin.com] with gcc, make, ar, and ln
  - cd qhull/src
  - make -f Makefile.txt

Qhull on Unix (Debian)

  The gzip file, qhull.tar.gz, contains documentation and source files for
  qhull and rbox.  It should compile on all Unix systems, including Debian.
  You may also use the source instructions below.
  
  To unpack the gzip file
  - tar zxf qhull.tar.gz
  - cd qhull
  
  Compile with the configure Makefile [R. Laboissiere]:
  - ./configure
  - make

Compiling the source distribution

  The gzip file, qhull-src.tgz, contains documentation and source files for
  qhull and rbox.  
  
  To unpack the gzip file
  - tar zxf qhull-src.tgz
  - cd qhull
  
  Compiling with Makefile (i.e., Makefile.txt)   
  - cd src
  - in Makefile, check the CC, CCOPTS1, PRINTMAN, and PRINTC defines
      - the defaults are gcc and enscript
      - CCOPTS1 should include the ANSI flag.  It defines __STDC__
  - in user.h, check the definitions of qh_SECticks and qh_CPUclock.
      - use '#define qh_CLOCKtype 2' for timing runs longer than 1 hour
  - type: make 
      - this builds: qhull qconvex qdelaunay qhalf qvoronoi rbox libqhull.a
  - type: make doc
      - this prints the man page
      - See also qhull/html/index.htm
  - if your compiler reports many errors, it is probably not a ANSI C compiler
      - you will need to set the -ansi switch or find another compiler
  - if your compiler warns about missing prototypes for fprintf() etc.
      - this is ok, your compiler should have these in stdio.h
  - if your compiler warns about missing prototypes for memset() etc.
      - include memory.h in qhull_a.h
  - if your compiler is gcc-2.95.1, you need to set flag -fno-strict-aliasing.  
      - This flag is set by default for other versions [Karas, Krishnaswami]
  - if your compiler reports "global.c: storage size of 'qh_qh' isn't known"
      - delete the initializer "={0}" in global.c, stat.c and mem.c
  - if your compiler warns about "stat.c: improper initializer"
      - this is ok, the initializer is not used
  - if you have trouble building libqhull.a with 'ar'
      - try 'make -f Makefile.txt qhullx' 
  - if the code compiles, the qhull test case will automatically execute
  - if an error occurs, there's an incompatibility between machines
      - For gcc-2.95.1, you need to set flag -fno-strict-aliasing.
	    It is set by default for other versions of gcc [Karas, Krishnaswami]
      - If you can, try a different compiler 
      - You can turn off the Qhull memory manager with qh_NOmem in mem.h
      - You can turn off compiler optimization (-O2 in Makefile)
      - If you find the source of the problem, please let us know
  - if you have Geomview (www.geomview.org)
       - try  'rbox 100 | qconvex G >a' and load 'a' into Geomview
       - run 'q_eg' for Geomview examples of Qhull output (see qh-eg.htm)
  - to install the programs and their man pages:
      - define MANDIR and BINDIR
      - type 'make install'

Compiling on Windows 95, 98, NT, 2000, XP

  Qhull compiles as a console application in Visual C++ 5.0 at warning 
  level 3.

  Visual C++ quickstart for qhull.exe only:
    - create a "Win32 console application" called "qhull"
	- add the following files:
	    geom.c geom2.c global.c io.c mem.c merge.c poly.c poly2.c qhull.c
		qset.c stat.c unix.c user.c
    - create a "Win32 console application" called "rbox" 
	- add rbox.c

  Visual C++ quickstart for qhull library, qhull.exe, qconvex.exe, etc.
    - To simplify setting up lots of projects, 
	- create a temporary "Win32 console application" called "source"
	- add all .c files from .../src/...
	- In Tools::Options::Tab
	  Set tab size to 8 and indent size to 2

    - create a "Win32 console application" called "rbox"
	- move rbox.c from "qhull source"
	- for Project:Settings..., Link
	  you only need the default libraries
	- build the project

    - create a "Win32 static library" called "library"
	- move these files from "qhull source"
	    geom.c geom2.c global.c io.c mem.c merge.c poly.c poly2.c qhull.c
		qset.c stat.c user.c
	- set the library file (use the same for debug and release)
	- build the project

    - create a "Win32 console application" called "qhull"
	- Move unix.c from "qhull source"
	- Add the library file created by "library"
	- Qhull does not use other libraries

    - create a "Win32 console application" called "qconvex"
	- Move qconvex.c from "qhull source"
	- Copy the library file from "qhull"

    - do the same for qdelaun.c, qhalf, qvoronoi.c, user_eg.c, user_eg2.c
	- delete "qhull sources" since it is no longer needed
	- use Project:Settings to make any changes
	- use batch build to rebuild everything
  
  Qhull compiles with Borland C++ 5.0 bcc32.  A Makefile is included.
  Execute 'make -f Mborland'.  If you use the Borland IDE, set the ANSI
  option in Options:Project:Compiler:Source:Language-compliance.
  
  Qhull compiles with Borland C++ 4.02 for Win32 and DOS Power Pack.  
  Use 'make -f Mborland -D_DPMI'.  Qhull 1.0 compiles with Borland 
  C++ 4.02.  For rbox 1.0, use "bcc32 -WX -w- -O2-e -erbox -lc rbox.c".  
  Use the same options for Qhull 1.0. [D. Zwick]
  
  Qhull compiles with Metrowerks C++ 1.7 with the ANSI option.

  If you turn on full warnings, the compiler will report a number of 
  unused variables, variables set but not used, and dead code.  These are
  intentional.  For example, variables may be initialized (unnecessarily)
  to prevent warnings about possible use of uninitialized variables.  

Compiling on the Power Macintosh

  Qhull compiles on the Power Macintosh with Metrowerk's C compiler.
  It uses the SIOUX interface to read point coordinates and return output.
  There is no graphical output.  For project files, see 'Compiling for
  Windows 95'.  Instead of using SIOUX, Qhull may be embedded within an 
  application.  

  Version 1 is available for Macintosh computers by download of qhull.sit.hqx
  It reads point coordinates from a standard file and returns output
  to a standard file.  There is no graphical output.


Compiling on other machines
 
  Some users have reported problems with compiling Qhull under Irix 5.1.  It
  compiles under other versions of Irix. 
  
  If you have troubles with the memory manager, you can turn it off by
  defining qh_NOmem in mem.h.

  You may compile Qhull with a C++ compiler.  


Distributed files

  README.txt           // instructions for installing Qhull 
  REGISTER.txt         // Qhull registration 
  COPYING.txt          // copyright notice 
  QHULL-GO.pif	       // Windows icon for qhull-go.bat
  Announce.txt         // announcement 
  Changes.txt          // change history for Qhull and rbox 
  File_id.diz	       // package descriptor
  index.htm            // Home page 
  html/qh-faq.htm      // Frequently asked questions
  html/qh-get.htm      // Download page
  html/index.htm       // Manual
  src/Makefile.txt     // Makefile for Unix or cygwin 'make' 
  src/Mborland         // Makefile for Borland C++/Win32
  src/Make-config.sh   // Create Debian configure and automake
 
src/      
  rbox consists of:
     rbox.exe          // Win32 executable (.zip only) 
     rbox.htm          // html manual 
     rbox.man          // Unix man page 
     rbox.txt
     rbox.c            // source program 
     
  qhull consists of:
     qhull.exe         // Win32 executables (.zip only) 
     qconvex.exe
     qdelaunay.exe
     qhalf.exe
     qvoronoi.exe
     qhull-go.bat      // DOS window
     qconvex.htm       // html manuals
     qdelaun.htm
     qdelau_f.htm        
     qhalf.htm
     qvoronoi.htm
     qvoron_f.htm
     qh-eg.htm
     qh-impre.htm
     qh-in.htm
     index.htm
     qh-opt*.htm
     qh-quick.htm
     qh--4d.gif,etc.   // images for manual 
     qhull.man         // Unix man page 
     qhull.txt
     q_eg              // shell script for Geomview examples
     q_egtest          // shell script for Geomview test examples
     q_test            // shell script to test qhull
  
  top-level source files:
     src/index.htm     // index to source files 
     qh-...htm         //   specific files
     user.h            // header file of user definable constants 
     qhull.h           // header file for qhull 
     unix.c            // Unix front end to qhull 
     qhull.c           // Quickhull algorithm with partitioning 
     user.c            // user re-definable functions 
     user_eg.c         // example of incorporating qhull into a user program 
     user_eg2.c        // more complex example 
     qhull_interface.cpp // call Qhull from C++

  other source files:
     qhull_a.h         // include file for *.c 
     geom.c            // geometric routines 
     geom2.c
     geom.h	
     global.c          // global variables 
     io.c              // input-output routines 
     io.h
     mem.c             // memory routines, this is stand-alone code 
     mem.h
     merge.c           // merging of non-convex facets 
     merge.h
     poly.c            // polyhedron routines 
     poly2.c
     poly.h 
     qset.c            // set routines, this only depends on mem.c 
     qset.h
     stat.c            // statistics 
     stat.h

Authors:

  C. Bradford Barber                    Hannu Huhdanpaa
  bradb@qhull.org                       hannu@qhull.org
  
                    The Geometry Center
                    University of Minnesota
  
  Qhull 1.0 and 2.0 were developed under NSF grants NSF/DMS-8920161 
  and NSF-CCR-91-15793 750-7504 at the Geometry Center and Harvard 
  University.  If you find Qhull useful, please let us know.
