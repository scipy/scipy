Name

      qhull, rbox         2010.1     2010/01/14
  
Convex hull, Delaunay triangulation, Voronoi diagrams, Halfspace intersection
 
      Documentation:
        html/index.htm
        http://www.qhull.org/html

      Available from:
        <http://www.qhull.org>
        <git@gitorious.org:qhull/qhull.git>
        <http://packages.debian.org/sid/libqhull5>
 
      News and a paper:
        <http://www.qhull.org/news>
        <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.117.405>

     Version 1 (simplicial only):
        <http://www.qhull.org/download/qhull-1.0.tar.gz>
       

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
  self-contained.  It comes with examples and test scripts.
  
  Qhull's C++ interface uses the STL.  The C++ test program uses QTestLib 
  from Nokia's Qt Framework.  For 2010, Qhull's C++ interface 
  may change without notice.
  
  Qhull is copyrighted software.  Please read COPYING.txt and REGISTER.txt
  before using or distributing Qhull.

To cite Qhull, please use

  Barber, C.B., Dobkin, D.P., and Huhdanpaa, H.T., "The Quickhull 
  algorithm for convex hulls," ACM Trans. on Mathematical Software,
  22(4):469-483, Dec 1996, http://www.qhull.org.

To contribute to Qhull

  Qhull is on Gitorious (http://gitorious.org:qhull, git@gitorious.org:qhull/qhull.git)
  
  For internal documentation, see html/qh-code.htm

-----------------
Installing Qhull on Windows

  The zip file contains rbox.exe, qhull.exe, qconvex.exe, qdelaunay.exe, 
  qhalf.exe, qvoronoi.exe, user_eg3.exe, documentation files, and source files.
  
  To install Qhull:
  - Unzip the files into a directory.  You may use WinZip32 <www.hotfiles.com>
  - Click on QHULL-GO
    
  - If you use qhull a lot, consider using 
    Road Bash (www.qhull.org/bash) or Cygwin (www.cygwin.com),
  - In Windows 95, the DOS window needs improvement.
      - Increase the size of the screen font to 8x12.
      - If the text is too dim, fix the screen colors with shareware (e.g., crt.exe)

  To learn about Qhull:
  - Execute 'qconvex' for a synopsis and examples.
  - Execute 'rbox 10 | qconvex' to compute the convex hull of 10 random points.
  - Execute 'rbox 10 | qconvex i TO file' to write results to 'file'. 
  - Browse the documentation: qhull\html\index.htm
  - If an error occurs, Windows sends the error to stdout instead of stderr. 
    Use 'TO xxx' to send normal output to xxx and error output to stdout

-----------------
Installing Qhull on Unix with gcc

  The tgz tarball contains documentation and source files.
  
  To install Qhull
  - Extract the files
  - cd qhull/src
  - make
  - or, make -f Makefile.txt

-----------------
Installing Qhull with Autoconf

  The tar.gz tarball contains documentation, source files, 
  and a config directory [R. Laboissiere].
  
  To install Qhull
  - Extract the files
  - ./configure
  - make
  - make install

-----------------
Compiling Qhull with Qt 

  Qt is a C++ framework for Windows, Linux, and Macintosh

  Qhull includes a C++ test using Qt's QTestLib

  To compile Qhull  
  - Download the Qt SDK from Nokia (http://qt.nokia.com/downloads)
  - Start Qt Creator
  - Load project/qhull-all.pro
  - Build all

-----------------
Compiling Qhull on Windows.

  To compile Qhull with Visual C++
  - Load solution project/qhull.sln
  - For project qhulltest, 
    install Qt for DevStudio (http://qt.nokia.com/downloads)
  - Build->Build all

  To compile Qhull with MINGW
  - Install Road Bash (http://www.qhull.org/bash)
  - or install MSYS (http://www.mingw.org/wiki/msys)
  - Install MINGW (http://www.mingw.org/)
  - cd src
  - make 
  - or, make -f Makefile.txt

  To compile Qhull with cygwin
  - Install cygwin (http://www.cygwin.com)
  - Include packages for gcc, make, ar, and ln
  - cd src
  - make
  - or, make -f Makefile.txt

-----------------
Compiling the source distribution

  The gzip file, qhull-src.tgz, contains documentation and source files for
  qhull and rbox.  
  
  To unpack the gzip file
  - tar zxf qhull-src.tgz
  - cd qhull
  
  Compiling qhull and rbox with Makefile (i.e., Makefile.txt)   
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
  - if your compiler is gcc-4.3, 4.2, 4.1 or 2.95.1, set flag -fno-strict-aliasing.  
      - This flag is set by default for other versions up to 4.0 [Karas, Krishnaswami]
      - strict-aliasing appears to work ok for gcc-4.4
      - See news/qhull-news.html#bugs for details
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

If you need to create a custom build, create projects as follows.
These instructions assume Visual C++, but similar ones apply to other
build systems.

  For qhull.exe and rbox.exe only:
    - create a "console application" called "qhull"
        - add the following files:
            geom.c geom2.c global.c io.c libqhull.c mem.c merge.c 
            poly.c poly2.c qset.c random.c stat.c unix.c 
            user.c usermem.c userprintf.c
    - create a "console application" called "rbox" 
        - add rbox.c rboxlib.c

  For all projects
    - To simplify setting up lots of projects, 
        - create a temporary "Win32 console application" called "source"
        - add all .c files from .../src/...
        - add all .cpp files from .../cpp/...
        - In Tools::Options::Tab
          Set tab size to 8 and indent size to 2

    - create a "Win32 console application" called "rbox"
        - move rbox.c rboxlib.c from "qhull source"
        - for Project:Settings..., Link
          you only need the default libraries
        - build the project

    - create a "Win32 static library" called "libqhull"
        - move these files from "qhull source"
            geom.c geom2.c global.c io.c mem.c merge.c poly.c poly2.c libqhull.c
                qset.c random.c rboxlib.c stat.c user.c usermem.c userprintf.c
        - set the library file (use the same for debug and release)
        - build the project

    - create a "Win32 static library" called "libqhullp"
        - move these files from "qhull source"
            geom.c geom2.c global.c io.c mem.c merge.c poly.c poly2.c libqhull.c
                qset.c random.c rboxlib.c stat.c user.c usermem.c userprintf.c
        - define qh_QHpointer=1
        - set the library file (use the same for debug and release)
        - build the project

    - create a "Win32 console application, empty project" called "qhull"
        - Move unix.c from "qhull source"
        - Project Dependency on libqhull (or copy the library file)
        - Linker->Input->Ignore LIBCMT
        - Qhull does not use other libraries

    - create a "Win32 console application, empty project" called "qconvex"
        - Move qconvex.c from "qhull source"
        - Project Dependency on libqhull (or copy the library file)
        - Linker->Input->Ignore LIBCMT
        - build the project

    - do the same for qdelaun.c, qhalf, qvoronoi.c, user_eg.c, user_eg2.c

    - create a "Win32 dynamic library" called "libqhullcpp"
        - Move cpp sources in cpp/*.cpp from "qhull source"
          Do not move cpp/qhulltest, road/RoadTest.cpp, or user_eg3.cpp
        - define qh_QHpointer=1
        - Add the library file created by "libqhullp"
        - build the project
    
    - create a "Win32 console application" called "user_eg3"
        - Move user_eg3 from "qhull source"
        - define qh_QHpointer=1
        - Add the library file created by "libqhullcpp"
        - build the project

    - create a "Win32 console application" called "qhulltest"
        - Install Qt for DevStudio (http://qt.nokia.com/downloads)
        - Add everything in cpp/qhulltest
        - Add road/RoadTest.cpp
        - define qh_QHpointer=1
        - Add the library file created by "libqhullcpp"
        - build the project

    - use Project:Settings to make any changes
    - use batch build to rebuild everything
    - delete "qhull sources" since it is no longer needed
  
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

Compiling on other machines
 
  Some users have reported problems with compiling Qhull under Irix 5.1.  It
  compiles under other versions of Irix. 
  
  If you have troubles with the memory manager, you can turn it off by
  defining qh_NOmem in mem.h.

  You may compile Qhull with a C++ compiler.  


-----------------
Distributed files

  README.txt           // instructions for installing Qhull 
  REGISTER.txt         // Qhull registration 
  COPYING.txt          // copyright notice 
  QHULL-GO.pif         // Windows icon for qhull-go.bat
  Announce.txt         // announcement 
  CMakeLists.txt       // CMake file
  File_id.diz          // package descriptor
  index.htm            // Home page 
  qhull*.md5sum        // md5sum for all files
  html/qh-faq.htm      // Frequently asked questions
  html/qh-get.htm      // Download page
  html/index.htm       // Manual
  src/Changes.txt      // change history for Qhull and rbox 
  src/CMakeLists.txt   // CMake file
  src/Makefile.txt     // Makefile for Unix or cygwin 'make' 
  src/Mborland         // Makefile for Borland C++/Win32
  src/Make-config.sh   // Create Debian configure and automake
  project/qhull-all.pro  // Qt project with Visual C++ vcproj files

eg/ 
  q_eg                 // shell script for Geomview examples
  q_egtest             // shell script for Geomview test examples
  q_test               // shell script to test qhull
  q_test-ok.txt        // output from q_test
  qhulltest-ok.txt     // output from qhulltest (Qt only)

src/      
  rbox consists of:
     rbox.exe          // Win32 executable (.zip only) 
     rbox.htm          // html manual 
     rbox.man          // Unix man page 
     rbox.txt
     rbox.c            // source program 
     rboxlib.c
     
  qhull consists of:
     qhull.exe         // Win32 executables (.zip only) 
     qconvex.exe
     qdelaunay.exe
     qhalf.exe
     qvoronoi.exe
     eg/qhull-go.bat   // DOS window
     qconvex.htm       // html manuals
     qdelaun.htm
     qdelau_f.htm        
     qhalf.htm
     qvoronoi.htm
     qvoron_f.htm
     qh-eg.htm
     qh-code.htm
     qh-impre.htm
     index.htm
     qh-opt*.htm
     qh-quick.htm
     qh--4d.gif,etc.   // images for manual 
     qhull.man         // Unix man page 
     qhull.txt
  
  top-level source files:
     src/index.htm     // index to source files 
     qh-...htm         //   specific files
     user.h            // header file of user definable constants 
     libqhull.h        // header file for qhull
     unix.c            // Unix front end to qhull 
     qconvex.c    
     qhalf.c
     qdelaunay.c
     qvoronoi.c
     libqhull.c        // Quickhull algorithm with partitioning
     user.c            // user re-definable functions 
     usermem.c
     userprintf.c
     user_eg.c         // example of incorporating qhull into a user program 
     user_eg2.c        // more complex example 

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

project/
  projects for building Qhull
     qhull.sln          // Solution for Visual C++ 2005 and higher
     qhull-all.pro      // Project for Qt
     */*.pro            // Qt projects for each component
     *.vcproj           // Visual C++ projects 
     
cpp/      
  cpp interface to Qhull
    Qhull.cpp           // calls libqhull.c
    Qhull.h
    user_eg3.cpp        // sample program
    
  input classes
    Coordinates.cpp
    Coordinates.h
    PointCoordinates.cpp
    PointCoordinates.h
    PointIterator.h
    RboxPoints.cpp      // calls rboxlib.c
    RboxPoints.h

  data structure classes  
    QhullFacet.cpp
    QhullFacet.h
    QhullHyperplane.cpp
    QhullHyperplane.h
    QhullPoint.cpp
    QhullPoint.h
    QhullQh.cpp
    QhullQh.h
    QhullStat.cpp
    QhullStat.h
    
  collection classes
    QhullFacetList.cpp
    QhullFacetList.h
    QhullFacetSet.cpp
    QhullFacetSet.h
    QhullIterator.h
    QhullLinkedList.h
    QhullPoints.cpp
    QhullPoints.h
    QhullPointSet.cpp
    QhullPointSet.h
    QhullRidge.cpp
    QhullRidge.h
    QhullSet.cpp
    QhullSet.h
    QhullSets.h

  supporting classes
    functionObjects.h
    QhullError.cpp
    QhullError.h
    QhullEvent.cpp
    QhullEvent.h
    QhullLog.cpp
    QhullLog.h
    UsingLibQhull.cpp
    UsingLibQhull.h

cpp/qhulltest/
  Tests for each class
    Coordinates_test.cpp
    PointCoordinates_test.cpp
    Point_test.cpp
    QhullFacetList_test.cpp
    QhullFacetSet_test.cpp
    QhullFacet_test.cpp
    QhullHyperplane_test.cpp
    QhullLinkedList_test.cpp
    QhullPointSet_test.cpp
    QhullPoints_test.cpp
    QhullPoint_test.cpp
    QhullRidge_test.cpp
    QhullSet_test.cpp
    qhulltest.cpp
    QhullVertexSet_test.cpp
    QhullVertex_test.cpp
    Qhull_test.cpp
    RboxPoints_test.cpp
    UsingLibQhull_test.cpp

cpp/road/
  Supporting base classes
    RoadError.cpp
    RoadError.h
    RoadLogEvent.cpp
    RoadLogEvent.h
    RoadTest.cpp
    RoadTest.h

-----------------
Authors:

  C. Bradford Barber                    Hannu Huhdanpaa
  bradb@shore.net                       hannu@qhull.org
  
  Qhull 1.0 and 2.0 were developed under NSF grants NSF/DMS-8920161 
  and NSF-CCR-91-15793 750-7504 at the Geometry Center and Harvard 
  University.  If you find Qhull useful, please let us know.
