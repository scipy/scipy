Name

      qhull, rbox         2012.1     2012/02/18
  
Convex hull, Delaunay triangulation, Voronoi diagrams, Halfspace intersection
 
      Documentation:
        html/index.htm
        http://www.qhull.org/html

      Available from:
        <http://www.qhull.org>
        <git@gitorious.org:qhull/qhull.git>
        <http://packages.debian.org/sid/libqhull5> [out-of-date]
 
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
  from Nokia's Qt Framework.  Qhull's C++ interface may change without 
  notice.  Eventually, it will move into the qhull shared library.
  
  Qhull is copyrighted software.  Please read COPYING.txt and REGISTER.txt
  before using or distributing Qhull.

To cite Qhull, please use

  Barber, C.B., Dobkin, D.P., and Huhdanpaa, H.T., "The Quickhull 
  algorithm for convex hulls," ACM Trans. on Mathematical Software,
  22(4):469-483, Dec 1996, http://www.qhull.org.

To contribute to Qhull

  Qhull is on Gitorious (http://gitorious.org:qhull, git@gitorious.org:qhull/qhull.git)
  
  For internal documentation, see html/qh-code.htm

To install Qhull

  Qhull is precompiled for Windows, otherwise it needs compilation.
  
  Besides makefiles for gcc, qhull includes CMakeLists.txt for CMake, 
  vcproj/sln files for Microsoft Visual Studio, and .pro files for Qt Creator.
  It compiles with mingw.
  
  Install and build instructions follow.  
  
  See the end of this document for a list of distributed files.

-----------------
Installing Qhull on Windows

  The zip file contains rbox.exe, qhull.exe, qconvex.exe, qdelaunay.exe, 
  qhalf.exe, qvoronoi.exe, testqset.exe, user_eg*.exe, documentation files, 
  and source files.
  
  To install Qhull:
  - Unzip the files into a directory.  You may use WinZip32 <www.hotfiles.com>
  - Click on QHULL-GO or open a command window into Qhull's bin directory.
    
  To uninstall Qhull
  - Delete the qhull directory
  
  To learn about Qhull:
  - Execute 'qconvex' for a synopsis and examples.
  - Execute 'rbox 10 | qconvex' to compute the convex hull of 10 random points.
  - Execute 'rbox 10 | qconvex i TO file' to write results to 'file'.
  - Browse the documentation: qhull\html\index.htm
  - If an error occurs, Windows sends the error to stdout instead of stderr.
    Use 'TO xxx' to send normal output to xxx and error output to stdout

  To improve the command window
  - Double-click the window bar to increase the size of the window
  - Right-click the window bar
  - Select Properties
  - Check QuickEdit Mode
    Select text with right-click or Enter
    Paste text with right-click
  - Change Font to Lucinda Console
  - Change Layout to Screen Buffer Height 999, Window Size Height 55
  - Change Colors to Screen Background White, Screen Text Black
  - Click OK
  - Select 'Modify shortcut that started this window', then OK

  If you use qhull a lot, install MSYS (www.mingw.org),
    Road Bash (www.qhull.org/bash), or Cygwin (www.cygwin.com).

-----------------
Installing Qhull on Unix with gcc

  To build Qhull, static libraries, shared library, and C++ interface
  - Extract Qhull from qhull...tgz or qhull...zip
  - make
  - export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH

  Or, to build Qhull and libqhullstatic.a
  - Extract Qhull from qhull...tgz or qhull...zip
  - cd src/libqhull
  - make

  The Makefiles may be edited for other compilers.
  If 'testqset' exits with an error, qhull is broken
  
-----------------
Installing Qhull with CMake 2.6 or later

  To build Qhull, static libraries, shared library, and C++ interface
  - Extract Qhull from qhull...tgz or qhull...zip
  - cd build
  - cmake ..
  - make
  - make install
  
  On Windows, CMake installs to C:/Program Files/qhull
  See CMakeLists.txt for further build instructions

-----------------
Installing Qhull with Qt

  To build Qhull, static libraries, shared library, C++ interface, and C++ test
  - Extract Qhull from qhull...tgz or qhull...zip
  - cd src
  - qmake
  - make
  
-----------------
Installing Qhull with Autoconf [WARNING out-of-date]

  The tar.gz tarball contains documentation, source files, 
  and a config directory [R. Laboissiere].

  [Nov 2011] Qhull 2009.1.2 does not include the C++ interface
  
  To install Qhull
  - Extract the files
  - ./configure
  - make
  - make install

-------------------
Working with Qhull's C++ interface

  Qhull's C++ interface is likely to change.  Stay current with Gitorious.

  To clone Qhull's next branch from http://gitorious.org/qhull
    git init
    git clone git://gitorious.org/qhull/qhull.git
    cd qhull
    git checkout next
    ...
    git pull origin next

------------------
Compiling Qhull with Microsoft Visual C++ 2005 or later

  To compile Qhull with Microsoft Visual C++
  - Extract Qhull from Gitorious, qhull...tgz, or qhull...zip
  - Load solution build/qhull.sln 
  - Build
  - Project qhulltest requires Qt for DevStudio (http://qt.nokia.com/downloads)
    Set the QTDIR environment variable to your Qt directory (e.g., c:/qt/4.7.4)
    If incorrect, precompile will fail with 'Can not locate the file specified'

-----------------
Compiling Qhull with Qt Creator

  Qt (http://qt.nokia.com) is a C++ framework for Windows, Linux, and Macintosh

  Qhull uses QTestLib to test qhull's C++ interface (qhulltest)
  
  To compile Qhull with Qt Creator
  - Extract Qhull from Gitorious, qhull...tgz, or qhull...zip
  - Download the Qt SDK from Nokia (http://qt.nokia.com/downloads)
  - Start Qt Creator
  - Load src/qhull-all.pro
  - Build

-----------------
Compiling Qhull with mingw on Windows

  To compile Qhull with MINGW
  - Extract Qhull from Gitorious, qhull...tgz, or qhull...zip
  - Install Road Bash (http://www.qhull.org/bash)
    or install MSYS (http://www.mingw.org/wiki/msys)
  - Install MINGW (http://www.mingw.org/).  Mingw is included with Qt SDK.  
  - make
  
-----------------
Compiling Qhull with cygwin on Windows

  To compile Qhull with cygwin
  - Extract Qhull from Gitorious, qhull...tgz, or qhull...zip
  - Install cygwin (http://www.cygwin.com)
  - Include packages for gcc, make, ar, and ln
  - make

-----------------
Compiling from Makfile without gcc

  The file, qhull-src.tgz, contains documentation and source files for
  qhull and rbox.  
  
  To unpack the gzip file
  - tar zxf qhull-src.tgz
  - cd qhull
  
  Compiling qhull and rbox with Makefile
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
  - if your compiler reports "global.c: storage size of 'qh_qh' isn't known"
      - delete the initializer "={0}" in global.c, stat.c and mem.c
  - if your compiler warns about "stat.c: improper initializer"
      - this is ok, the initializer is not used
  - if you have trouble building libqhull.a with 'ar'
      - try 'make -f Makefile.txt qhullx' 
  - if the code compiles, the qhull test case will automatically execute
  - if an error occurs, there's an incompatibility between machines
      - If you can, try a different compiler 
      - You can turn off the Qhull memory manager with qh_NOmem in mem.h
      - You can turn off compiler optimization (-O2 in Makefile)
      - If you find the source of the problem, please let us know
  - to install the programs and their man pages:
      - define MANDIR and BINDIR
      - type 'make install'

  - if you have Geomview (www.geomview.org)
       - try  'rbox 100 | qconvex G >a' and load 'a' into Geomview
       - run 'q_eg' for Geomview examples of Qhull output (see qh-eg.htm)

------------------
Compiling on other machines and compilers

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

  Qhull compiles on the Power Macintosh with Metrowerk's C compiler.
  It uses the SIOUX interface to read point coordinates and return output.
  There is no graphical output.  For project files, see 'Compiling a
  custom build'.  Instead of using SIOUX, Qhull may be embedded within an
  application.  

  Some users have reported problems with compiling Qhull under Irix 5.1.  It
  compiles under other versions of Irix. 
  
  If you have troubles with the memory manager, you can turn it off by
  defining qh_NOmem in mem.h.

-----------------
Distributed files

  README.txt           // Instructions for installing Qhull 
  REGISTER.txt         // Qhull registration 
  COPYING.txt          // Copyright notice 
  QHULL-GO.lnk         // Windows icon for eg/qhull-go.bat
  Announce.txt         // Announcement 
  CMakeLists.txt       // CMake build file (2.6 or later)
  File_id.diz          // Package descriptor
  index.htm            // Home page 
  Makefile             // Makefile for gcc and other compilers
  qhull*.md5sum        // md5sum for all files

  bin/*		       // Qhull executables and dll (.zip only)
  build/qhull.sln      // DevStudio solution and project files (2005 or later)
  build/*.vcproj
  config/*             // Autoconf files for creating configure (Unix only)
  eg/*                 // Test scripts and geomview files from q_eg
  html/index.htm       // Manual
  html/qh-faq.htm      // Frequently asked questions
  html/qh-get.htm      // Download page
  html/qhull-cpp.xml   // C++ style notes as a Road FAQ (www.qhull.org/road)
  src/Changes.txt      // Change history for Qhull and rbox 
  src/qhull-all.pro    // Qt project

eg/ 
  q_eg                 // shell script for Geomview examples (eg.01.cube)
  q_egtest             // shell script for Geomview test examples
  q_test               // shell script to test qhull
  q_test-ok.txt        // output from q_test
  qhulltest-ok.txt     // output from qhulltest (Qt only)

rbox consists of (bin, html):
  rbox.exe             // Win32 executable (.zip only) 
  rbox.htm             // html manual 
  rbox.man             // Unix man page 
  rbox.txt

qhull consists of (bin, html):
  qhull.exe            // Win32 executables and dlls (.zip only) 
  qconvex.exe
  qdelaunay.exe
  qhalf.exe
  qvoronoi.exe
  qhull.dll
  qhull_p.dll
  qhull-go.bat         // command window
  qconvex.htm          // html manual
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
  qh--*.gif            // images for manual
  normal_voronoi_knauss_oesterle.jpg
  qhull.man            // Unix man page 
  qhull.txt

bin/
  msvcr80.dll          // Visual C++ redistributable file (.zip only)

src/
  qhull/unix.c         // Qhull and rbox applications
  qconvex/qconvex.c    
  qhalf/qhalf.c
  qdelaunay/qdelaunay.c
  qvoronoi/qvoronoi.c
  rbox/rbox.c

  user_eg/user_eg.c    // example of using qhull_p.dll (requires -Dqh_QHpointer)
  user_eg2/user_eg2.c  // example of using qhull.dll from a user program
  user_eg3/user_eg3.cpp // example of Qhull's C++ interface with libqhullstatic_p.a
  qhulltest/qhulltest.cpp // Test of Qhull's C++ interface using Qt's QTestLib
  qhull-*.pri          // Include files for Qt projects

src/libqhull
  libqhull.pro         // Qt project for shared library (qhull.dll)
  index.htm            // design documentation for libqhull
  qh-*.htm
  qhull-exports.def    // Export Definition file for Visual C++
  Makefile             // Simple gcc Makefile for qhull and libqhullstatic.a
  Mborland             // Makefile for Borland C++ 5.0

  libqhull.h           // header file for qhull
  user.h               // header file of user definable constants 
  libqhull.c           // Quickhull algorithm with partitioning
  user.c               // user re-definable functions 
  usermem.c
  userprintf.c
  userprintf_rbox.c

  qhull_a.h            // include files for libqhull/*.c 
  geom.c               // geometric routines 
  geom2.c
  geom.h     
  global.c             // global variables 
  io.c                 // input-output routines 
  io.h   
  mem.c                // memory routines, this is stand-alone code 
  mem.h
  merge.c              // merging of non-convex facets 
  merge.h
  poly.c               // polyhedron routines 
  poly2.c
  poly.h 
  qset.c               // set routines, this only depends on mem.c 
  qset.h
  random.c             // utilities w/ Park & Miller's random number generator
  random.h
  rboxlib.c            // point set generator for rbox
  stat.c               // statistics 
  stat.h

src/libqhullp
  libqhullp.pro        // Qt project for shared library (qhull_p.dll)
  qhull_p-exports.def  // Export Definition file for Visual C++

src/libqhullstatic/
  libqhullstatic.pro   // Qt project for static library     
     
src/libqhullstaticp/
  libqhullstaticp.pro  // Qt project for static library with qh_QHpointer
     
src/libqhullcpp/
  libqhullcpp.pro      // Qt project for static C++ library     
  Qhull.cpp            // Call libqhull.c from C++
  Qhull.h
  qt-qhull.cpp         // Supporting methods for Qt
  qhull_interface.cpp  // Another approach to C++
    
  Coordinates.cpp      // input classes
  Coordinates.h
  PointCoordinates.cpp
  PointCoordinates.h
  RboxPoints.cpp       // call rboxlib.c from C++
  RboxPoints.h

  QhullFacet.cpp       // data structure classes
  QhullFacet.h
  QhullHyperplane.cpp
  QhullHyperplane.h
  QhullPoint.cpp
  QhullPoint.h
  QhullQh.cpp
  QhullStat.cpp
  QhullStat.h
  QhullVertex.cpp
  QhullVertex.h
  
  QhullFacetList.cpp   // collection classes
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
  QhullVertexSet.cpp
  QhullVertexSet.h

  functionObjects.h    // supporting classes
  QhullError.cpp
  QhullError.h
  QhullQh.cpp
  QhullQh.h
  UsingLibQhull.cpp
  UsingLibQhull.h

src/qhulltest/
  qhulltest.pro        // Qt project for test of C++ interface     
  Coordinates_test.cpp // Test of each class
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
  QhullVertexSet_test.cpp
  QhullVertex_test.cpp
  Qhull_test.cpp
  RboxPoints_test.cpp
  UsingLibQhull_test.cpp

src/road/
  RoadError.cpp        // Supporting base classes
  RoadError.h
  RoadLogEvent.cpp
  RoadLogEvent.h
  RoadTest.cpp         // Run multiple test files with QTestLib
  RoadTest.h

src/testqset/
  testqset.pro		// Qt project for test qset.c with mem.c
  testqset.c
  
-----------------
Authors:

  C. Bradford Barber                  Hannu Huhdanpaa (Version 1.0)
  bradb@shore.net                     hannu@qhull.org
  
  Qhull 1.0 and 2.0 were developed under NSF grants NSF/DMS-8920161 
  and NSF-CCR-91-15793 750-7504 at the Geometry Center and Harvard 
  University.  If you find Qhull useful, please let us know.
