Before regenerating the wrappers with SWIG, ensure that you
are using a SWIG version working with Python 3.

- Swig 2.0.1 with the patch from
  http://sourceforge.net/tracker/?func=detail&aid=3047039&group_id=1645&atid=301645

- Also Swig >= 2.0.4 will probably work.

The wrappers are generated with the following commands:
   swig -c++ -python csr.i
   swig -c++ -python csc.i
   swig -c++ -python coo.i
   swig -c++ -python dia.i
   swig -c++ -python bsr.i
   swig -c++ -python csgraph.i
