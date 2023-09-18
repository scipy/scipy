/* Helper code for graphviz output.
   Copyright (C) 2019-2022 Free Software Foundation, Inc.
   Contributed by David Malcolm <dmalcolm@redhat.com>.

This file is part of GCC.

GCC is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3, or (at your option)
any later version.

GCC is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with GCC; see the file COPYING3.  If not see
<http://www.gnu.org/licenses/>.  */

#ifndef GCC_GRAPHVIZ_H
#define GCC_GRAPHVIZ_H

#include "pretty-print.h" /* for ATTRIBUTE_GCC_PPDIAG.  */

/* A class for writing .dot output to a pretty_printer with
   indentation to show nesting.  */

class graphviz_out {
 public:
  graphviz_out (pretty_printer *pp);

  void print (const char *fmt, ...)
    ATTRIBUTE_GCC_PPDIAG(2,3);
  void println (const char *fmt, ...)
    ATTRIBUTE_GCC_PPDIAG(2,3);

  void indent () { m_indent++; }
  void outdent () { m_indent--; }

  void write_indent ();

  void begin_tr ();
  void end_tr ();

  void begin_td ();
  void end_td ();

  void begin_trtd ();
  void end_tdtr ();

  pretty_printer *get_pp () const { return m_pp; }

 private:
  pretty_printer *m_pp;
  int m_indent;
};

#endif /* GCC_GRAPHVIZ_H */
