/* Header file for the GIMPLE range tracing/debugging facilties.
   Copyright (C) 2021-2022 Free Software Foundation, Inc.
   Contributed by Andrew MacLeod <amacleod@redhat.com>
   and Aldy Hernandez <aldyh@redhat.com>.

This file is part of GCC.

GCC is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation; either version 3, or (at your option) any later
version.

GCC is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 for more details.

You should have received a copy of the GNU General Public License
along with GCC; see the file COPYING3.  If not see
<http://www.gnu.org/licenses/>.  */

#ifndef GCC_GIMPLE_RANGE_TRACE_H
#define GCC_GIMPLE_RANGE_TRACE_H

// This class manages range tracing for the ranger and gori components.
// Tracing will provide a unique integer index whenever a new trace
// is started. This can be used to identify where a calculation has gone wrong.

class range_tracer
{
public:
  range_tracer (const char *name = "");
  unsigned header (const char *str);
  void trailer (unsigned counter, const char *caller, bool result, tree name,
		const irange &r);
  void print (unsigned counter, const char *str);
  inline void enable_trace () { tracing = true; }
  inline void disable_trace () { tracing = false; }
  virtual void breakpoint (unsigned index);
private:
  unsigned do_header (const char *str);
  void print_prefix (unsigned idx, bool blanks);
  static const unsigned bump = 2;
  unsigned indent;
  static const unsigned name_len = 100;
  char component[name_len];
  bool tracing;
};


// If tracing is enabled, start a new trace header, returning the trace index.
// Otherwise return 0.

inline unsigned
range_tracer::header (const char *str)
{
  if (tracing)
    return do_header (str);
  return 0;
}

// RAII class to change current dump_file and dump_flags, and restore
// when the object goes out of scope.

class push_dump_file
{
public:
  push_dump_file (FILE *, dump_flags_t);
  ~push_dump_file ();
private:
  FILE *old_dump_file;
  dump_flags_t old_dump_flags;
};

void dump_ranger (FILE *);
void dump_ranger (FILE *, const vec<basic_block> &path);

#endif // GCC_GIMPLE_RANGE_TRACE_H
