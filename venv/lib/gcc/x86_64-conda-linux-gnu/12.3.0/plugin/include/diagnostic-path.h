/* Paths through the code associated with a diagnostic.
   Copyright (C) 2019-2022 Free Software Foundation, Inc.
   Contributed by David Malcolm <dmalcolm@redhat.com>

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

#ifndef GCC_DIAGNOSTIC_PATH_H
#define GCC_DIAGNOSTIC_PATH_H

#include "diagnostic.h" /* for ATTRIBUTE_GCC_DIAG.  */
#include "diagnostic-event-id.h"

/* A diagnostic_path is an optional additional piece of metadata associated
   with a diagnostic (via its rich_location).

   It describes a sequence of events predicted by the compiler that
   lead to the problem occurring, with their locations in the user's source,
   and text descriptions.

   For example, the following error has a 3-event path:

     test.c: In function 'demo':
     test.c:29:5: error: passing NULL as argument 1 to 'PyList_Append' which
       requires a non-NULL parameter
        29 |     PyList_Append(list, item);
           |     ^~~~~~~~~~~~~~~~~~~~~~~~~
       'demo': events 1-3
          |
          |   25 |   list = PyList_New(0);
          |      |          ^~~~~~~~~~~~~
          |      |          |
          |      |          (1) when 'PyList_New' fails, returning NULL
          |   26 |
          |   27 |   for (i = 0; i < count; i++) {
          |      |   ~~~
          |      |   |
          |      |   (2) when 'i < count'
          |   28 |     item = PyLong_FromLong(random());
          |   29 |     PyList_Append(list, item);
          |      |     ~~~~~~~~~~~~~~~~~~~~~~~~~
          |      |     |
          |      |     (3) when calling 'PyList_Append', passing NULL from (1) as argument 1
          |

    The diagnostic-printing code has consolidated the path into a single
    run of events, since all the events are near each other and within the same
    function; more complicated examples (such as interprocedural paths)
    might be printed as multiple runs of events.  */

/* Abstract base classes, describing events within a path, and the paths
   themselves.  */

/* One event within a diagnostic_path.  */

class diagnostic_event
{
 public:
  virtual ~diagnostic_event () {}

  virtual location_t get_location () const = 0;

  virtual tree get_fndecl () const = 0;

  /* Stack depth, so that consumers can visualizes the interprocedural
     calls, returns, and frame nesting.  */
  virtual int get_stack_depth () const = 0;

  /* Get a localized (and possibly colorized) description of this event.  */
  virtual label_text get_desc (bool can_colorize) const = 0;
};

/* Abstract base class for getting at a sequence of events.  */

class diagnostic_path
{
 public:
  virtual ~diagnostic_path () {}
  virtual unsigned num_events () const = 0;
  virtual const diagnostic_event & get_event (int idx) const = 0;

  bool interprocedural_p () const;
};

/* Concrete subclasses.  */

/* A simple implementation of diagnostic_event.  */

class simple_diagnostic_event : public diagnostic_event
{
 public:
  simple_diagnostic_event (location_t loc, tree fndecl, int depth,
			   const char *desc);
  ~simple_diagnostic_event ();

  location_t get_location () const FINAL OVERRIDE { return m_loc; }
  tree get_fndecl () const FINAL OVERRIDE { return m_fndecl; }
  int get_stack_depth () const FINAL OVERRIDE { return m_depth; }
  label_text get_desc (bool) const FINAL OVERRIDE
  {
    return label_text::borrow (m_desc);
  }

 private:
  location_t m_loc;
  tree m_fndecl;
  int m_depth;
  char *m_desc; // has been i18n-ed and formatted
};

/* A simple implementation of diagnostic_path, as a vector of
   simple_diagnostic_event instances.  */

class simple_diagnostic_path : public diagnostic_path
{
 public:
  simple_diagnostic_path (pretty_printer *event_pp)
  : m_event_pp (event_pp) {}

  unsigned num_events () const FINAL OVERRIDE;
  const diagnostic_event & get_event (int idx) const FINAL OVERRIDE;

  diagnostic_event_id_t add_event (location_t loc, tree fndecl, int depth,
				   const char *fmt, ...)
    ATTRIBUTE_GCC_DIAG(5,6);

 private:
  auto_delete_vec<simple_diagnostic_event> m_events;

  /* (for use by add_event).  */
  pretty_printer *m_event_pp;
};

extern void debug (diagnostic_path *path);

#endif /* ! GCC_DIAGNOSTIC_PATH_H */
