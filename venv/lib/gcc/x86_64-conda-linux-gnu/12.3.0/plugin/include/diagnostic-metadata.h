/* Additional metadata for a diagnostic.
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

#ifndef GCC_DIAGNOSTIC_METADATA_H
#define GCC_DIAGNOSTIC_METADATA_H

/* A bundle of additional metadata that can be associated with a
   diagnostic.

   Currently this only supports associating a CWE identifier with a
   diagnostic.  */

class diagnostic_metadata
{
 public:
  diagnostic_metadata () : m_cwe (0) {}

  void add_cwe (int cwe) { m_cwe = cwe; }
  int get_cwe () const { return m_cwe; }

 private:
  int m_cwe;
};

#endif /* ! GCC_DIAGNOSTIC_METADATA_H */
