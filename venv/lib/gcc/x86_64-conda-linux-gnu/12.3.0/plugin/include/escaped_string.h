/* Shared escaped string class.
   Copyright (C) 1999-2022 Free Software Foundation, Inc.

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

#ifndef GCC_ESCAPED_STRING_H
#define GCC_ESCAPED_STRING_H

#include <cstdlib>

/* A class to handle converting a string that might contain
   control characters, (eg newline, form-feed, etc), into one
   in which contains escape sequences instead.  */

class escaped_string
{
 public:
  escaped_string () { m_owned = false; m_str = NULL; };
  ~escaped_string () { if (m_owned) free (m_str); }
  operator const char *() const { return m_str; }
  void escape (const char *);
 private:
  escaped_string(const escaped_string&) {}
  escaped_string& operator=(const escaped_string&) { return *this; }
  char *m_str;
  bool  m_owned;
};

#endif /* ! GCC_ESCAPED_STRING_H */
