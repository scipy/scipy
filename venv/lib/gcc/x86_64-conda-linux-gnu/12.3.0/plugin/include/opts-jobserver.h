/* GNU make's jobserver related functionality.
   Copyright (C) 2022 Free Software Foundation, Inc.

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
<http://www.gnu.org/licenses/>.

See dbgcnt.def for usage information.  */

#ifndef GCC_JOBSERVER_H
#define GCC_JOBSERVER_H

using namespace std;

struct jobserver_info
{
  /* Default constructor.  */
  jobserver_info ();

  /* Error message if there is a problem.  */
  string error_msg = "";
  /* Skipped MAKEFLAGS where --jobserver-auth is skipped.  */
  string skipped_makeflags = "";
  /* File descriptor for reading used for jobserver communication.  */
  int rfd = -1;
  /* File descriptor for writing used for jobserver communication.  */
  int wfd = -1;
  /* Named pipe path.  */
  string pipe_path = "";
  /* Return true if jobserver is active.  */
  bool is_active = false;
};

#endif /* GCC_JOBSERVER_H */
