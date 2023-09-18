/* Dummy definitions of VxWorks-related macros
   Copyright (C) 2007-2022 Free Software Foundation, Inc.

This file is part of GCC.

GCC is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation; either version 3, or (at your option) any later
version.

GCC is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

Under Section 7 of GPL version 3, you are granted additional
permissions described in the GCC Runtime Library Exception, version
3.1, as published by the Free Software Foundation.

You should have received a copy of the GNU General Public License and
a copy of the GCC Runtime Library Exception along with this program;
see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
<http://www.gnu.org/licenses/>.  */

/* True if we're targeting VxWorks, VxWorks7 and/or 64bit.  */
#ifndef TARGET_VXWORKS
#define TARGET_VXWORKS 0
#endif

#ifndef TARGET_VXWORKS7
#define TARGET_VXWORKS7 0
#endif

#ifndef TARGET_VXWORKS64
#define TARGET_VXWORKS64 0
#endif

/* True if generating code for a VxWorks RTP.  */
#ifndef TARGET_VXWORKS_RTP
#define TARGET_VXWORKS_RTP false
#endif

/* The symbol that points to an RTP's table of GOTs.  */
#define VXWORKS_GOTT_BASE (gcc_unreachable (), "")

/* The symbol that holds the index of the current module's GOT in
   VXWORKS_GOTT_BASE.  */
#define VXWORKS_GOTT_INDEX (gcc_unreachable (), "")
