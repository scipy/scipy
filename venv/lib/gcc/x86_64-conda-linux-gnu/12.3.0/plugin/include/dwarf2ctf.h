/* dwarf2ctf.h - DWARF interface for CTF/BTF generation.
   Copyright (C) 2021-2022 Free Software Foundation, Inc.

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

/* This file contains declarations and prototypes to define an interface
   between DWARF and CTF/BTF generation.  */

#ifndef GCC_DWARF2CTF_H
#define GCC_DWARF2CTF_H 1

#include "dwarf2out.h"
#include "flags.h"

/* Debug Format Interface.  Used in dwarf2out.cc.  */

extern void ctf_debug_init (void);
extern void ctf_debug_init_postprocess (bool);
extern bool ctf_do_die (dw_die_ref);
extern void ctf_debug_early_finish (const char *);
extern void ctf_debug_finish (const char *);

/* Wrappers for CTF/BTF to fetch information from GCC DWARF DIE.  Used in
   ctfc.cc.

   A CTF container does not store all debug information internally.  Some of
   the info is fetched indirectly via the DIE reference available in each CTF
   container entry.

   These functions will be used by the CTF container to give access to its
   consumers (CTF/BTF) to various debug information available in DWARF DIE.
   Direct access to debug information in GCC dwarf structures by the consumers
   of CTF/BTF information is not ideal.  */

/* Source location information.  */

extern const char * ctf_get_die_loc_file (dw_die_ref);
extern unsigned int ctf_get_die_loc_line (dw_die_ref);
extern unsigned int ctf_get_die_loc_col (dw_die_ref);

#endif /* GCC_DWARF2CTF_H */
