/* Definitions for AMD x86-64 running Linux-based GNU systems with ELF format.
   Copyright (C) 2001-2022 Free Software Foundation, Inc.
   Contributed by Jan Hubicka <jh@suse.cz>, based on linux.h.

This file is part of GCC.

GCC is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3, or (at your option)
any later version.

GCC is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Under Section 7 of GPL version 3, you are granted additional
permissions described in the GCC Runtime Library Exception, version
3.1, as published by the Free Software Foundation.

You should have received a copy of the GNU General Public License and
a copy of the GCC Runtime Library Exception along with this program;
see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
<http://www.gnu.org/licenses/>.  */

#define GNU_USER_LINK_EMULATION32 "elf_i386"
#define GNU_USER_LINK_EMULATION64 "elf_x86_64"
#define GNU_USER_LINK_EMULATIONX32 "elf32_x86_64"

#define GLIBC_DYNAMIC_LINKER32 "/lib/ld-linux.so.2"
#define GLIBC_DYNAMIC_LINKER64 "/lib64/ld-linux-x86-64.so.2"
#define GLIBC_DYNAMIC_LINKERX32 "/libx32/ld-linux-x32.so.2"

#undef MUSL_DYNAMIC_LINKER32
#define MUSL_DYNAMIC_LINKER32 "/lib/ld-musl-i386.so.1"
#undef MUSL_DYNAMIC_LINKER64
#define MUSL_DYNAMIC_LINKER64 "/lib/ld-musl-x86_64.so.1"
#undef MUSL_DYNAMIC_LINKERX32
#define MUSL_DYNAMIC_LINKERX32 "/lib/ld-musl-x32.so.1"
