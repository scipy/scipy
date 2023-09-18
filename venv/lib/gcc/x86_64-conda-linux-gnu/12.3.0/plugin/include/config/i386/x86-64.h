/* OS independent definitions for AMD x86-64.
   Copyright (C) 2001-2022 Free Software Foundation, Inc.
   Contributed by Bo Thorsen <bo@suse.de>.

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

#undef ASM_COMMENT_START
#define ASM_COMMENT_START "#"

#undef DBX_REGISTER_NUMBER
#define DBX_REGISTER_NUMBER(n) \
  (TARGET_64BIT ? dbx64_register_map[n] : svr4_dbx_register_map[n])

/* Output assembler code to FILE to call the profiler.  */
#define NO_PROFILE_COUNTERS 1

#undef MCOUNT_NAME
#define MCOUNT_NAME "mcount"

#undef SIZE_TYPE
#define SIZE_TYPE (TARGET_LP64 ? "long unsigned int" : "unsigned int")

#undef PTRDIFF_TYPE
#define PTRDIFF_TYPE (TARGET_LP64 ? "long int" : "int")

#undef WCHAR_TYPE
#define WCHAR_TYPE "int"

#undef WCHAR_TYPE_SIZE
#define WCHAR_TYPE_SIZE 32

#undef ASM_SPEC
#define ASM_SPEC "%{m16|m32:--32} %{m64:--64} %{mx32:--x32}"

#undef ASM_OUTPUT_ALIGNED_BSS
#define ASM_OUTPUT_ALIGNED_BSS(FILE, DECL, NAME, SIZE, ALIGN) \
  x86_output_aligned_bss (FILE, DECL, NAME, SIZE, ALIGN)

#undef  ASM_OUTPUT_ALIGNED_DECL_COMMON
#define ASM_OUTPUT_ALIGNED_DECL_COMMON(FILE, DECL, NAME, SIZE, ALIGN)		\
  x86_elf_aligned_decl_common (FILE, DECL, NAME, SIZE, ALIGN);

#undef  ASM_OUTPUT_ALIGNED_DECL_LOCAL
#define ASM_OUTPUT_ALIGNED_DECL_LOCAL(FILE, DECL, NAME, SIZE, ALIGN)  \
  do								      \
    {								      \
      fprintf ((FILE), "%s", LOCAL_ASM_OP);			      \
      assemble_name ((FILE), (NAME));				      \
      fprintf ((FILE), "\n");					      \
      ASM_OUTPUT_ALIGNED_DECL_COMMON (FILE, DECL, NAME, SIZE, ALIGN); \
    }								      \
  while (0)

#define SUBALIGN_LOG 3

/* i386 System V Release 4 uses DWARF debugging info.
   x86-64 ABI specifies DWARF2.  */

#define DWARF2_DEBUGGING_INFO 1
#define DWARF2_UNWIND_INFO 1

#undef PREFERRED_DEBUGGING_TYPE
#define PREFERRED_DEBUGGING_TYPE DWARF2_DEBUG

#undef TARGET_ASM_SELECT_SECTION
#define TARGET_ASM_SELECT_SECTION  x86_64_elf_select_section

#undef TARGET_ASM_UNIQUE_SECTION
#define TARGET_ASM_UNIQUE_SECTION  x86_64_elf_unique_section

#undef TARGET_SECTION_TYPE_FLAGS
#define TARGET_SECTION_TYPE_FLAGS  x86_64_elf_section_type_flags
