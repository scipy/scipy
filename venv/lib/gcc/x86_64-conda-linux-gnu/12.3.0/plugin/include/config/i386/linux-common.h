/* Definitions for Intel 386 running Linux-based GNU systems with ELF format.
   Copyright (C) 2012-2022 Free Software Foundation, Inc.
   Contributed by Ilya Enkovich.

This file is part of GCC.

GCC is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3, or (at your option)
any later version.

GCC is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GCC; see the file COPYING3.  If not see
<http://www.gnu.org/licenses/>.  */

#undef TARGET_OS_CPP_BUILTINS
#define TARGET_OS_CPP_BUILTINS()               \
  do                                           \
    {                                          \
      GNU_USER_TARGET_OS_CPP_BUILTINS();       \
      ANDROID_TARGET_OS_CPP_BUILTINS();	       \
    }                                          \
  while (0)

#define EXTRA_TARGET_D_OS_VERSIONS()		\
  ANDROID_TARGET_D_OS_VERSIONS();

#undef CC1_SPEC
#define CC1_SPEC \
  LINUX_OR_ANDROID_CC (GNU_USER_TARGET_CC1_SPEC, \
		       GNU_USER_TARGET_CC1_SPEC " " ANDROID_CC1_SPEC)

#undef	LINK_SPEC
#define LINK_SPEC \
  LINUX_OR_ANDROID_LD (GNU_USER_TARGET_LINK_SPEC, \
		       GNU_USER_TARGET_LINK_SPEC " " ANDROID_LINK_SPEC)

#undef  LIB_SPEC
#define LIB_SPEC \
  LINUX_OR_ANDROID_LD (GNU_USER_TARGET_LIB_SPEC, \
		    GNU_USER_TARGET_NO_PTHREADS_LIB_SPEC " " ANDROID_LIB_SPEC)

#undef  STARTFILE_SPEC
#define STARTFILE_SPEC \
  LINUX_OR_ANDROID_LD (GNU_USER_TARGET_STARTFILE_SPEC, \
		       ANDROID_STARTFILE_SPEC)

#undef  ENDFILE_SPEC
#define ENDFILE_SPEC \
  LINUX_OR_ANDROID_LD (GNU_USER_TARGET_MATHFILE_SPEC " " \
		       GNU_USER_TARGET_ENDFILE_SPEC,	 \
		       GNU_USER_TARGET_MATHFILE_SPEC " " \
		       ANDROID_ENDFILE_SPEC)

#ifdef HAVE_LD_PUSHPOPSTATE_SUPPORT
#define MPX_LD_AS_NEEDED_GUARD_PUSH "--push-state --no-as-needed"
#define MPX_LD_AS_NEEDED_GUARD_POP "--pop-state"
#else
#define MPX_LD_AS_NEEDED_GUARD_PUSH ""
#define MPX_LD_AS_NEEDED_GUARD_POP ""
#endif

extern void file_end_indicate_exec_stack_and_gnu_property (void);

#undef TARGET_ASM_FILE_END
#define TARGET_ASM_FILE_END file_end_indicate_exec_stack_and_gnu_property
