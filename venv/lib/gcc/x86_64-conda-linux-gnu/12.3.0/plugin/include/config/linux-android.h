/* Configuration file for Linux Android targets.
   Copyright (C) 2008-2022 Free Software Foundation, Inc.
   Contributed by Doug Kwan (dougkwan@google.com)
   Rewritten by CodeSourcery, Inc.

   This file is part of GCC.

   GCC is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published
   by the Free Software Foundation; either version 3, or (at your
   option) any later version.

   GCC is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
   License for more details.

   You should have received a copy of the GNU General Public License
   along with GCC; see the file COPYING3.  If not see
   <http://www.gnu.org/licenses/>.  */

#define ANDROID_TARGET_OS_CPP_BUILTINS()			\
    do {							\
	if (TARGET_ANDROID)					\
	  builtin_define ("__ANDROID__");			\
    } while (0)

#define ANDROID_TARGET_D_OS_VERSIONS()				\
    do {							\
	if (TARGET_ANDROID)					\
	  builtin_version ("Android");				\
    } while (0)

#if ANDROID_DEFAULT
# define NOANDROID "mno-android"
#else
# define NOANDROID "!mandroid"
#endif

#define LINUX_OR_ANDROID_CC(LINUX_SPEC, ANDROID_SPEC) \
  "%{" NOANDROID "|tno-android-cc:" LINUX_SPEC ";:" ANDROID_SPEC "}"

#define LINUX_OR_ANDROID_LD(LINUX_SPEC, ANDROID_SPEC) \
  "%{" NOANDROID "|tno-android-ld:" LINUX_SPEC ";:" ANDROID_SPEC "}"

#define ANDROID_LINK_SPEC \
  "%{shared: -Bsymbolic}"

#define ANDROID_CC1_SPEC						\
  "%{!mglibc:%{!muclibc:%{!mbionic: -mbionic}}} "			\
  "%{!fno-pic:%{!fno-PIC:%{!fpic:%{!fPIC: -fPIC}}}}"

#define ANDROID_CC1PLUS_SPEC						\
  "%{!fexceptions:%{!fno-exceptions: -fno-exceptions}} "		\
  "%{!frtti:%{!fno-rtti: -fno-rtti}}"

#define ANDROID_LIB_SPEC \
  "%{!static: -ldl}"

#define ANDROID_STARTFILE_SPEC						\
  "%{shared: crtbegin_so%O%s;:"						\
  "  %{static: crtbegin_static%O%s;: crtbegin_dynamic%O%s}}"

#define ANDROID_ENDFILE_SPEC \
  "%{shared: crtend_so%O%s;: crtend_android%O%s}"
