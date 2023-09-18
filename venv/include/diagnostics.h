/* Copyright (C) 2017-2023 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H

/* If at all possible, fix the source rather than using these macros
   to silence warnings.  If you do use these macros be aware that
   you'll need to condition their use on particular compiler versions,
   which can be done for gcc using ansidecl.h's GCC_VERSION macro.

   gcc versions between 4.2 and 4.6 do not allow pragma control of
   diagnostics inside functions, giving a hard error if you try to use
   the finer control available with later versions.
   gcc prior to 4.2 warns about diagnostic push and pop.

   The other macros have restrictions too, for example gcc-5, gcc-6
   and gcc-7 warn that -Wstringop-truncation is unknown, unless you
   also add DIAGNOSTIC_IGNORE ("-Wpragma").  */

#ifdef __GNUC__
# define DIAGNOSTIC_PUSH _Pragma ("GCC diagnostic push")
# define DIAGNOSTIC_POP _Pragma ("GCC diagnostic pop")

/* Stringification.  */
# define DIAGNOSTIC_STRINGIFY_1(x) #x
# define DIAGNOSTIC_STRINGIFY(x) DIAGNOSTIC_STRINGIFY_1 (x)

# define DIAGNOSTIC_IGNORE(option) \
  _Pragma (DIAGNOSTIC_STRINGIFY (GCC diagnostic ignored option))
# define DIAGNOSTIC_ERROR(option) \
  _Pragma (DIAGNOSTIC_STRINGIFY (GCC diagnostic error option))
#else
# define DIAGNOSTIC_PUSH
# define DIAGNOSTIC_POP
# define DIAGNOSTIC_IGNORE(option)
#endif

#if defined (__clang__) /* clang */

# define DIAGNOSTIC_IGNORE_SELF_MOVE DIAGNOSTIC_IGNORE ("-Wself-move")
# define DIAGNOSTIC_IGNORE_DEPRECATED_DECLARATIONS \
  DIAGNOSTIC_IGNORE ("-Wdeprecated-declarations")
# define DIAGNOSTIC_IGNORE_DEPRECATED_REGISTER \
  DIAGNOSTIC_IGNORE ("-Wdeprecated-register")
# if __has_warning ("-Wenum-compare-switch")
#  define DIAGNOSTIC_IGNORE_SWITCH_DIFFERENT_ENUM_TYPES \
   DIAGNOSTIC_IGNORE ("-Wenum-compare-switch")
# endif

# define DIAGNOSTIC_IGNORE_FORMAT_NONLITERAL \
  DIAGNOSTIC_IGNORE ("-Wformat-nonliteral")

# if __has_warning ("-Wuser-defined-warnings")
#  define DIAGNOSTIC_IGNORE_USER_DEFINED_WARNINGS \
   DIAGNOSTIC_IGNORE ("-Wuser-defined-warnings")
# endif

# if __has_warning ("-Wunused-but-set-variable")
#  define DIAGNOSTIC_IGNORE_UNUSED_BUT_SET_VARIABLE \
   DIAGNOSTIC_IGNORE ("-Wunused-but-set-variable")
# endif

# define DIAGNOSTIC_ERROR_SWITCH \
  DIAGNOSTIC_ERROR ("-Wswitch")

#elif defined (__GNUC__) /* GCC */

# define DIAGNOSTIC_IGNORE_DEPRECATED_DECLARATIONS \
  DIAGNOSTIC_IGNORE ("-Wdeprecated-declarations")

# if __GNUC__ >= 7
#  define DIAGNOSTIC_IGNORE_DEPRECATED_REGISTER \
   DIAGNOSTIC_IGNORE ("-Wregister")
# endif

# define DIAGNOSTIC_IGNORE_STRINGOP_TRUNCATION \
  DIAGNOSTIC_IGNORE ("-Wstringop-truncation")

# if __GNUC__ >= 11
# define DIAGNOSTIC_IGNORE_STRINGOP_OVERREAD \
  DIAGNOSTIC_IGNORE ("-Wstringop-overread")
#endif

# define DIAGNOSTIC_IGNORE_FORMAT_NONLITERAL \
  DIAGNOSTIC_IGNORE ("-Wformat-nonliteral")

# if __GNUC__ >= 5
#  define DIAGNOSTIC_IGNORE_UNUSED_BUT_SET_VARIABLE \
   DIAGNOSTIC_IGNORE ("-Wunused-but-set-variable")
# endif

# if __GNUC__ >= 13
#  define DIAGNOSTIC_IGNORE_SELF_MOVE DIAGNOSTIC_IGNORE ("-Wself-move")
# endif

/* GCC 4.8's "diagnostic push/pop" seems broken when using this, -Wswitch
   remains enabled at the error level even after a pop.  Therefore, don't
   use it for GCC < 5.  */
# if __GNUC__ >= 5
#  define DIAGNOSTIC_ERROR_SWITCH DIAGNOSTIC_ERROR ("-Wswitch")
# endif

#endif

#ifndef DIAGNOSTIC_IGNORE_SELF_MOVE
# define DIAGNOSTIC_IGNORE_SELF_MOVE
#endif

#ifndef DIAGNOSTIC_IGNORE_DEPRECATED_DECLARATIONS
# define DIAGNOSTIC_IGNORE_DEPRECATED_DECLARATIONS
#endif

#ifndef DIAGNOSTIC_IGNORE_DEPRECATED_REGISTER
# define DIAGNOSTIC_IGNORE_DEPRECATED_REGISTER
#endif

#ifndef DIAGNOSTIC_IGNORE_SWITCH_DIFFERENT_ENUM_TYPES
# define DIAGNOSTIC_IGNORE_SWITCH_DIFFERENT_ENUM_TYPES
#endif

#ifndef DIAGNOSTIC_IGNORE_STRINGOP_TRUNCATION
# define DIAGNOSTIC_IGNORE_STRINGOP_TRUNCATION
#endif

#ifndef DIAGNOSTIC_IGNORE_STRINGOP_OVERREAD
# define DIAGNOSTIC_IGNORE_STRINGOP_OVERREAD
#endif

#ifndef DIAGNOSTIC_IGNORE_FORMAT_NONLITERAL
# define DIAGNOSTIC_IGNORE_FORMAT_NONLITERAL
#endif

#ifndef DIAGNOSTIC_IGNORE_USER_DEFINED_WARNINGS
# define DIAGNOSTIC_IGNORE_USER_DEFINED_WARNINGS
#endif

#ifndef DIAGNOSTIC_IGNORE_UNUSED_BUT_SET_VARIABLE
# define DIAGNOSTIC_IGNORE_UNUSED_BUT_SET_VARIABLE
#endif

#ifndef DIAGNOSTIC_ERROR_SWITCH
# define DIAGNOSTIC_ERROR_SWITCH
#endif

#endif /* DIAGNOSTICS_H */
