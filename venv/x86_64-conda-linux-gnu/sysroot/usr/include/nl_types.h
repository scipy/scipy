/* Copyright (C) 1996, 1997, 1999, 2003, 2004 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */

#ifndef _NL_TYPES_H
#define _NL_TYPES_H 1

#include <features.h>

/* The default message set used by the gencat program.  */
#define NL_SETD 1

/* Value for FLAG parameter of `catgets' to say we want XPG4 compliance.  */
#define NL_CAT_LOCALE 1


__BEGIN_DECLS

/* Message catalog descriptor type.  */
typedef void *nl_catd;

/* Type used by `nl_langinfo'.  */
typedef int nl_item;

/* Open message catalog for later use, returning descriptor.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern nl_catd catopen (__const char *__cat_name, int __flag) __nonnull ((1));

/* Return translation with NUMBER in SET of CATALOG; if not found
   return STRING.  */
extern char *catgets (nl_catd __catalog, int __set, int __number,
		      __const char *__string) __THROW __nonnull ((1));

/* Close message CATALOG.  */
extern int catclose (nl_catd __catalog) __THROW __nonnull ((1));

__END_DECLS

#endif /* nl_types.h  */
