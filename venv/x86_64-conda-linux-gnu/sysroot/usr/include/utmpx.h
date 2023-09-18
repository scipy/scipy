/* Copyright (C) 1997, 1998, 1999, 2003 Free Software Foundation, Inc.
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

#ifndef	_UTMPX_H
#define	_UTMPX_H	1

#include <features.h>
#include <sys/time.h>

/* Required according to Unix98.  */
#ifndef __pid_t_defined
typedef __pid_t pid_t;
# define __pid_t_defined
#endif

/* Get system dependent values and data structures.  */
#include <bits/utmpx.h>

#ifdef __USE_GNU
/* Compatibility names for the strings of the canonical file names.  */
# define UTMPX_FILE	_PATH_UTMPX
# define UTMPX_FILENAME	_PATH_UTMPX
# define WTMPX_FILE	_PATH_WTMPX
# define WTMPX_FILENAME	_PATH_WTMPX
#endif

/* For the getutmp{,x} functions we need the `struct utmp'.  */
#ifdef __USE_GNU
struct utmp;
#endif


__BEGIN_DECLS

/* Open user accounting database.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern void setutxent (void);

/* Close user accounting database.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern void endutxent (void);

/* Get the next entry from the user accounting database.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern struct utmpx *getutxent (void);

/* Get the user accounting database entry corresponding to ID.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern struct utmpx *getutxid (__const struct utmpx *__id);

/* Get the user accounting database entry corresponding to LINE.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern struct utmpx *getutxline (__const struct utmpx *__line);

/* Write the entry UTMPX into the user accounting database.

   This function is a possible cancellation point and therefore not
   marked with __THROW.  */
extern struct utmpx *pututxline (__const struct utmpx *__utmpx);


#ifdef __USE_GNU
/* Change name of the utmpx file to be examined.

   This function is not part of POSIX and therefore no official
   cancellation point.  But due to similarity with an POSIX interface
   or due to the implementation it is a cancellation point and
   therefore not marked with __THROW.  */
extern int utmpxname (__const char *__file);

/* Append entry UTMP to the wtmpx-like file WTMPX_FILE.

   This function is not part of POSIX and therefore no official
   cancellation point.  But due to similarity with an POSIX interface
   or due to the implementation it is a cancellation point and
   therefore not marked with __THROW.  */
extern void updwtmpx (__const char *__wtmpx_file,
		      __const struct utmpx *__utmpx);


/* Copy the information in UTMPX to UTMP.

   This function is not part of POSIX and therefore no official
   cancellation point.  But due to similarity with an POSIX interface
   or due to the implementation it is a cancellation point and
   therefore not marked with __THROW.  */
extern void getutmp (__const struct utmpx *__utmpx,
		     struct utmp *__utmp);

/* Copy the information in UTMP to UTMPX.

   This function is not part of POSIX and therefore no official
   cancellation point.  But due to similarity with an POSIX interface
   or due to the implementation it is a cancellation point and
   therefore not marked with __THROW.  */
extern void getutmpx (__const struct utmp *__utmp,
		      struct utmpx *__utmpx);
#endif

__END_DECLS

#endif /* utmpx.h  */
