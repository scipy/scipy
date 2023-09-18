/* Copyright (C) 1996, 1997, 1998, 1999, 2004, 2008
   Free Software Foundation, Inc.
   This file is part of the GNU C Library.
   Contributed by Ulrich Drepper <drepper@cygnus.com>, 1996.

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

#ifndef _REGEXP_H
#define _REGEXP_H	1

/* The contents of this header file was first standardized in X/Open
   System Interface and Headers Issue 2, originally coming from SysV.
   In issue 4, version 2, it is marked as TO BE WITDRAWN, and it has
   been withdrawn in SUSv3.

   This code shouldn't be used in any newly written code.  It is
   included only for compatibility reasons.  Use the POSIX definition
   in <regex.h> for portable applications and a reasonable interface.  */

#include <features.h>
#include <alloca.h>
#include <regex.h>
#include <stdlib.h>
#include <string.h>

/* The implementation provided here emulates the needed functionality
   by mapping to the POSIX regular expression matcher.  The interface
   for the here included function is weird (this really is a harmless
   word).

   The user has to provide six macros before this header file can be
   included:

   INIT		Declarations vor variables which can be used by the
		other macros.

   GETC()	Return the value of the next character in the regular
		expression pattern.  Successive calls should return
		successive characters.

   PEEKC()	Return the value of the next character in the regular
		expression pattern.  Immediately successive calls to
		PEEKC() should return the same character which should
		also be the next character returned by GETC().

   UNGETC(c)	Cause `c' to be returned by the next call to GETC() and
		PEEKC().

   RETURN(ptr)	Used for normal exit of the `compile' function.  `ptr'
		is a pointer to the character after the last character of
		the compiled regular expression.

   ERROR(val)	Used for abnormal return from `compile'.  `val' is the
		error number.  The error codes are:
		11	Range endpoint too large.
		16	Bad number.
		25	\digit out of range.
		36	Illegal or missing delimiter.
		41	No remembered search string.
		42	\( \) imbalance.
		43	Too many \(.
		44	More tan two numbers given in \{ \}.
		45	} expected after \.
		46	First number exceeds second in \{ \}.
		49	[ ] imbalance.
		50	Regular expression overflow.

  */

__BEGIN_DECLS

/* Interface variables.  They contain the results of the successful
   calls to `setp' and `advance'.  */
extern char *loc1;
extern char *loc2;

/* The use of this variable in the `advance' function is not
   supported.  */
extern char *locs;


#ifndef __DO_NOT_DEFINE_COMPILE
/* Get and compile the user supplied pattern up to end of line or
   string or until EOF is seen, whatever happens first.  The result is
   placed in the buffer starting at EXPBUF and delimited by ENDBUF.

   This function cannot be defined in the libc itself since it depends
   on the macros.  */
char *
compile (char *__restrict instring, char *__restrict expbuf,
	 __const char *__restrict endbuf, int eof)
{
  char *__input_buffer = NULL;
  size_t __input_size = 0;
  size_t __current_size = 0;
  int __ch;
  int __error;
  INIT

  /* Align the expression buffer according to the needs for an object
     of type `regex_t'.  Then check for minimum size of the buffer for
     the compiled regular expression.  */
  regex_t *__expr_ptr;
# if defined __GNUC__ && __GNUC__ >= 2
  const size_t __req = __alignof__ (regex_t *);
# else
  /* How shall we find out?  We simply guess it and can change it is
     this really proofs to be wrong.  */
  const size_t __req = 8;
# endif
  expbuf += __req;
  expbuf -= (expbuf - ((char *) 0)) % __req;
  if (endbuf < expbuf + sizeof (regex_t))
    {
      ERROR (50);
    }
  __expr_ptr = (regex_t *) expbuf;
  /* The remaining space in the buffer can be used for the compiled
     pattern.  */
  __expr_ptr->__REPB_PREFIX (buffer) = expbuf + sizeof (regex_t);
  __expr_ptr->__REPB_PREFIX (allocated)
    = endbuf - (char *) __expr_ptr->__REPB_PREFIX (buffer);

  while ((__ch = (GETC ())) != eof)
    {
      if (__ch == '\0' || __ch == '\n')
	{
	  UNGETC (__ch);
	  break;
	}

      if (__current_size + 1 >= __input_size)
	{
	  size_t __new_size = __input_size ? 2 * __input_size : 128;
	  char *__new_room = (char *) alloca (__new_size);
	  /* See whether we can use the old buffer.  */
	  if (__new_room + __new_size == __input_buffer)
	    {
	      __input_size += __new_size;
	      __input_buffer = (char *) memcpy (__new_room, __input_buffer,
					       __current_size);
	    }
	  else if (__input_buffer + __input_size == __new_room)
	    __input_size += __new_size;
	  else
	    {
	      __input_size = __new_size;
	      __input_buffer = (char *) memcpy (__new_room, __input_buffer,
						__current_size);
	    }
	}
      __input_buffer[__current_size++] = __ch;
    }
  if (__current_size)
    __input_buffer[__current_size++] = '\0';
  else
    __input_buffer = "";

  /* Now compile the pattern.  */
  __error = regcomp (__expr_ptr, __input_buffer, REG_NEWLINE);
  if (__error != 0)
    /* Oh well, we have to translate POSIX error codes.  */
    switch (__error)
      {
      case REG_BADPAT:
      case REG_ECOLLATE:
      case REG_ECTYPE:
      case REG_EESCAPE:
      case REG_BADRPT:
      case REG_EEND:
      case REG_ERPAREN:
      default:
	/* There is no matching error code.  */
	RETURN (36);
      case REG_ESUBREG:
	RETURN (25);
      case REG_EBRACK:
	RETURN (49);
      case REG_EPAREN:
	RETURN (42);
      case REG_EBRACE:
	RETURN (44);
      case REG_BADBR:
	RETURN (46);
      case REG_ERANGE:
	RETURN (11);
      case REG_ESPACE:
      case REG_ESIZE:
	ERROR (50);
      }

  /* Everything is ok.  */
  RETURN ((char *) (__expr_ptr->__REPB_PREFIX (buffer)
		    + __expr_ptr->__REPB_PREFIX (used)));
}
#endif


/* Find the next match in STRING.  The compiled regular expression is
   found in the buffer starting at EXPBUF.  `loc1' will return the
   first character matched and `loc2' points to the next unmatched
   character.  */
extern int step (__const char *__restrict __string,
		 __const char *__restrict __expbuf) __THROW;

/* Match the beginning of STRING with the compiled regular expression
   in EXPBUF.  If the match is successful `loc2' will contain the
   position of the first unmatched character.  */
extern int advance (__const char *__restrict __string,
		    __const char *__restrict __expbuf) __THROW;


__END_DECLS

#endif /* regexp.h */
