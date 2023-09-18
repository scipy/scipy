/* Copyright (C) 1991-2003,2006,2009 Free Software Foundation, Inc.
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

/*
 *	ISO C99 Standard: 7.23 Date and time	<time.h>
 */

#ifndef	_TIME_H

#if (! defined __need_time_t && !defined __need_clock_t && \
     ! defined __need_timespec)
# define _TIME_H	1
# include <features.h>

__BEGIN_DECLS

#endif

#ifdef	_TIME_H
/* Get size_t and NULL from <stddef.h>.  */
# define __need_size_t
# define __need_NULL
# include <stddef.h>

/* This defines CLOCKS_PER_SEC, which is the number of processor clock
   ticks per second.  */
# include <bits/time.h>

/* This is the obsolete POSIX.1-1988 name for the same constant.  */
# if !defined __STRICT_ANSI__ && !defined __USE_XOPEN2K
#  ifndef CLK_TCK
#   define CLK_TCK	CLOCKS_PER_SEC
#  endif
# endif

#endif /* <time.h> included.  */

#if !defined __clock_t_defined && (defined _TIME_H || defined __need_clock_t)
# define __clock_t_defined	1

# include <bits/types.h>

__BEGIN_NAMESPACE_STD
/* Returned by `clock'.  */
typedef __clock_t clock_t;
__END_NAMESPACE_STD
#if defined __USE_XOPEN || defined __USE_POSIX || defined __USE_MISC
__USING_NAMESPACE_STD(clock_t)
#endif

#endif /* clock_t not defined and <time.h> or need clock_t.  */
#undef	__need_clock_t

#if !defined __time_t_defined && (defined _TIME_H || defined __need_time_t)
# define __time_t_defined	1

# include <bits/types.h>

__BEGIN_NAMESPACE_STD
/* Returned by `time'.  */
typedef __time_t time_t;
__END_NAMESPACE_STD
#if defined __USE_POSIX || defined __USE_MISC || defined __USE_SVID
__USING_NAMESPACE_STD(time_t)
#endif

#endif /* time_t not defined and <time.h> or need time_t.  */
#undef	__need_time_t

#if !defined __clockid_t_defined && \
   ((defined _TIME_H && defined __USE_POSIX199309) || defined __need_clockid_t)
# define __clockid_t_defined	1

# include <bits/types.h>

/* Clock ID used in clock and timer functions.  */
typedef __clockid_t clockid_t;

#endif /* clockid_t not defined and <time.h> or need clockid_t.  */
#undef	__clockid_time_t

#if !defined __timer_t_defined && \
    ((defined _TIME_H && defined __USE_POSIX199309) || defined __need_timer_t)
# define __timer_t_defined	1

# include <bits/types.h>

/* Timer ID returned by `timer_create'.  */
typedef __timer_t timer_t;

#endif /* timer_t not defined and <time.h> or need timer_t.  */
#undef	__need_timer_t


#if !defined __timespec_defined &&				\
    ((defined _TIME_H &&					\
      (defined __USE_POSIX199309 || defined __USE_MISC)) ||	\
      defined __need_timespec)
# define __timespec_defined	1

# include <bits/types.h>	/* This defines __time_t for us.  */

/* POSIX.1b structure for a time value.  This is like a `struct timeval' but
   has nanoseconds instead of microseconds.  */
struct timespec
  {
    __time_t tv_sec;		/* Seconds.  */
    long int tv_nsec;		/* Nanoseconds.  */
  };

#endif /* timespec not defined and <time.h> or need timespec.  */
#undef	__need_timespec


#ifdef	_TIME_H
__BEGIN_NAMESPACE_STD
/* Used by other time functions.  */
struct tm
{
  int tm_sec;			/* Seconds.	[0-60] (1 leap second) */
  int tm_min;			/* Minutes.	[0-59] */
  int tm_hour;			/* Hours.	[0-23] */
  int tm_mday;			/* Day.		[1-31] */
  int tm_mon;			/* Month.	[0-11] */
  int tm_year;			/* Year	- 1900.  */
  int tm_wday;			/* Day of week.	[0-6] */
  int tm_yday;			/* Days in year.[0-365]	*/
  int tm_isdst;			/* DST.		[-1/0/1]*/

#ifdef	__USE_BSD
  long int tm_gmtoff;		/* Seconds east of UTC.  */
  __const char *tm_zone;	/* Timezone abbreviation.  */
#else
  long int __tm_gmtoff;		/* Seconds east of UTC.  */
  __const char *__tm_zone;	/* Timezone abbreviation.  */
#endif
};
__END_NAMESPACE_STD
#if defined __USE_XOPEN || defined __USE_POSIX || defined __USE_MISC
__USING_NAMESPACE_STD(tm)
#endif


#ifdef __USE_POSIX199309
/* POSIX.1b structure for timer start values and intervals.  */
struct itimerspec
  {
    struct timespec it_interval;
    struct timespec it_value;
  };

/* We can use a simple forward declaration.  */
struct sigevent;

#endif	/* POSIX.1b */

#ifdef __USE_XOPEN2K
# ifndef __pid_t_defined
typedef __pid_t pid_t;
#  define __pid_t_defined
# endif
#endif


__BEGIN_NAMESPACE_STD
/* Time used by the program so far (user time + system time).
   The result / CLOCKS_PER_SECOND is program time in seconds.  */
extern clock_t clock (void) __THROW;

/* Return the current time and put it in *TIMER if TIMER is not NULL.  */
extern time_t time (time_t *__timer) __THROW;

/* Return the difference between TIME1 and TIME0.  */
extern double difftime (time_t __time1, time_t __time0)
     __THROW __attribute__ ((__const__));

/* Return the `time_t' representation of TP and normalize TP.  */
extern time_t mktime (struct tm *__tp) __THROW;


/* Format TP into S according to FORMAT.
   Write no more than MAXSIZE characters and return the number
   of characters written, or 0 if it would exceed MAXSIZE.  */
extern size_t strftime (char *__restrict __s, size_t __maxsize,
			__const char *__restrict __format,
			__const struct tm *__restrict __tp) __THROW;
__END_NAMESPACE_STD

# ifdef __USE_XOPEN
/* Parse S according to FORMAT and store binary time information in TP.
   The return value is a pointer to the first unparsed character in S.  */
extern char *strptime (__const char *__restrict __s,
		       __const char *__restrict __fmt, struct tm *__tp)
     __THROW;
# endif

# ifdef __USE_XOPEN2K8
/* Similar to the two functions above but take the information from
   the provided locale and not the global locale.  */
# include <xlocale.h>

extern size_t strftime_l (char *__restrict __s, size_t __maxsize,
			  __const char *__restrict __format,
			  __const struct tm *__restrict __tp,
			  __locale_t __loc) __THROW;
# endif

# ifdef __USE_GNU
extern char *strptime_l (__const char *__restrict __s,
			 __const char *__restrict __fmt, struct tm *__tp,
			 __locale_t __loc) __THROW;
# endif


__BEGIN_NAMESPACE_STD
/* Return the `struct tm' representation of *TIMER
   in Universal Coordinated Time (aka Greenwich Mean Time).  */
extern struct tm *gmtime (__const time_t *__timer) __THROW;

/* Return the `struct tm' representation
   of *TIMER in the local timezone.  */
extern struct tm *localtime (__const time_t *__timer) __THROW;
__END_NAMESPACE_STD

# if defined __USE_POSIX || defined __USE_MISC
/* Return the `struct tm' representation of *TIMER in UTC,
   using *TP to store the result.  */
extern struct tm *gmtime_r (__const time_t *__restrict __timer,
			    struct tm *__restrict __tp) __THROW;

/* Return the `struct tm' representation of *TIMER in local time,
   using *TP to store the result.  */
extern struct tm *localtime_r (__const time_t *__restrict __timer,
			       struct tm *__restrict __tp) __THROW;
# endif	/* POSIX or misc */

__BEGIN_NAMESPACE_STD
/* Return a string of the form "Day Mon dd hh:mm:ss yyyy\n"
   that is the representation of TP in this format.  */
extern char *asctime (__const struct tm *__tp) __THROW;

/* Equivalent to `asctime (localtime (timer))'.  */
extern char *ctime (__const time_t *__timer) __THROW;
__END_NAMESPACE_STD

# if defined __USE_POSIX || defined __USE_MISC
/* Reentrant versions of the above functions.  */

/* Return in BUF a string of the form "Day Mon dd hh:mm:ss yyyy\n"
   that is the representation of TP in this format.  */
extern char *asctime_r (__const struct tm *__restrict __tp,
			char *__restrict __buf) __THROW;

/* Equivalent to `asctime_r (localtime_r (timer, *TMP*), buf)'.  */
extern char *ctime_r (__const time_t *__restrict __timer,
		      char *__restrict __buf) __THROW;
# endif	/* POSIX or misc */


/* Defined in localtime.c.  */
extern char *__tzname[2];	/* Current timezone names.  */
extern int __daylight;		/* If daylight-saving time is ever in use.  */
extern long int __timezone;	/* Seconds west of UTC.  */


# ifdef	__USE_POSIX
/* Same as above.  */
extern char *tzname[2];

/* Set time conversion information from the TZ environment variable.
   If TZ is not defined, a locale-dependent default is used.  */
extern void tzset (void) __THROW;
# endif

# if defined __USE_SVID || defined __USE_XOPEN
extern int daylight;
extern long int timezone;
# endif

# ifdef __USE_SVID
/* Set the system time to *WHEN.
   This call is restricted to the superuser.  */
extern int stime (__const time_t *__when) __THROW;
# endif


/* Nonzero if YEAR is a leap year (every 4 years,
   except every 100th isn't, and every 400th is).  */
# define __isleap(year)	\
  ((year) % 4 == 0 && ((year) % 100 != 0 || (year) % 400 == 0))


# ifdef __USE_MISC
/* Miscellaneous functions many Unices inherited from the public domain
   localtime package.  These are included only for compatibility.  */

/* Like `mktime', but for TP represents Universal Time, not local time.  */
extern time_t timegm (struct tm *__tp) __THROW;

/* Another name for `mktime'.  */
extern time_t timelocal (struct tm *__tp) __THROW;

/* Return the number of days in YEAR.  */
extern int dysize (int __year) __THROW  __attribute__ ((__const__));
# endif


# ifdef __USE_POSIX199309
/* Pause execution for a number of nanoseconds.

   This function is a cancellation point and therefore not marked with
   __THROW.  */
extern int nanosleep (__const struct timespec *__requested_time,
		      struct timespec *__remaining);


/* Get resolution of clock CLOCK_ID.  */
extern int clock_getres (clockid_t __clock_id, struct timespec *__res) __THROW;

/* Get current value of clock CLOCK_ID and store it in TP.  */
extern int clock_gettime (clockid_t __clock_id, struct timespec *__tp) __THROW;

/* Set clock CLOCK_ID to value TP.  */
extern int clock_settime (clockid_t __clock_id, __const struct timespec *__tp)
     __THROW;

#  ifdef __USE_XOPEN2K
/* High-resolution sleep with the specified clock.

   This function is a cancellation point and therefore not marked with
   __THROW.  */
extern int clock_nanosleep (clockid_t __clock_id, int __flags,
			    __const struct timespec *__req,
			    struct timespec *__rem);

/* Return clock ID for CPU-time clock.  */
extern int clock_getcpuclockid (pid_t __pid, clockid_t *__clock_id) __THROW;
#  endif


/* Create new per-process timer using CLOCK_ID.  */
extern int timer_create (clockid_t __clock_id,
			 struct sigevent *__restrict __evp,
			 timer_t *__restrict __timerid) __THROW;

/* Delete timer TIMERID.  */
extern int timer_delete (timer_t __timerid) __THROW;

/* Set timer TIMERID to VALUE, returning old value in OVLAUE.  */
extern int timer_settime (timer_t __timerid, int __flags,
			  __const struct itimerspec *__restrict __value,
			  struct itimerspec *__restrict __ovalue) __THROW;

/* Get current value of timer TIMERID and store it in VLAUE.  */
extern int timer_gettime (timer_t __timerid, struct itimerspec *__value)
     __THROW;

/* Get expiration overrun for timer TIMERID.  */
extern int timer_getoverrun (timer_t __timerid) __THROW;
# endif


# ifdef __USE_XOPEN_EXTENDED
/* Set to one of the following values to indicate an error.
     1  the DATEMSK environment variable is null or undefined,
     2  the template file cannot be opened for reading,
     3  failed to get file status information,
     4  the template file is not a regular file,
     5  an error is encountered while reading the template file,
     6  memory allication failed (not enough memory available),
     7  there is no line in the template that matches the input,
     8  invalid input specification Example: February 31 or a time is
        specified that can not be represented in a time_t (representing
	the time in seconds since 00:00:00 UTC, January 1, 1970) */
extern int getdate_err;

/* Parse the given string as a date specification and return a value
   representing the value.  The templates from the file identified by
   the environment variable DATEMSK are used.  In case of an error
   `getdate_err' is set.

   This function is a possible cancellation points and therefore not
   marked with __THROW.  */
extern struct tm *getdate (__const char *__string);
# endif

# ifdef __USE_GNU
/* Since `getdate' is not reentrant because of the use of `getdate_err'
   and the static buffer to return the result in, we provide a thread-safe
   variant.  The functionality is the same.  The result is returned in
   the buffer pointed to by RESBUFP and in case of an error the return
   value is != 0 with the same values as given above for `getdate_err'.

   This function is not part of POSIX and therefore no official
   cancellation point.  But due to similarity with an POSIX interface
   or due to the implementation it is a cancellation point and
   therefore not marked with __THROW.  */
extern int getdate_r (__const char *__restrict __string,
		      struct tm *__restrict __resbufp);
# endif

__END_DECLS

#endif /* <time.h> included.  */

#endif /* <time.h> not already included.  */
