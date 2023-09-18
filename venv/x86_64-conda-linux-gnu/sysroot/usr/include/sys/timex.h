/* Copyright (C) 1995-1997, 1999, 2007, 2009 Free Software Foundation, Inc.
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

#ifndef	_SYS_TIMEX_H
#define	_SYS_TIMEX_H	1

#include <features.h>
#include <sys/time.h>

/* These definitions from linux/timex.h as of 2.6.30.  */

#define NTP_API	4	/* NTP API version */

struct ntptimeval
{
  struct timeval time;	/* current time (ro) */
  long int maxerror;	/* maximum error (us) (ro) */
  long int esterror;	/* estimated error (us) (ro) */
  long int tai;		/* TAI offset (ro) */

  long int __unused1;
  long int __unused2;
  long int __unused3;
  long int __unused4;
};

struct timex
{
  unsigned int modes;	/* mode selector */
  long int offset;	/* time offset (usec) */
  long int freq;	/* frequency offset (scaled ppm) */
  long int maxerror;	/* maximum error (usec) */
  long int esterror;	/* estimated error (usec) */
  int status;		/* clock command/status */
  long int constant;	/* pll time constant */
  long int precision;	/* clock precision (usec) (read only) */
  long int tolerance;	/* clock frequency tolerance (ppm) (read only) */
  struct timeval time;	/* (read only) */
  long int tick;	/* (modified) usecs between clock ticks */

  long int ppsfreq;	/* pps frequency (scaled ppm) (ro) */
  long int jitter;	/* pps jitter (us) (ro) */
  int shift;		/* interval duration (s) (shift) (ro) */
  long int stabil;	/* pps stability (scaled ppm) (ro) */
  long int jitcnt;	/* jitter limit exceeded (ro) */
  long int calcnt;	/* calibration intervals (ro) */
  long int errcnt;	/* calibration errors (ro) */
  long int stbcnt;	/* stability limit exceeded (ro) */

  int tai;		/* TAI offset (ro) */

  /* ??? */
  int  :32; int  :32; int  :32; int  :32;
  int  :32; int  :32; int  :32; int  :32;
  int  :32; int  :32; int  :32;
};

/* Mode codes (timex.mode) */
#define ADJ_OFFSET		0x0001	/* time offset */
#define ADJ_FREQUENCY		0x0002	/* frequency offset */
#define ADJ_MAXERROR		0x0004	/* maximum time error */
#define ADJ_ESTERROR		0x0008	/* estimated time error */
#define ADJ_STATUS		0x0010	/* clock status */
#define ADJ_TIMECONST		0x0020	/* pll time constant */
#define ADJ_TAI			0x0080	/* set TAI offset */
#define ADJ_MICRO		0x1000	/* select microsecond resolution */
#define ADJ_NANO		0x2000	/* select nanosecond resolution */
#define ADJ_TICK		0x4000	/* tick value */
#define ADJ_OFFSET_SINGLESHOT	0x8001	/* old-fashioned adjtime */
#define ADJ_OFFSET_SS_READ	0xa001	/* read-only adjtime */

/* xntp 3.4 compatibility names */
#define MOD_OFFSET	ADJ_OFFSET
#define MOD_FREQUENCY	ADJ_FREQUENCY
#define MOD_MAXERROR	ADJ_MAXERROR
#define MOD_ESTERROR	ADJ_ESTERROR
#define MOD_STATUS	ADJ_STATUS
#define MOD_TIMECONST	ADJ_TIMECONST
#define MOD_CLKB	ADJ_TICK
#define MOD_CLKA	ADJ_OFFSET_SINGLESHOT /* 0x8000 in original */
#define MOD_TAI		ADJ_TAI
#define MOD_MICRO	ADJ_MICRO
#define MOD_NANO	ADJ_NANO


/* Status codes (timex.status) */
#define STA_PLL		0x0001	/* enable PLL updates (rw) */
#define STA_PPSFREQ	0x0002	/* enable PPS freq discipline (rw) */
#define STA_PPSTIME	0x0004	/* enable PPS time discipline (rw) */
#define STA_FLL		0x0008	/* select frequency-lock mode (rw) */

#define STA_INS		0x0010	/* insert leap (rw) */
#define STA_DEL		0x0020	/* delete leap (rw) */
#define STA_UNSYNC	0x0040	/* clock unsynchronized (rw) */
#define STA_FREQHOLD	0x0080	/* hold frequency (rw) */

#define STA_PPSSIGNAL	0x0100	/* PPS signal present (ro) */
#define STA_PPSJITTER	0x0200	/* PPS signal jitter exceeded (ro) */
#define STA_PPSWANDER	0x0400	/* PPS signal wander exceeded (ro) */
#define STA_PPSERROR	0x0800	/* PPS signal calibration error (ro) */

#define STA_CLOCKERR	0x1000	/* clock hardware fault (ro) */
#define STA_NANO	0x2000	/* resolution (0 = us, 1 = ns) (ro) */
#define STA_MODE	0x4000	/* mode (0 = PLL, 1 = FLL) (ro) */
#define STA_CLK		0x8000	/* clock source (0 = A, 1 = B) (ro) */

/* Read-only bits */
#define STA_RONLY (STA_PPSSIGNAL | STA_PPSJITTER | STA_PPSWANDER | \
    STA_PPSERROR | STA_CLOCKERR | STA_NANO | STA_MODE | STA_CLK)

/* Clock states (time_state) */
#define TIME_OK		0	/* clock synchronized, no leap second */
#define TIME_INS	1	/* insert leap second */
#define TIME_DEL	2	/* delete leap second */
#define TIME_OOP	3	/* leap second in progress */
#define TIME_WAIT	4	/* leap second has occurred */
#define TIME_ERROR	5	/* clock not synchronized */
#define TIME_BAD	TIME_ERROR /* bw compat */

/* Maximum time constant of the PLL.  */
#define MAXTC		6

__BEGIN_DECLS

extern int __adjtimex (struct timex *__ntx) __THROW;
extern int adjtimex (struct timex *__ntx) __THROW;

#ifdef __REDIRECT_NTH
extern int __REDIRECT_NTH (ntp_gettime, (struct ntptimeval *__ntv),
			   ntp_gettimex);
#else
extern int ntp_gettimex (struct ntptimeval *__ntv) __THROW;
# define ntp_gettime ntp_gettimex
#endif
extern int ntp_adjtime (struct timex *__tntx) __THROW;

__END_DECLS

#endif /* sys/timex.h */
