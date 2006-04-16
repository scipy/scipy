/*
 * config.c -- $Id$
 * configuration tester for UNIX machines
 *
 * Copyright (c) 1999.  See accompanying LEGAL file for details.
 */

#ifdef TEST_GCC
# ifdef __GNUC__
#define MAIN_BODY \
  value = 0;
# else
#error not gcc
# endif
#endif

#ifdef TEST_UTIME
/* check settings of: USE_GETRUSAGE USE_TIMES */
#include "timeu.c"
#define MAIN_BODY \
  double s; \
  double t = p_cpu_secs(&s); \
  value = 0;
#endif

#ifdef TEST_WTIME
/* check settings of: USE_GETTIMEOFDAY */
#include "timew.c"
#define MAIN_BODY \
  double t = p_wall_secs(); \
  value = 0;
#endif

#ifdef TEST_USERNM
/* check settings of: NO_PASSWD */
#include "usernm.c"
#define MAIN_BODY \
  char *u = p_getuser(); \
  value = (u!=0);
#endif

#ifdef TEST_TIOCGPGRP
/* check settings of: USE_TIOCGPGRP_IOCTL */
#include "uinbg.c"
#define MAIN_BODY \
  value = u_in_background();
#endif

#ifdef TEST_GETCWD
/* check settings of: USE_GETWD */
#include <unistd.h>
static char dirbuf[1024];
#ifdef USE_GETWD
#define getcwd(x,y) getwd(x)
#endif
#define MAIN_BODY \
  char *u = getcwd(dirbuf, 1024); \
  value = (u!=0);
#endif

#ifdef TEST_DIRENT
/* check settings of: DIRENT_HEADER USE_GETWD */
#include "dir.c"
p_twkspc p_wkspc;
#define MAIN_BODY \
  p_dir *d = p_dopen("no/such/thing"); \
  char *l = p_dnext(d, &value); \
  value = p_chdir(l) || p_rmdir(l) || p_mkdir(l);
#endif

#ifdef TEST_POLL
/* check settings of: USE_SYS_POLL_H USE_SELECT HAVE_SYS_SELECT_H
                      NO_SYS_TIME_H NEED_SELECT_PROTO */
#include "uevent.c"
#define MAIN_BODY \
  int p = u_poll(1000); \
  value = 0;
#endif

int
main(int argc, char *argv[])
{
  int value;
  MAIN_BODY
  return value;
}
