/*
 * dir.c -- $Id$
 * UNIX version of play directory operations
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#ifndef _POSIX_SOURCE
#define _POSIX_SOURCE 1
#endif
#ifndef _XOPEN_SOURCE
/* Digital UNIX needs _XOPEN_SOURCE to define S_IFDIR bit */
#define _XOPEN_SOURCE 1
#endif

#include "config.h"

#include "play.h"
#include "playu.h"
#include "pstdlib.h"
#include "pstdio.h"

#include <stdio.h>
#include <string.h>
/* chdir, rmdir in unistd.h; mkdir, stat in sys/stat.h */
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

/* DIRENT_HEADER might be: <sys/dir.h>, <sys/ndir.h>, <ndir.h>
 * (see autoconf AC_HEADER_DIRENT) */
#ifndef DIRENT_HEADER
# include <dirent.h>
#else
# define dirent direct
# include DIRENT_HEADER
#endif
#ifndef S_ISDIR
 #ifndef S_IFDIR
  #define S_IFDIR _S_IFDIR
  # ifndef _S_IFDIR
  #define _S_IFDIR _IFDIR
  # endif
 #endif
 #define S_ISDIR(m) ((m)&S_IFDIR)
#endif

#ifdef TEST_DIRENT
#define p_strcpy
#define u_pathname
#define p_malloc malloc
#define p_free free
#define getcwd(x,y) (x)
#endif

#ifdef USE_GETWD
#define getcwd(x,y) getwd(x)
#endif

struct p_dir {
  DIR *dir;
  char *dirname;
  int namelen;
};

p_dir *
p_dopen(const char *unix_name)
{
  const char *name = u_pathname(unix_name);
  DIR *dir = opendir(name);
  p_dir *pdir = 0;
  if (dir) {
    pdir = p_malloc(sizeof(p_dir));
    if (pdir) {
      int len = strlen(name);
      pdir->dir = dir;
      pdir->dirname = p_malloc(len+2);
      strcpy(pdir->dirname, name);
      if (len>0 && name[len-1]!='/')
        pdir->dirname[len]='/', pdir->dirname[++len]='\0';
      pdir->namelen = len;
    }
  }
  return pdir;
}

int
p_dclose(p_dir *dir)
{
  int flag = closedir(dir->dir);
  p_free(dir->dirname);
  p_free(dir);
  return flag;
}

char *
p_dnext(p_dir *dir, int *is_dir)
{
  struct dirent *entry;
  char *name;
  do {
    entry = readdir(dir->dir);
    name = entry? entry->d_name : 0;
  } while (name && name[0]=='.' && ((name[1]=='.' && !name[2]) || !name[1]));
  if (name) {
    struct stat buf;
    strcpy(p_wkspc.c, dir->dirname);
    strncat(p_wkspc.c+dir->namelen, name, P_WKSIZ-dir->namelen);
    *is_dir = (!stat(p_wkspc.c,&buf) && S_ISDIR(buf.st_mode));
  }
  return name;
}

int
p_chdir(const char *dirname)
{
  return chdir(u_pathname(dirname));
}

int
p_rmdir(const char *dirname)
{
  return rmdir(u_pathname(dirname));
}

int
p_mkdir(const char *dirname)
{
  /* if process umask is 022, directory mode will drop to 755 */
  return mkdir(u_pathname(dirname), 0777);
}

char *
p_getcwd(void)
{
  char *dir = getcwd(p_wkspc.c, P_WKSIZ);
  struct stat buf;
  if (dir && !strncmp(dir, "/tmp_mnt/", 9) && !stat(dir, &buf)) {
    /* demangle automounter names */
    dev_t device = buf.st_dev;
    ino_t inode = buf.st_ino;
    if (!stat(dir+8, &buf) &&
        buf.st_dev==device && buf.st_ino==inode) {
      dir += 8;
    } else {
      int i;
      for (i=9 ; dir[i] ; i++) if (dir[i]=='/') break;
      if (dir[i]=='/' && !stat(dir+i, &buf) &&
          buf.st_dev==device && buf.st_ino==inode)
        dir += i;
    }
  }
  return dir;
}
