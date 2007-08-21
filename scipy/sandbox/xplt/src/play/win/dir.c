/*
 * dir.c -- $Id$
 * MS Windows version of plib directory operations
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

#include "config.h"
#include "pstdlib.h"
#include "pstdio.h"
#include "playw.h"

#include <string.h>

struct p_dir {
  HANDLE hd;
  int count;
  WIN32_FIND_DATA data;
};

p_dir *
p_dopen(const char *unix_name)
{
  char *name = w_pathname(unix_name);
  p_dir *dir = p_malloc(sizeof(p_dir));
  int i = 0;
  while (name[i]) i++;
  if (i>0 && name[i-1]!='\\') name[i++] = '\\';
  name[i++] = '*';
  name[i++] = '\0';
  dir->hd = FindFirstFile(name, &dir->data);
  if (dir->hd != INVALID_HANDLE_VALUE) {
    /* even empty directories contain . and .. */
    dir->count = 0;
  } else {
    p_free(dir);
    dir = 0;
  }
  return dir;
}

int
p_dclose(p_dir *dir)
{
  int flag = -(!FindClose(dir->hd));
  p_free(dir);
  return flag;
}

char *
p_dnext(p_dir *dir, int *is_dir)
{
  for (;;) {
    if (dir->count++) {
      if (!FindNextFile(dir->hd, &dir->data))
        return 0;   /* GetLastError()==ERROR_NO_MORE_FILES or error */
    }
    if (dir->data.cFileName[0]!='.' ||
        (dir->data.cFileName[1] && (dir->data.cFileName[1]!='.' ||
                                    dir->data.cFileName[2]))) break;
  }
  *is_dir = ((dir->data.dwFileAttributes&FILE_ATTRIBUTE_DIRECTORY) != 0);
  return (char *)dir->data.cFileName;
}

int
p_chdir(const char *dirname)
{
  const char *path = w_pathname(dirname);
  /* is this necessary?
  if (path[0] && path[1]==':' &&
      ((path[0]>='a' && path[0]<='z') || (path[0]>='A' && path[0]<='Z'))) {
    char drive[8];
    drive[0] = *path++;
    drive[1] = *path++;
    drive[2] = '\0';
    if (!SetCurrentDrive(drive)) return -1;
  }
  */
  return -(!SetCurrentDirectory(path));
}

int
p_rmdir(const char *dirname)
{
  return -(!RemoveDirectory(w_pathname(dirname)));
}

int
p_mkdir(const char *dirname)
{
  return -(!CreateDirectory(w_pathname(dirname), 0));
}

char *
p_getcwd(void)
{
  DWORD n = GetCurrentDirectory(P_WKSIZ, p_wkspc.c);
  if (n>P_WKSIZ || n==0) return 0;
  return w_unixpath(p_wkspc.c);
}
