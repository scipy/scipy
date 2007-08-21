/*
 * pstdio.h -- $Id$
 * portability layer I/O wrappers
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

/* filesystem services (mostly ANSI or POSIX)
 * - necessary to implement UNIX-like filenaming semantics universally
 * - p_free frees result of p_getcwd, p_frall for p_lsdir result */

#ifndef PSTDIO_H
#define PSTDIO_H

#include "extern_c.h"

typedef struct p_file p_file;
typedef struct p_dir p_dir;

extern p_file *p_fopen(const char *unix_name, const char *mode);
extern p_file *p_popen(const char *command, const char *mode);

extern unsigned long p_fsize(p_file *file);
extern unsigned long p_ftell(p_file *file);
extern int p_fseek(p_file *file, unsigned long addr);

extern char *p_fgets(p_file *file, char *buf, int buflen);
extern int p_fputs(p_file *file, const char *buf);
extern unsigned long p_fread(p_file *file,
                             void *buf, unsigned long nbytes);
extern unsigned long p_fwrite(p_file *file,
                              const void *buf, unsigned long nbytes);

extern int p_feof(p_file *file);
extern int p_ferror(p_file *file);
extern int p_fflush(p_file *file);
extern int p_fclose(p_file *file);

extern int p_remove(const char *unix_name);
extern int p_rename(const char *unix_old, const char *unix_new);

extern int p_chdir(const char *unix_name);
extern int p_rmdir(const char *unix_name);
extern int p_mkdir(const char *unix_name);
extern char *p_getcwd(void);

extern p_dir *p_dopen(const char *unix_name);
extern int p_dclose(p_dir *dir);
/* returned filename does not need to be freed, but
 * value may be clobbered by dclose, next dnext, or p_wkspc use
 * . and .. do not appear in returned list */
extern char *p_dnext(p_dir *dir, int *is_dir);


END_EXTERN_C

#endif
