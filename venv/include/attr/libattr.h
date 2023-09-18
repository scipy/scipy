/*
  Copyright (C) 2009  Andreas Gruenbacher <agruen@suse.de>

  This program is free software: you can redistribute it and/or modify it
  under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 2.1 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __LIBATTR_H
#define __LIBATTR_H

#ifdef __cplusplus
extern "C" {
#endif

struct error_context;

extern int attr_copy_file (const char *, const char *,
			   int (*) (const char *, struct error_context *),
			   struct error_context *);
extern int attr_copy_fd (const char *, int, const char *, int,
			 int (*) (const char *, struct error_context *),
			 struct error_context *);

/* Keep this function for backwards compatibility. */
extern int attr_copy_check_permissions(const char *, struct error_context *);

#define ATTR_ACTION_SKIP	1
#define ATTR_ACTION_PERMISSIONS	2

extern int attr_copy_action(const char *, struct error_context *);

#ifdef __cplusplus
}
#endif

#endif
