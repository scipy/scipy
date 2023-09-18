/* Copyright (C) 2007 Free Software Foundation, Inc.
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

#ifndef __NETIUCV_IUCV_H
#define __NETIUCV_IUCV_H	1

#include <features.h>
#include <bits/sockaddr.h>

__BEGIN_DECLS

struct sockaddr_iucv
  {
    __SOCKADDR_COMMON (siucv_);
    unsigned short	siucv_port;		/* Reserved */
    unsigned int	siucv_addr;		/* Reserved */
    char		siucv_nodeid[8];	/* Reserved */
    char		siucv_user_id[8];	/* Guest User Id */
    char		siucv_name[8];		/* Application Name */
  };

__END_DECLS

#endif
