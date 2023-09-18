/* Warn on violations of the restrict qualifier.
   Copyright (C) 2017-2022 Free Software Foundation, Inc.
   Contributed by Martin Sebor <msebor@redhat.com>.

   This file is part of GCC.

   GCC is free software; you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free
   Software Foundation; either version 3, or (at your option) any later
   version.

   GCC is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   for more details.

   You should have received a copy of the GNU General Public License
   along with GCC; see the file COPYING3.  If not see
   <http://www.gnu.org/licenses/>.  */

#ifndef GIMPLE_SSA_WARN_RESTRICT_H

extern opt_code check_bounds_or_overlap (gimple *, tree, tree, tree, tree,
					 bool = false, bool = true);
extern opt_code check_bounds_or_overlap (class pointer_query &, gimple *,
					 tree, tree, tree, tree,
					 bool = false, bool = true);

#endif /* GIMPLE_SSA_WARN_RESTRICT_H */
