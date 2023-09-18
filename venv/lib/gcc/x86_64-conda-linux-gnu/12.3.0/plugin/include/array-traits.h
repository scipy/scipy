/* Descriptions of array-like objects.
   Copyright (C) 2019-2022 Free Software Foundation, Inc.

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

#ifndef GCC_ARRAY_TRAITS_H
#define GCC_ARRAY_TRAITS_H

/* Implementation for single integers (and similar types).  */
template<typename T, T zero = T (0)>
struct scalar_array_traits
{
  typedef T element_type;
  static const bool has_constant_size = true;
  static const size_t constant_size = 1;
  static const T *base (const T &x) { return &x; }
  static size_t size (const T &) { return 1; }
};

template<typename T>
struct array_traits : scalar_array_traits<T> {};

/* Implementation for arrays with a static size.  */
template<typename T, size_t N>
struct array_traits<T[N]>
{
  typedef T element_type;
  static const bool has_constant_size = true;
  static const size_t constant_size = N;
  static const T *base (const T (&x)[N]) { return x; }
  static size_t size (const T (&)[N]) { return N; }
};

#endif
