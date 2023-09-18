/* Language-independent APIs to enable/disable per-location warnings.

   Copyright (C) 2021-2022 Free Software Foundation, Inc.
   Contributed by Martin Sebor <msebor@redhat.com>

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

#ifndef DIAGNOSTIC_SPEC_H_INCLUDED
#define DIAGNOSTIC_SPEC_H_INCLUDED

#include "hash-map.h"

/* A "bitset" of warning groups.  */

class nowarn_spec_t
{
public:
  enum
    {
     /* Middle end warnings about invalid accesses.  */
     NW_ACCESS = 1 << 0,
     /* Front end/lexical warnings.  */
     NW_LEXICAL = 1 << 1,
     /* Warnings about null pointers.  */
     NW_NONNULL = 1 << 2,
     /* Warnings about uninitialized reads.  */
     NW_UNINIT = 1 << 3,
     /* Warnings about arithmetic overflow.  */
     NW_VFLOW = 1 << 4,
     /* Warnings about dangling pointers.  */
     NW_DANGLING = 1 << 5,
     /* All other unclassified warnings.  */
     NW_OTHER = 1 << 6,
     /* All groups of warnings.  */
     NW_ALL = (NW_ACCESS | NW_LEXICAL | NW_NONNULL
	       | NW_UNINIT | NW_VFLOW | NW_DANGLING | NW_OTHER)
   };

  nowarn_spec_t (): m_bits () { }

  nowarn_spec_t (opt_code);

  /* Return the raw bitset.  */
  operator unsigned() const
  {
    return m_bits;
  }

  /* Return true if the bitset is clear.  */
  bool operator!() const
  {
    return !m_bits;
  }

  /* Return the inverse of the bitset.  */
  nowarn_spec_t operator~() const
  {
    nowarn_spec_t res (*this);
    res.m_bits &= ~NW_ALL;
    return res;
  }

  /* Set *THIS to the bitwise OR of *THIS and RHS.  */
  nowarn_spec_t& operator|= (const nowarn_spec_t &rhs)
  {
    m_bits |= rhs.m_bits;
    return *this;
  }

  /* Set *THIS to the bitwise AND of *THIS and RHS.  */
  nowarn_spec_t& operator&= (const nowarn_spec_t &rhs)
  {
    m_bits &= rhs.m_bits;
    return *this;
  }

  /* Set *THIS to the bitwise exclusive OR of *THIS and RHS.  */
  nowarn_spec_t& operator^= (const nowarn_spec_t &rhs)
  {
    m_bits ^= rhs.m_bits;
    return *this;
  }

private:
  /* Bitset of warning groups.  */
  unsigned m_bits;
};

/* Return the bitwise OR of LHS and RHS.  */

inline nowarn_spec_t
operator| (const nowarn_spec_t &lhs, const nowarn_spec_t &rhs)
{
  return nowarn_spec_t (lhs) |= rhs;
}

/* Return the bitwise AND of LHS and RHS.  */

inline nowarn_spec_t
operator& (const nowarn_spec_t &lhs, const nowarn_spec_t &rhs)
{
  return nowarn_spec_t (lhs) &= rhs;
}

/* Return true if LHS is equal RHS.  */

inline bool
operator== (const nowarn_spec_t &lhs, const nowarn_spec_t &rhs)
{
  return static_cast<unsigned>(lhs) == static_cast<unsigned>(rhs);
}

/* Return true if LHS is not equal RHS.  */

inline bool
operator!= (const nowarn_spec_t &lhs, const nowarn_spec_t &rhs)
{
  return !(lhs == rhs);
}

typedef hash_map<location_hash, nowarn_spec_t> nowarn_map_t;

/* A mapping from a 'location_t' to the warning spec set for it.  */
extern GTY(()) nowarn_map_t *nowarn_map;

#endif // DIAGNOSTIC_SPEC_H_INCLUDED
