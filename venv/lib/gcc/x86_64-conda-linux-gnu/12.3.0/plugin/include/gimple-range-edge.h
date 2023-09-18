/* Gimple range edge header file.
   Copyright (C) 2020-2022 Free Software Foundation, Inc.
   Contributed by Andrew MacLeod <amacleod@redhat.com>
   and Aldy Hernandez <aldyh@redhat.com>.

This file is part of GCC.

GCC is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3, or (at your option)
any later version.

GCC is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GCC; see the file COPYING3.  If not see
<http://www.gnu.org/licenses/>.  */

#ifndef GIMPLE_RANGE_EDGE_H
#define GIMPLE_RANGE_EDGE_H

// This class is used to query ranges on constant edges in GIMPLE.
//
// For a COND_EXPR, the TRUE edge will return [1,1] and the false edge a [0,0].
//
// For SWITCH_EXPR, it is awkward to calculate ranges.  When a request
// is made, the entire switch is evalauted and the results cached. 
// Any future requests to that switch will use the cached value, providing
// dramatic decrease in computation time.
//
// The API is simple, just ask for the range on the edge.
// The return value is NULL for no range, or the branch statement which the
// edge gets the range from, along with the range.

class gimple_outgoing_range
{
public:
  gimple_outgoing_range (int max_sw_edges = INT_MAX);
  ~gimple_outgoing_range ();
  gimple *edge_range_p (irange &r, edge e);
private:
  void calc_switch_ranges (gswitch *sw);
  bool get_edge_range (irange &r, gimple *s, edge e);

  int m_max_edges;
  hash_map<edge, irange *> *m_edge_table;
  irange_allocator m_range_allocator;
};

// If there is a range control statement at the end of block BB, return it.
gimple *gimple_outgoing_range_stmt_p (basic_block bb);
// Return the range on edge E if it is from a GCOND.  Either TRUE or FALSE.
void gcond_edge_range (irange &r, edge e);

#endif  // GIMPLE_RANGE_EDGE_H
