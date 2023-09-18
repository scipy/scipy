/* Template class for Dijkstra's algorithm on directed graphs.
   Copyright (C) 2019-2022 Free Software Foundation, Inc.
   Contributed by David Malcolm <dmalcolm@redhat.com>.

This file is part of GCC.

GCC is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3, or (at your option)
any later version.

GCC is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with GCC; see the file COPYING3.  If not see
<http://www.gnu.org/licenses/>.  */

#ifndef GCC_SHORTEST_PATHS_H
#define GCC_SHORTEST_PATHS_H

#include "timevar.h"

enum shortest_path_sense
{
  /* Find the shortest path from the given origin node to each
     node in the graph.  */
  SPS_FROM_GIVEN_ORIGIN,

  /* Find the shortest path from each node in the graph to the
     given target node.  */
  SPS_TO_GIVEN_TARGET
};

/* A record of the shortest path for each node relative to a special
   "given node", either:
   SPS_FROM_GIVEN_ORIGIN:
     from the given origin node to each node in a graph, or
   SPS_TO_GIVEN_TARGET:
     from each node in a graph to the given target node.

   The constructor runs Dijkstra's algorithm, and the results are
   stored in this class.  */

template <typename GraphTraits, typename Path_t>
class shortest_paths
{
public:
  typedef typename GraphTraits::graph_t graph_t;
  typedef typename GraphTraits::node_t node_t;
  typedef typename GraphTraits::edge_t edge_t;
  typedef Path_t path_t;

  shortest_paths (const graph_t &graph, const node_t *given_node,
		  enum shortest_path_sense sense);

  path_t get_shortest_path (const node_t *other_node) const;
  int get_shortest_distance (const node_t *other_node) const;

private:
  const graph_t &m_graph;

  enum shortest_path_sense m_sense;

  /* For each node (by index), the minimal distance between that node
     and the given node (with direction depending on m_sense).  */
  auto_vec<int> m_dist;

  /* For each node (by index):
     SPS_FROM_GIVEN_ORIGIN:
       the previous edge in the shortest path from the origin,
     SPS_TO_GIVEN_TARGET:
       the next edge in the shortest path to the target.  */
  auto_vec<const edge_t *> m_best_edge;
};

/* shortest_paths's constructor.

   Use Dijkstra's algorithm relative to GIVEN_NODE to populate m_dist and
   m_best_edge with enough information to be able to generate Path_t instances
   to give the shortest path...
   SPS_FROM_GIVEN_ORIGIN: to each node in a graph from the origin node, or
   SPS_TO_GIVEN_TARGET: from each node in a graph to the target node.  */

template <typename GraphTraits, typename Path_t>
inline
shortest_paths<GraphTraits, Path_t>::
shortest_paths (const graph_t &graph,
		const node_t *given_node,
		enum shortest_path_sense sense)
: m_graph (graph),
  m_sense (sense),
  m_dist (graph.m_nodes.length ()),
  m_best_edge (graph.m_nodes.length ())
{
  auto_timevar tv (TV_ANALYZER_SHORTEST_PATHS);

  auto_vec<int> queue (graph.m_nodes.length ());

  for (unsigned i = 0; i < graph.m_nodes.length (); i++)
    {
      m_dist.quick_push (INT_MAX);
      m_best_edge.quick_push (NULL);
      queue.quick_push (i);
    }
  m_dist[given_node->m_index] = 0;

  while (queue.length () > 0)
    {
      /* Get minimal distance in queue.
	 FIXME: this is O(N^2); replace with a priority queue.  */
      int idx_with_min_dist = -1;
      int idx_in_queue_with_min_dist = -1;
      int min_dist = INT_MAX;
      for (unsigned i = 0; i < queue.length (); i++)
	{
	  int idx = queue[i];
	  if (m_dist[queue[i]] < min_dist)
	    {
	      min_dist = m_dist[idx];
	      idx_with_min_dist = idx;
	      idx_in_queue_with_min_dist = i;
	    }
	}
      if (idx_with_min_dist == -1)
	break;
      gcc_assert (idx_in_queue_with_min_dist != -1);

      // FIXME: this is confusing: there are two indices here

      queue.unordered_remove (idx_in_queue_with_min_dist);

      node_t *n
	= static_cast <node_t *> (m_graph.m_nodes[idx_with_min_dist]);

      if (m_sense == SPS_FROM_GIVEN_ORIGIN)
	{
	  int i;
	  edge_t *succ;
	  FOR_EACH_VEC_ELT (n->m_succs, i, succ)
	    {
	      // TODO: only for dest still in queue
	      node_t *dest = succ->m_dest;
	      int alt = m_dist[n->m_index] + 1;
	      if (alt < m_dist[dest->m_index])
		{
		  m_dist[dest->m_index] = alt;
		  m_best_edge[dest->m_index] = succ;
		}
	    }
	}
      else
	{
	  int i;
	  edge_t *pred;
	  FOR_EACH_VEC_ELT (n->m_preds, i, pred)
	    {
	      // TODO: only for dest still in queue
	      node_t *src = pred->m_src;
	      int alt = m_dist[n->m_index] + 1;
	      if (alt < m_dist[src->m_index])
		{
		  m_dist[src->m_index] = alt;
		  m_best_edge[src->m_index] = pred;
		}
	    }
	}
   }
}

/* Generate an Path_t instance giving the shortest path between OTHER_NODE
   and the given node.

   SPS_FROM_GIVEN_ORIGIN: shortest path from given origin node to OTHER_NODE
   SPS_TO_GIVEN_TARGET: shortest path from OTHER_NODE to given target node.

   If no such path exists, return an empty path.  */

template <typename GraphTraits, typename Path_t>
inline Path_t
shortest_paths<GraphTraits, Path_t>::
get_shortest_path (const node_t *other_node) const
{
  Path_t result;

  while (m_best_edge[other_node->m_index])
    {
      result.m_edges.safe_push (m_best_edge[other_node->m_index]);
      if (m_sense == SPS_FROM_GIVEN_ORIGIN)
	other_node = m_best_edge[other_node->m_index]->m_src;
      else
	other_node = m_best_edge[other_node->m_index]->m_dest;
    }

  if (m_sense == SPS_FROM_GIVEN_ORIGIN)
    result.m_edges.reverse ();

  return result;
}

/* Get the shortest distance...
   SPS_FROM_GIVEN_ORIGIN: ...from given origin node to OTHER_NODE
   SPS_TO_GIVEN_TARGET: ...from OTHER_NODE to given target node.  */

template <typename GraphTraits, typename Path_t>
inline int
shortest_paths<GraphTraits, Path_t>::
get_shortest_distance (const node_t *other_node) const
{
  return m_dist[other_node->m_index];
}

#endif /* GCC_SHORTEST_PATHS_H */
