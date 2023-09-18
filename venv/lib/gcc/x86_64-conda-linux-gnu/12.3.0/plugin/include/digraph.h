/* Template classes for directed graphs.
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

#ifndef GCC_DIGRAPH_H
#define GCC_DIGRAPH_H

#include "diagnostic.h"
#include "tree-diagnostic.h" /* for default_tree_printer.  */
#include "graphviz.h"

/* Templates for a family of classes: digraph, node, edge, and cluster.
   This assumes a traits type with the following typedefs:
   node_t: the node class
   edge_t: the edge class
   dump_args_t: additional args for dot-dumps
   cluster_t: the cluster class (for use when generating .dot files).

   Using a template allows for typesafe nodes and edges: a node's
   predecessor and successor edges can be of a node-specific edge
   subclass, without needing casting.  */

/* Abstract base class for a node in a directed graph.  */

template <typename GraphTraits>
class dnode
{
 public:
  typedef typename GraphTraits::edge_t edge_t;
  typedef typename GraphTraits::dump_args_t dump_args_t;

  virtual ~dnode () {}
  virtual void dump_dot (graphviz_out *gv, const dump_args_t &args) const = 0;

  auto_vec<edge_t *> m_preds;
  auto_vec<edge_t *> m_succs;
};

/* Abstract base class for an edge in a directed graph.  */

template <typename GraphTraits>
class dedge
{
 public:
  typedef typename GraphTraits::node_t node_t;
  typedef typename GraphTraits::dump_args_t dump_args_t;

  dedge (node_t *src, node_t *dest)
  : m_src (src), m_dest (dest) {}

  virtual ~dedge () {}

  virtual void dump_dot (graphviz_out *gv, const dump_args_t &args) const = 0;

  node_t *const m_src;
  node_t *const m_dest;
};

/* Abstract base class for a directed graph.
   This class maintains the vectors of nodes and edges,
   and owns the nodes and edges.  */

template <typename GraphTraits>
class digraph
{
 public:
  typedef typename GraphTraits::node_t node_t;
  typedef typename GraphTraits::edge_t edge_t;
  typedef typename GraphTraits::dump_args_t dump_args_t;
  typedef typename GraphTraits::cluster_t cluster_t;

  digraph () {}
  virtual ~digraph () {}

  void dump_dot_to_pp (pretty_printer *pp,
		       cluster_t *root_cluster,
		       const dump_args_t &args) const;
  void dump_dot_to_file (FILE *fp,
			 cluster_t *root_cluster,
			 const dump_args_t &args) const;
  void dump_dot (const char *path,
		 cluster_t *root_cluster,
		 const dump_args_t &args) const;

  void add_node (node_t *node);
  void add_edge (edge_t *edge);

  auto_delete_vec<node_t> m_nodes;
  auto_delete_vec<edge_t> m_edges;
};

/* Abstract base class for splitting dnodes into hierarchical clusters
   in the generated .dot file.

   See "Subgraphs and Clusters" within
     https://www.graphviz.org/doc/info/lang.html
   and e.g.
     https://graphviz.gitlab.io/_pages/Gallery/directed/cluster.html

   If a root_cluster is passed to dump_dot*, then all nodes will be
   added to it at the start of dumping, via calls to add_node.

   The root cluster can organize the nodes into a hierarchy of
   child clusters.

   After all nodes are added to the root cluster, dump_dot will then
   be called on it (and not on the nodes themselves).  */

template <typename GraphTraits>
class cluster
{
 public:
  typedef typename GraphTraits::node_t node_t;
  typedef typename GraphTraits::dump_args_t dump_args_t;

  virtual ~cluster () {}

  virtual void add_node (node_t *node) = 0;

  /* Recursively dump the cluster, all nodes, and child clusters.  */
  virtual void dump_dot (graphviz_out *gv, const dump_args_t &) const = 0;
};

/* Write .dot information for this graph to PP, passing ARGS to the nodes
   and edges.
   If ROOT_CLUSTER is non-NULL, use it to organize the nodes into clusters.  */

template <typename GraphTraits>
inline void
digraph<GraphTraits>::dump_dot_to_pp (pretty_printer *pp,
				      cluster_t *root_cluster,
				      const dump_args_t &args) const
{
  graphviz_out gv (pp);

  pp_string (pp, "digraph \"");
  pp_string (pp, "base");
  pp_string (pp, "\" {\n");

  gv.indent ();

  pp_string (pp, "overlap=false;\n");
  pp_string (pp, "compound=true;\n");

  /* If using clustering, emit all nodes via clusters.  */
  if (root_cluster)
    {
      int i;
      node_t *n;
      FOR_EACH_VEC_ELT (m_nodes, i, n)
	root_cluster->add_node (n);
      root_cluster->dump_dot (&gv, args);
    }
  else
    {
      /* Otherwise, display all nodes at top level.  */
      int i;
      node_t *n;
      FOR_EACH_VEC_ELT (m_nodes, i, n)
	n->dump_dot (&gv, args);
    }

  /* Edges.  */
  int i;
  edge_t *e;
  FOR_EACH_VEC_ELT (m_edges, i, e)
    e->dump_dot (&gv, args);

  /* Terminate "digraph" */
  gv.outdent ();
  pp_string (pp, "}");
  pp_newline (pp);
}

/* Write .dot information for this graph to FP, passing ARGS to the nodes
   and edges.
   If ROOT_CLUSTER is non-NULL, use it to organize the nodes into clusters.  */

template <typename GraphTraits>
inline void
digraph<GraphTraits>::dump_dot_to_file (FILE *fp,
					cluster_t *root_cluster,
					const dump_args_t &args) const
{
  pretty_printer pp;
  // TODO:
  pp_format_decoder (&pp) = default_tree_printer;
  pp.buffer->stream = fp;
  dump_dot_to_pp (&pp, root_cluster, args);
  pp_flush (&pp);
}

/* Write .dot information for this graph to a file at PATH, passing ARGS
   to the nodes and edges.
   If ROOT_CLUSTER is non-NULL, use it to organize the nodes into clusters.  */

template <typename GraphTraits>
inline void
digraph<GraphTraits>::dump_dot (const char *path,
				cluster_t *root_cluster,
				const dump_args_t &args) const
{
  FILE *fp = fopen (path, "w");
  dump_dot_to_file (fp, root_cluster, args);
  fclose (fp);
}

/* Add NODE to this DIGRAPH, taking ownership.  */

template <typename GraphTraits>
inline void
digraph<GraphTraits>::add_node (node_t *node)
{
  m_nodes.safe_push (node);
}

/* Add EDGE to this digraph, and to the preds/succs of its endpoints.
   Take ownership of EDGE.  */

template <typename GraphTraits>
inline void
digraph<GraphTraits>::add_edge (edge_t *edge)
{
  m_edges.safe_push (edge);
  edge->m_dest->m_preds.safe_push (edge);
  edge->m_src->m_succs.safe_push (edge);

}

#endif /* GCC_DIGRAPH_H */
