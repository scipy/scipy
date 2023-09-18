/* Inlining decision heuristics.
   Copyright (C) 2003-2022 Free Software Foundation, Inc.
   Contributed by Jan Hubicka

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

#ifndef GCC_IPA_INLINE_H
#define GCC_IPA_INLINE_H

/* Data we cache about callgraph edges during inlining to avoid expensive
   re-computations during the greedy algorithm.  */
class edge_growth_cache_entry
{
public:
  sreal time, nonspec_time;
  int size;
  ipa_hints hints;

  edge_growth_cache_entry()
    : size (0), hints (0) {}

  edge_growth_cache_entry(int64_t time, int64_t nonspec_time,
			  int size, ipa_hints hints)
    : time (time), nonspec_time (nonspec_time), size (size),
      hints (hints) {}
};

extern fast_call_summary<edge_growth_cache_entry *, va_heap> *edge_growth_cache;

/* In ipa-inline-analysis.cc  */
int estimate_size_after_inlining (struct cgraph_node *, struct cgraph_edge *);
int estimate_growth (struct cgraph_node *);
bool growth_positive_p (struct cgraph_node *, struct cgraph_edge *, int);
int do_estimate_edge_size (struct cgraph_edge *edge);
sreal do_estimate_edge_time (struct cgraph_edge *edge, sreal *nonspec_time = NULL);
ipa_hints do_estimate_edge_hints (struct cgraph_edge *edge);
void reset_node_cache (struct cgraph_node *node);
void initialize_growth_caches ();
void free_growth_caches (void);

/* In ipa-inline.cc  */
unsigned int early_inliner (function *fun);
bool inline_account_function_p (struct cgraph_node *node);


/* In ipa-inline-transform.cc  */
bool inline_call (struct cgraph_edge *, bool, vec<cgraph_edge *> *, int *, bool,
		  bool *callee_removed = NULL);
unsigned int inline_transform (struct cgraph_node *);
void clone_inlined_nodes (struct cgraph_edge *e, bool, bool, int *);

extern int ncalls_inlined;
extern int nfunctions_inlined;
extern function_summary <tree *> *ipa_saved_clone_sources;

/* Return estimated size of the inline sequence of EDGE.  */

static inline int
estimate_edge_size (struct cgraph_edge *edge)
{
  edge_growth_cache_entry *entry;
  if (edge_growth_cache == NULL
      || (entry = edge_growth_cache->get (edge)) == NULL
      || entry->size == 0)
    return do_estimate_edge_size (edge);
  return entry->size - (entry->size > 0);
}

/* Return lower bound on estimated callee growth after inlining EDGE.  */

static inline int
estimate_min_edge_growth (struct cgraph_edge *edge)
{
  ipa_call_summary *s = ipa_call_summaries->get (edge);
  struct cgraph_node *callee = edge->callee->ultimate_alias_target ();
  return (ipa_fn_summaries->get (callee)->min_size - s->call_stmt_size);
}

/* Return estimated callee growth after inlining EDGE.  */

static inline int
estimate_edge_growth (struct cgraph_edge *edge)
{
  ipa_call_summary *s = ipa_call_summaries->get (edge);
  gcc_checking_assert (s->call_stmt_size || !edge->callee->analyzed);
  return (estimate_edge_size (edge) - s->call_stmt_size);
}

/* Return estimated callee runtime increase after inlining
   EDGE.  */

static inline sreal
estimate_edge_time (struct cgraph_edge *edge, sreal *nonspec_time = NULL)
{
  edge_growth_cache_entry *entry;
  if (edge_growth_cache == NULL
      || (entry = edge_growth_cache->get (edge)) == NULL
      || entry->time == 0)
    return do_estimate_edge_time (edge, nonspec_time);
  if (nonspec_time)
    *nonspec_time = edge_growth_cache->get (edge)->nonspec_time;
  return entry->time;
}


/* Return estimated callee runtime increase after inlining
   EDGE.  */

static inline ipa_hints
estimate_edge_hints (struct cgraph_edge *edge)
{
  edge_growth_cache_entry *entry;
  if (edge_growth_cache == NULL
      || (entry = edge_growth_cache->get (edge)) == NULL
      || entry->hints == 0)
    return do_estimate_edge_hints (edge);
  return entry->hints - 1;
}

#endif /* GCC_IPA_INLINE_H */
