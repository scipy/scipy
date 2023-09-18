/* IPA function body analysis.
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

#ifndef GCC_IPA_SUMMARY_H
#define GCC_IPA_SUMMARY_H

#include "sreal.h"
#include "ipa-predicate.h"


/* Hints are reasons why IPA heuristics should prefer specializing given
   function.  They are represented as bitmap of the following values.  */
enum ipa_hints_vals {
  /* When specialization turns indirect call into a direct call,
     it is good idea to do so.  */
  INLINE_HINT_indirect_call = 1,
  /* Inlining may make loop iterations or loop stride known.  It is good idea
     to do so because it enables loop optimizations.  */
  INLINE_HINT_loop_iterations = 2,
  INLINE_HINT_loop_stride = 4,
  /* Inlining within same strongly connected component of callgraph is often
     a loss due to increased stack frame usage and prologue setup costs.  */
  INLINE_HINT_same_scc = 8,
  /* Inlining functions in strongly connected component is not such a great
     win.  */
  INLINE_HINT_in_scc = 16,
  /* If function is declared inline by user, it may be good idea to inline
     it.  Set by simple_edge_hints in ipa-inline-analysis.cc.  */
  INLINE_HINT_declared_inline = 32,
  /* Programs are usually still organized for non-LTO compilation and thus
     if functions are in different modules, inlining may not be so important. 
     Set by simple_edge_hints in ipa-inline-analysis.cc.   */
  INLINE_HINT_cross_module = 64,
  /* We know that the callee is hot by profile.  */
  INLINE_HINT_known_hot = 128,
  /* There is builtin_constant_p dependent on parameter which is usually
     a strong hint to inline.  */
  INLINE_HINT_builtin_constant_p = 256
};

typedef int ipa_hints;

/* Simple description of whether a memory load or a condition refers to a load
   from an aggregate and if so, how and where from in the aggregate.
   Individual fields have the same meaning like fields with the same name in
   struct condition.  */

struct agg_position_info
{
  HOST_WIDE_INT offset;
  bool agg_contents;
  bool by_ref;
};

/* Representation of function body size and time depending on the call
   context.  We keep simple array of record, every containing of predicate
   and time/size to account.  */
class size_time_entry
{
public:
  /* Predicate for code to be executed.  */
  ipa_predicate exec_predicate;
  /* Predicate for value to be constant and optimized out in a specialized copy.
     When deciding on specialization this makes it possible to see how much
     the executed code paths will simplify.  */
  ipa_predicate nonconst_predicate;
  int size;
  sreal time;
};

/* Summary about function and stack frame sizes.  We keep this info 
   for inline clones and also for WPA streaming. For this reason this is not
   part of ipa_fn_summary which exists only for offline functions.  */
class ipa_size_summary
{
public:
  /* Estimated stack frame consumption by the function.  */
  HOST_WIDE_INT estimated_self_stack_size;
  /* Size of the function body.  */
  int self_size;
  /* Estimated size of the function after inlining.  */
  int size;

  ipa_size_summary ()
  : estimated_self_stack_size (0), self_size (0), size (0)
  {
  }
};

/* Structure to capture how frequently some interesting events occur given a
   particular predicate.  The structure is used to estimate how often we
   encounter loops with known iteration count or stride in various
   contexts.  */

struct GTY(()) ipa_freqcounting_predicate
{
  /* The described event happens with this frequency... */
  sreal freq;
  /* ...when this predicate evaluates to false. */
  ipa_predicate * GTY((skip)) predicate;
};

/* Function inlining information.  */
class GTY(()) ipa_fn_summary
{
public:
  /* Keep all field empty so summary dumping works during its computation.
     This is useful for debugging.  */
  ipa_fn_summary ()
    : min_size (0),
      inlinable (false), single_caller (false),
      fp_expressions (false), target_info (0),
      estimated_stack_size (false),
      time (0), conds (NULL),
      size_time_table (), call_size_time_table (vNULL),
      loop_iterations (NULL), loop_strides (NULL),
      builtin_constant_p_parms (vNULL),
      growth (0), scc_no (0)
  {
  }

  /* Copy constructor.  */
  ipa_fn_summary (const ipa_fn_summary &s)
    : min_size (s.min_size),
    inlinable (s.inlinable), single_caller (s.single_caller),
    fp_expressions (s.fp_expressions),
    target_info (s.target_info),
    estimated_stack_size (s.estimated_stack_size),
    time (s.time), conds (s.conds), size_time_table (),
    call_size_time_table (vNULL),
    loop_iterations (s.loop_iterations), loop_strides (s.loop_strides),
    builtin_constant_p_parms (s.builtin_constant_p_parms),
    growth (s.growth), scc_no (s.scc_no)
  {}

  /* Default constructor.  */
  ~ipa_fn_summary ();

  /* Information about the function body itself.  */

  /* Minimal size increase after inlining.  */
  int min_size;

  /* False when there something makes inlining impossible (such as va_arg).  */
  unsigned inlinable : 1;
  /* True wen there is only one caller of the function before small function
     inlining.  */
  unsigned int single_caller : 1;
  /* True if function contains any floating point expressions.  */
  unsigned int fp_expressions : 1;
  /* Like fp_expressions field above, but it's to hold some target specific
     information, such as some target specific isa flags.  Note that for
     offloading target compilers, this field isn't streamed.  */
  unsigned int target_info;

  /* Information about function that will result after applying all the
     inline decisions present in the callgraph.  Generally kept up to
     date only for functions that are not inline clones. */

  /* Estimated stack frame consumption by the function.  */
  HOST_WIDE_INT estimated_stack_size;
  /* Estimated runtime of function after inlining.  */
  sreal GTY((skip)) time;

  /* Conditional size/time information.  The summaries are being
     merged during inlining.  */
  conditions conds;
  /* Normal code is accounted in size_time_table, while calls are
     accounted in call_size_time_table.  This is because calls
     are often adjusted by IPA optimizations and thus this summary
     is generated from call summary information when needed.  */
  auto_vec<size_time_entry> GTY((skip)) size_time_table;
  /* Unlike size_time_table that is initialized for all summaries
     call_size_time_table is allocated only for functions with
     many calls.  Use effecient vl_ptr storage.  */
  vec<size_time_entry, va_heap, vl_ptr> GTY((skip)) call_size_time_table;

  /* Predicates on when some loops in the function can have known bounds.  */
  vec<ipa_freqcounting_predicate, va_gc> *loop_iterations;
  /* Predicates on when some loops in the function can have known strides.  */
  vec<ipa_freqcounting_predicate, va_gc> *loop_strides;
  /* Parameters tested by builtin_constant_p.  */
  vec<int, va_heap, vl_ptr> GTY((skip)) builtin_constant_p_parms;
  /* Estimated growth for inlining all copies of the function before start
     of small functions inlining.
     This value will get out of date as the callers are duplicated, but
     using up-to-date value in the badness metric mean a lot of extra
     expenses.  */
  int growth;
  /* Number of SCC on the beginning of inlining process.  */
  int scc_no;

  /* Record time and size under given predicates.  */
  void account_size_time (int, sreal, const ipa_predicate &,
			  const ipa_predicate &,
		  	  bool call = false);

  /* We keep values scaled up, so fractional sizes can be accounted.  */
  static const int size_scale = 2;
  /* Maximal size of size_time_table before we start to be conservative.  */
  static const int max_size_time_table_size = 256;
};

class GTY((user)) ipa_fn_summary_t:
  public fast_function_summary <ipa_fn_summary *, va_gc>
{
public:
  ipa_fn_summary_t (symbol_table *symtab):
    fast_function_summary <ipa_fn_summary *, va_gc> (symtab) {}

  static ipa_fn_summary_t *create_ggc (symbol_table *symtab)
  {
    class ipa_fn_summary_t *summary
      = new (ggc_alloc_no_dtor<ipa_fn_summary_t> ()) ipa_fn_summary_t (symtab);
    summary->disable_insertion_hook ();
    return summary;
  }

  /* Remove ipa_fn_summary for all callees of NODE.  */
  void remove_callees (cgraph_node *node);

  virtual void insert (cgraph_node *, ipa_fn_summary *);
  virtual void remove (cgraph_node *node, ipa_fn_summary *)
  {
    remove_callees (node);
  }

  virtual void duplicate (cgraph_node *src, cgraph_node *dst,
			  ipa_fn_summary *src_data, ipa_fn_summary *dst_data);
};

extern GTY(()) fast_function_summary <ipa_fn_summary *, va_gc>
  *ipa_fn_summaries;

class ipa_size_summary_t:
  public fast_function_summary <ipa_size_summary *, va_heap>
{
public:
  ipa_size_summary_t (symbol_table *symtab):
    fast_function_summary <ipa_size_summary *, va_heap> (symtab)
  {
    disable_insertion_hook ();
  }

  virtual void duplicate (cgraph_node *, cgraph_node *,
			  ipa_size_summary *src_data,
			  ipa_size_summary *dst_data)
  {
    *dst_data = *src_data;
  }
};
extern fast_function_summary <ipa_size_summary *, va_heap>
  *ipa_size_summaries;

/* Information kept about callgraph edges.  */
class ipa_call_summary
{
public:
  /* Keep all field empty so summary dumping works during its computation.
     This is useful for debugging.  */
  ipa_call_summary ()
    : predicate (NULL), param (vNULL), call_stmt_size (0), call_stmt_time (0),
      loop_depth (0), is_return_callee_uncaptured (false)
    {
    }

  /* Copy constructor.  */
  ipa_call_summary (const ipa_call_summary &s):
    predicate (s.predicate), param (s.param), call_stmt_size (s.call_stmt_size),
    call_stmt_time (s.call_stmt_time), loop_depth (s.loop_depth),
    is_return_callee_uncaptured (s.is_return_callee_uncaptured)
  {
  }

  /* Default destructor.  */
  ~ipa_call_summary ();

  ipa_predicate *predicate;
  /* Vector indexed by parameters.  */
  vec<inline_param_summary> param;
  /* Estimated size and time of the call statement.  */
  int call_stmt_size;
  int call_stmt_time;
  /* Depth of loop nest, 0 means no nesting.  */
  unsigned int loop_depth;
  /* Indicates whether the caller returns the value of it's callee.  */
  bool is_return_callee_uncaptured;
};

class ipa_call_summary_t: public fast_call_summary <ipa_call_summary *, va_heap>
{
public:
  ipa_call_summary_t (symbol_table *symtab):
    fast_call_summary <ipa_call_summary *, va_heap> (symtab) {}

  /* Hook that is called by summary when an edge is duplicated.  */
  virtual void duplicate (cgraph_edge *src, cgraph_edge *dst,
			  ipa_call_summary *src_data,
			  ipa_call_summary *dst_data);
};

/* Estimated execution times, code sizes and other information about the
   code executing a call described by ipa_call_context.  */

struct ipa_call_estimates
{
  /* Estimated size needed to execute call in the given context. */
  int size;

  /* Minimal size needed for the call that is + independent on the call context
     and can be used for fast estimates.  */
  int min_size;

  /* Estimated time needed to execute call in the given context. */
  sreal time;

  /* Estimated time needed to execute the function when not ignoring
     computations known to be constant in this context.  */
  sreal nonspecialized_time;

  /* Further discovered reasons why to inline or specialize the give calls.  */
  ipa_hints hints;

  /* Frequency how often a loop with known number of iterations is encountered.
     Calculated with hints.  */
  sreal loops_with_known_iterations;

  /* Frequency how often a loop with known strides is encountered.  Calculated
     with hints.  */
  sreal loops_with_known_strides;
};

class ipa_cached_call_context;

/* This object describe a context of call.  That is a summary of known
   information about its parameters.  Main purpose of this context is
   to give more realistic estimations of function runtime, size and
   inline hints.  */
class ipa_call_context
{
public:
  ipa_call_context (cgraph_node *node,
      		    clause_t possible_truths,
		    clause_t nonspec_possible_truths,
		    vec<inline_param_summary> inline_param_summary,
		    ipa_auto_call_arg_values *arg_values);
  ipa_call_context ()
  : m_node(NULL)
  {
  }
  void estimate_size_and_time (ipa_call_estimates *estimates,
			       bool est_times = true, bool est_hints = true);
  bool equal_to (const ipa_call_context &);
  bool exists_p ()
  {
    return m_node != NULL;
  }
private:
  /* Called function.  */
  cgraph_node *m_node;
  /* Clause describing what predicate conditionals can be satisfied
     in this context if function is inlined/specialized.  */
  clause_t m_possible_truths;
  /* Clause describing what predicate conditionals can be satisfied
     in this context if function is kept offline.  */
  clause_t m_nonspec_possible_truths;
  /* Inline summary maintains info about change probabilities.  */
  vec<inline_param_summary> m_inline_param_summary;

  /* Even after having calculated clauses, the information about argument
     values is used to resolve indirect calls.  */
  ipa_call_arg_values m_avals;

  friend ipa_cached_call_context;
};

/* Variant of ipa_call_context that is stored in a cache over a longer period
   of time.  */

class ipa_cached_call_context : public ipa_call_context
{
public:
  void duplicate_from (const ipa_call_context &ctx);
  void release ();
};

extern fast_call_summary <ipa_call_summary *, va_heap> *ipa_call_summaries;

/* In ipa-fnsummary.cc  */
void ipa_debug_fn_summary (struct cgraph_node *);
void ipa_dump_fn_summaries (FILE *f);
void ipa_dump_fn_summary (FILE *f, struct cgraph_node *node);
void ipa_dump_hints (FILE *f, ipa_hints);
void ipa_free_fn_summary (void);
void ipa_free_size_summary (void);
void inline_analyze_function (struct cgraph_node *node);
void estimate_ipcp_clone_size_and_time (struct cgraph_node *node,
					ipa_auto_call_arg_values *avals,
					ipa_call_estimates *estimates);
void ipa_merge_fn_summary_after_inlining (struct cgraph_edge *edge);
void ipa_update_overall_fn_summary (struct cgraph_node *node, bool reset = true);
void compute_fn_summary (struct cgraph_node *, bool);
bool refs_local_or_readonly_memory_p (tree);
bool points_to_local_or_readonly_memory_p (tree);


void evaluate_properties_for_edge (struct cgraph_edge *e,
	       		           bool inline_p,
				   clause_t *clause_ptr,
				   clause_t *nonspec_clause_ptr,
				   ipa_auto_call_arg_values *avals,
				   bool compute_contexts);

void ipa_fnsummary_cc_finalize (void);
HOST_WIDE_INT ipa_get_stack_frame_offset (struct cgraph_node *node);
void ipa_remove_from_growth_caches (struct cgraph_edge *edge);

/* Return true if EDGE is a cross module call.  */

static inline bool
cross_module_call_p (struct cgraph_edge *edge)
{
  /* Here we do not want to walk to alias target becuase ICF may create
     cross-unit aliases.  */
  if (edge->caller->unit_id == edge->callee->unit_id)
    return false;
  /* If the call is to a (former) comdat function or s symbol with mutiple
     extern inline definitions then treat is as in-module call.  */
  if (edge->callee->merged_extern_inline || edge->callee->merged_comdat
      || DECL_COMDAT (edge->callee->decl))
    return false;
  return true;
}

#endif /* GCC_IPA_FNSUMMARY_H */
