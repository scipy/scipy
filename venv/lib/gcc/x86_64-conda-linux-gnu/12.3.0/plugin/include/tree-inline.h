/* Tree inlining hooks and declarations.
   Copyright (C) 2001-2022 Free Software Foundation, Inc.
   Contributed by Alexandre Oliva  <aoliva@redhat.com>

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

#ifndef GCC_TREE_INLINE_H
#define GCC_TREE_INLINE_H


struct cgraph_edge;

/* Indicate the desired behavior wrt call graph edges.  We can either
   duplicate the edge (inlining, cloning), move the edge (versioning,
   parallelization), or move the edges of the clones (saving).  */

enum copy_body_cge_which
{
  CB_CGE_DUPLICATE,
  CB_CGE_MOVE,
  CB_CGE_MOVE_CLONES
};

typedef int_hash <unsigned short, 0> dependence_hash;

/* Data required for function body duplication.  */

struct copy_body_data
{
  /* FUNCTION_DECL for function being inlined, or in general the
     source function providing the original trees.  */
  tree src_fn;

  /* FUNCTION_DECL for function being inlined into, or in general
     the destination function receiving the new trees.  */
  tree dst_fn;

  /* Callgraph node of the source function.  */
  struct cgraph_node *src_node;

  /* Callgraph node of the destination function.  */
  struct cgraph_node *dst_node;

  /* struct function for function being inlined.  Usually this is the same
     as DECL_STRUCT_FUNCTION (src_fn), but can be different if saved_cfg
     and saved_eh are in use.  */
  struct function *src_cfun;

  /* The VAR_DECL for the return value.  */
  tree retvar;

  /* The map from local declarations in the inlined function to
     equivalents in the function into which it is being inlined.  */
  hash_map<tree, tree> *decl_map;

  /* Create a new decl to replace DECL in the destination function.  */
  tree (*copy_decl) (tree, struct copy_body_data *);

  /* Current BLOCK.  */
  tree block;

  /* GIMPLE_CALL if va arg parameter packs should be expanded or NULL
     is not.  */
  gcall *call_stmt;

  /* > 0 if we are remapping a type currently.  */
  int remapping_type_depth;

  /* Exception landing pad the inlined call lies in.  */
  int eh_lp_nr;

  /* Maps region and landing pad structures from the function being copied
     to duplicates created within the function we inline into.  */
  hash_map<void *, void *> *eh_map;

  /* We use the same mechanism do all sorts of different things.  Rather
     than enumerating the different cases, we categorize the behavior
     in the various situations.  */

  /* What to do with call graph edges.  */
  enum copy_body_cge_which transform_call_graph_edges;

  /* True if a new CFG should be created.  False for inlining, true for
     everything else.  */
  bool transform_new_cfg;

  /* True if RETURN_EXPRs should be transformed to just the contained
     MODIFY_EXPR.  The branch semantics of the return will be handled
     by manipulating the CFG rather than a statement.  */
  bool transform_return_to_modify;

  /* True if the parameters of the source function are transformed.
     Only true for inlining.  */
  bool transform_parameter;

  /* True if this statement will need to be regimplified.  */
  bool regimplify;

  /* True if trees may not be unshared.  */
  bool do_not_unshare;

  /* True if trees should not be folded during the copying.  */
  bool do_not_fold;

  /* True if new declarations may not be created during type remapping.  */
  bool prevent_decl_creation_for_types;

  /* True if the location information will need to be reset.  */
  bool reset_location;

  /* Replace error_mark_node as upper bound of array types with
     an uninitialized VAR_DECL temporary.  */
  bool adjust_array_error_bounds;

  /* Usually copy_decl callback always creates new decls, in that case
     we want to remap all variably_modified_type_p types.  If this flag
     is set, remap_type will do further checks to see if remap_decl
     of any decls mentioned in the type will remap to anything but itself
     and only in that case will actually remap the type.  */
  bool dont_remap_vla_if_no_change;

  /* Statements that might be possibly folded.  */
  hash_set<gimple *> *statements_to_fold;

  /* Entry basic block to currently copied body.  */
  basic_block entry_bb;

  /* For partial function versioning, bitmap of bbs to be copied,
     otherwise NULL.  */
  bitmap blocks_to_copy;

  /* Debug statements that need processing.  */
  vec<gdebug *> debug_stmts;

  /* A map from local declarations in the inlined function to
     equivalents in the function into which it is being inlined,
     where the originals have been mapped to a value rather than
     to a variable.  */
  hash_map<tree, tree> *debug_map;

  /* A map from the inlined functions dependence info cliques to
     equivalents in the function into which it is being inlined.  */
  hash_map<dependence_hash, unsigned short> *dependence_map;

  /* A list of addressable local variables remapped into the caller
     when inlining a call within an OpenMP SIMD-on-SIMT loop.  */
  vec<tree> *dst_simt_vars;

  /* Basic block to which clobbers for local variables from the inline
     function that need to live in memory should be added.  */
  basic_block eh_landing_pad_dest;

  /* If clobbers for local variables from the inline function
     that need to live in memory should be added to EH landing pads
     outside of the inlined function, this should be the number
     of basic blocks in the caller before inlining.  Zero otherwise.  */
  int add_clobbers_to_eh_landing_pads;

  /* Class managing changes to function parameters and return value planned
     during IPA stage.  */
  class ipa_param_body_adjustments *param_body_adjs;

  /* Hash set of SSA names that have been killed during call graph edge
   redirection and should not be introduced into debug statements or NULL if no
   SSA_NAME was deleted during redirections happened.  */
  hash_set <tree> *killed_new_ssa_names;
};

/* Weights of constructions for estimate_num_insns.  */

struct eni_weights
{
  /* Cost per call.  */
  unsigned call_cost;

  /* Cost per indirect call.  */
  unsigned indirect_call_cost;

  /* Cost per call to a target specific builtin */
  unsigned target_builtin_call_cost;

  /* Cost of "expensive" div and mod operations.  */
  unsigned div_mod_cost;

  /* Cost for omp construct.  */
  unsigned omp_cost;

  /* Cost for tm transaction.  */
  unsigned tm_cost;

  /* Cost of return.  */
  unsigned return_cost;

  /* True when time of statement should be estimated.  Thus, the
     cost of a switch statement is logarithmic rather than linear in number
     of cases.  */
  bool time_based;
};

/* Weights that estimate_num_insns uses for heuristics in inlining.  */

extern eni_weights eni_inlining_weights;

/* Weights that estimate_num_insns uses to estimate the size of the
   produced code.  */

extern eni_weights eni_size_weights;

/* Weights that estimate_num_insns uses to estimate the time necessary
   to execute the produced code.  */

extern eni_weights eni_time_weights;

/* Function prototypes.  */
void init_inline_once (void);
extern tree copy_tree_body_r (tree *, int *, void *);
extern void insert_decl_map (copy_body_data *, tree, tree);
unsigned int optimize_inline_calls (tree);
tree maybe_inline_call_in_expr (tree);
bool tree_inlinable_function_p (tree);
tree copy_tree_r (tree *, int *, void *);
tree copy_decl_no_change (tree decl, copy_body_data *id);
int estimate_move_cost (tree type, bool);
int estimate_num_insns (gimple *, eni_weights *);
int estimate_num_insns_fn (tree, eni_weights *);
int estimate_num_insns_seq (gimple_seq, eni_weights *);
bool tree_versionable_function_p (tree);
extern tree remap_decl (tree decl, copy_body_data *id);
extern tree remap_type (tree type, copy_body_data *id);
extern gimple_seq copy_gimple_seq_and_replace_locals (gimple_seq seq);
extern bool debug_find_tree (tree, tree);
extern tree copy_fn (tree, tree&, tree&);
extern const char *copy_forbidden (struct function *fun);
extern tree copy_decl_for_dup_finish (copy_body_data *id, tree decl, tree copy);
extern tree copy_decl_to_var (tree, copy_body_data *);
extern tree force_value_to_type (tree type, tree value);

/* This is in tree-inline.cc since the routine uses
   data structures from the inliner.  */
extern tree build_duplicate_type (tree);

#endif /* GCC_TREE_INLINE_H */
