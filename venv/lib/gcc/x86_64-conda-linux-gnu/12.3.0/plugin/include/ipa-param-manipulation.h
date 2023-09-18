/* Manipulation of formal and actual parameters of functions and function
   calls.
   Copyright (C) 2017-2022 Free Software Foundation, Inc.

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
<http://www.gnu.org/licenses/>.



This file defines classes and other data structures that are used to manipulate
the prototype of a function, especially to create, remove or split its formal
parameters, but also to remove its return value, and also its call statements
correspondingly.

The most basic one is a vector of structures ipa_adjusted_param.  It is simply
a description how the new parameters should look like after the transformation
in what way they relate to the previous ones (if in any).  Such relation to an
old parameter can be an outright copy or an IPA-SRA replacement. If an old
parameter is not listed or otherwise mentioned, it is removed as unused or at
least unnecessary.  Note that this most basic structure does not work for
modifying calls of functions with variable number of arguments.

Class ipa_param_adjustments is only a little more than a thin encapsulation of
a vector of ipa_param_adjustments.  Along with this vector it contains an index
of the first potential vararg argument and a boolean flag whether the return
value should be removed or not.  Moreover, the class contains method
modify_call which can transform a call statement so that it correctly calls a
modified function.  These two data structures were designed to have a small
memory footprint because they are allocated for each clone of a call graph node
that has its prototype changed and live until the end of IPA clone
materialization and call redirection phase.

On the other hand, class ipa_param_body_adjustments can afford to allocate more
data because its life span is much smaller, it is allocated and destroyed in
the course of materialization of each single clone that needs it or only when a
particular pass needs to change a function it is operating on.  This class has
various methods required to change function declaration and the body of the
function according to instructions given either by class ipa_param_adjustments
or only a vector of ipa_adjusted_params.

When these classes are used in the context of call graph clone materialization
and subsequent call statement redirection - which is the point at which we
modify arguments in call statements - they need to cooperate with each other in
order to handle what we refer to as pass-through (IPA-SRA) splits.  These are
situations when a formal parameter of one function is split into several
smaller ones and some of them are then passed on in a call to another function
because the formal parameter of this callee has also been split.

Consider a simple example:

struct S {int a, b, c;};
struct Z {int x; S s;};

foo (S s)
{
  use (s.b);
}

bar (Z z)
{
  use (z.s.a);
  foo (z.s);
}

baz ()
{
  bar (*global);
}

Both bar and foo would have their parameter split.  Foo would receive one
replacement representing s.b.  Function bar would see its parameter split into
one replacement representing z.s.a and another representing z.s.b which would
be passed on to foo.  It would be a so called pass-through split IPA-SRA
replacement, one which is passed in a call as an actual argument to another
IPA-SRA replacement in another function.

Note that the call chain the example can be arbitrarily long and recursive and
that any function in it can be cloned by another IPA pass and any number of
adjacent functions in the call chain can be inlined into each other.  Call
redirection takes place only after bodies of the function have been modified by
all of the above.

Call redirection has to be able to find the right decl or SSA_NAME that
corresponds to the transitive split in the caller.  The SSA names are assigned
right after clone materialization/ modification and cannot be "added" to call
arguments at any later point.  Moreover, if the caller has been inlined the
SSA_NAMEs in question no longer belong to PARM_DECLs but to VAR_DECLs,
indistinguishable from any others.

Therefore, when clone materialization finds a call statement which it knows is
a part of a transitive split, it will simply add as arguments all new "split"
replacements (that have grater or equal offset than the original call
argument):

  foo (repl_for_a, repl_for_b, <rest of original arguments>);

It will also store into ipa_edge_modification_info (which is internal to
ipa-param-modification.c) information about which replacement is which and
where original arguments are.  Call redirection will then invoke
ipa_param_adjustments::modify_call which will access this information and
eliminate all replacements which the callee does not expect (repl_for_a in our
example above).  In between these two steps, however, a call statement might
have extraneous arguments.  */

#ifndef IPA_PARAM_MANIPULATION_H
#define IPA_PARAM_MANIPULATION_H

/* Indices into ipa_param_prefixes to identify a human-readable prefix for newly
   synthesized parameters.  Keep in sync with the array.  */
enum ipa_param_name_prefix_indices
  {
   IPA_PARAM_PREFIX_SYNTH,
   IPA_PARAM_PREFIX_ISRA,
   IPA_PARAM_PREFIX_SIMD,
   IPA_PARAM_PREFIX_MASK,
   IPA_PARAM_PREFIX_COUNT
};

/* We do not support manipulating functions with more than
   1<<IPA_PARAM_MAX_INDEX_BITS parameters.  */
#define IPA_PARAM_MAX_INDEX_BITS 16

/* Operation to be performed for the parameter in ipa_parm_adjustment
   below.  */

enum ipa_parm_op
{
  /* Do not use or you will trigger an assert.  */
  IPA_PARAM_OP_UNDEFINED,

  /* This new parameter is an unmodified parameter at index base_index. */
  IPA_PARAM_OP_COPY,

  /* This describes a brand new parameter.  If it somehow relates to any
     original parameters, the user needs to manage the transition itself.  */
  IPA_PARAM_OP_NEW,

    /* Split parameter as indicated by fields base_index, offset and type.  */
  IPA_PARAM_OP_SPLIT
};

/* Structure that describes one parameter of a function after transformation.
   Omitted parameters will be removed.  */

struct GTY(()) ipa_adjusted_param
{
  /* Type of the new parameter.  Required for all operations except
     IPA_PARM_OP_COPY when the original type will be preserved.  */
  tree type;

  /* Alias reference type to be used in MEM_REFs when adjusting caller
     arguments.  Required for IPA_PARM_OP_SPLIT operation.  */
  tree alias_ptr_type;

  /* Offset into the original parameter (for the cases when the new parameter
     is a component of an original one).  Required for IPA_PARM_OP_SPLIT
     operation.  */
  unsigned unit_offset;

  /* Zero based index of the original parameter this one is based on.  Required
     for IPA_PARAM_OP_COPY and IPA_PARAM_OP_SPLIT, users of IPA_PARAM_OP_NEW
     only need to specify it if they use replacement lookup provided by
     ipa_param_body_adjustments.  */
  unsigned base_index : IPA_PARAM_MAX_INDEX_BITS;

  /* Zero based index of the parameter this one is based on in the previous
     clone.  If there is no previous clone, it must be equal to base_index.  */
  unsigned prev_clone_index : IPA_PARAM_MAX_INDEX_BITS;

  /* Specify the operation, if any, to be performed on the parameter.  */
  enum ipa_parm_op op : 2;

  /* If set, this structure describes a parameter copied over from a previous
     IPA clone, any transformations are thus not to be re-done.  */
  unsigned prev_clone_adjustment : 1;

  /* Index into ipa_param_prefixes specifying a prefix to be used with
     DECL_NAMEs of newly synthesized parameters.  */
  unsigned param_prefix_index : 2;

  /* Storage order of the original parameter (for the cases when the new
     parameter is a component of an original one).  */
  unsigned reverse : 1;

  /* A bit free for the user.  */
  unsigned user_flag : 1;
};

void ipa_dump_adjusted_parameters (FILE *f,
				   vec<ipa_adjusted_param, va_gc> *adj_params);

/* Class used to record planned modifications to parameters of a function and
   also to perform necessary modifications at the caller side at the gimple
   level.  Used to describe all cgraph node clones that have their parameters
   changed, therefore the class should only have a small memory footprint.  */

class GTY(()) ipa_param_adjustments
{
public:
  /* Constructor from NEW_PARAMS showing how new parameters should look like
      plus copying any pre-existing actual arguments starting from argument
      with index ALWAYS_COPY_START (if non-negative, negative means do not copy
      anything beyond what is described in NEW_PARAMS), and SKIP_RETURN, which
      indicates that the function should return void after transformation.  */

  ipa_param_adjustments (vec<ipa_adjusted_param, va_gc> *new_params,
			 int always_copy_start, bool skip_return)
    : m_adj_params (new_params), m_always_copy_start (always_copy_start),
    m_skip_return (skip_return)
    {}

  /* Modify a call statement arguments (and possibly remove the return value)
     as described in the data fields of this class.  */
  gcall *modify_call (cgraph_edge *cs, bool update_references);
  /* Return if the first parameter is left intact.  */
  bool first_param_intact_p ();
  /* Build a function type corresponding to the modified call.  */
  tree build_new_function_type (tree old_type, bool type_is_original_p);
  /* Build a declaration corresponding to the target of the modified call.  */
  tree adjust_decl (tree orig_decl);
  /* Fill a vector marking which parameters are intact by the described
     modifications. */
  void get_surviving_params (vec<bool> *surviving_params);
  /* Fill a vector with new indices of surviving original parameters.  */
  void get_updated_indices (vec<int> *new_indices);
  /* If a parameter with original INDEX has survived intact, return its new
     index.  Otherwise return -1.  In that case, if it has been split and there
     is a new parameter representing a portion at UNIT_OFFSET for which a value
     of a TYPE can be substituted, store its new index into SPLIT_INDEX,
     otherwise store -1 there.  */
  int get_updated_index_or_split (int index, unsigned unit_offset, tree type,
				  int *split_index);
  /* Return the original index for the given new parameter index.  Return a
     negative number if not available.  */
  int get_original_index (int newidx);

  void dump (FILE *f);
  void debug ();

  /* How the known part of arguments should look like.  */
  vec<ipa_adjusted_param, va_gc> *m_adj_params;

  /* If non-negative, copy any arguments starting at this offset without any
     modifications so that functions with variable number of arguments can be
     modified. This number should be equal to the number of original forma
     parameters.  */
  int m_always_copy_start;
  /* If true, make the function not return any value.  */
  bool m_skip_return;

  static bool type_attribute_allowed_p (tree);
private:
  ipa_param_adjustments () {}

  void init (vec<tree> *cur_params);
  int get_max_base_index ();
  bool method2func_p (tree orig_type);
};

/* Structure used to map expressions accessing split or replaced parameters to
   new PARM_DECLs.  */

struct ipa_param_body_replacement
{
  /* The old decl of the original parameter.   */
  tree base;
  /* The new decl it should be replaced with.  */
  tree repl;
  /* Users of ipa_param_body_adjustments that modify standalone functions
     outside of IPA clone materialization can use the following field for their
     internal purposes.  */
  tree dummy;
  /* The offset within BASE that REPL represents.  */
  unsigned unit_offset;
};

struct ipa_replace_map;

/* Class used when actually performing adjustments to formal parameters of a
   function to map accesses that need to be replaced to replacements.  The
   class attempts to work in two very different sets of circumstances: as a
   part of tree-inine.c's tree_function_versioning machinery to clone functions
   (when M_ID is not NULL) and in s standalone fashion, modifying an existing
   function in place (when M_ID is NULL).  While a lot of stuff handled in a
   unified way in both modes, there are many aspects of the processs that
   requires distinct paths.  */

class ipa_param_body_adjustments
{
public:
  /* Constructor to use from within tree-inline.  */
  ipa_param_body_adjustments (ipa_param_adjustments *adjustments,
			      tree fndecl, tree old_fndecl,
			      struct copy_body_data *id, tree *vars,
			      vec<ipa_replace_map *, va_gc> *tree_map);
  /* Constructor to use for modifying a function outside of tree-inline from an
     instance of ipa_param_adjustments.  */
  ipa_param_body_adjustments (ipa_param_adjustments *adjustments,
			      tree fndecl);
  /* Constructor to use for modifying a function outside of tree-inline from a
     simple vector of desired parameter modification.  */
  ipa_param_body_adjustments (vec<ipa_adjusted_param, va_gc> *adj_params,
			      tree fndecl);

  /* The do-it-all function for modifying a function outside of
     tree-inline.  */
  bool perform_cfun_body_modifications ();

  /* Change the PARM_DECLs.  */
  void modify_formal_parameters ();
  /* Register a replacement decl for the transformation done in APM.  */
  void register_replacement (ipa_adjusted_param *apm, tree replacement);
  /* Lookup a replacement for a given offset within a given parameter.  */
  tree lookup_replacement (tree base, unsigned unit_offset);
  /* Lookup a replacement for an expression, if there is one.  */
  ipa_param_body_replacement *get_expr_replacement (tree expr,
						    bool ignore_default_def);
  /* Lookup the new base for surviving names previously belonging to a
     parameter. */
  tree get_replacement_ssa_base (tree old_decl);
  /* Modify a statement.  */
  bool modify_gimple_stmt (gimple **stmt, gimple_seq *extra_stmts,
			   gimple *orig_stmt);
  /* Return the new chain of parameters.  */
  tree get_new_param_chain ();
  /* Replace all occurances of SSAs in m_dead_ssa_debug_equiv in t with what
     they are mapped to.  */
  void remap_with_debug_expressions (tree *t);

  /* Pointers to data structures defining how the function should be
     modified.  */
  vec<ipa_adjusted_param, va_gc> *m_adj_params;
  ipa_param_adjustments *m_adjustments;

  /* Vector of old parameter declarations that must have their debug bind
     statements re-mapped and debug decls created.  */

  auto_vec<tree, 16> m_reset_debug_decls;

  /* Set to true if there are any IPA_PARAM_OP_SPLIT adjustments among stored
     adjustments.  */
  bool m_split_modifications_p;

  /* Sets of statements and SSA_NAMEs that only manipulate data from parameters
     removed because they are not necessary.  */
  hash_set<gimple *> m_dead_stmts;
  hash_set<tree> m_dead_ssas;

  /* Mapping from DCEd SSAs to what their potential debug_binds should be.  */
  hash_map<tree, tree> m_dead_ssa_debug_equiv;
  /* Mapping from DCEd statements to debug expressions that will be placed on
     the RHS of debug statement that will replace this one.  */
  hash_map<gimple *, tree> m_dead_stmt_debug_equiv;

private:
  void common_initialization (tree old_fndecl, tree *vars,
			      vec<ipa_replace_map *, va_gc> *tree_map);
  tree carry_over_param (tree t);
  unsigned get_base_index (ipa_adjusted_param *apm);
  ipa_param_body_replacement *lookup_replacement_1 (tree base,
						    unsigned unit_offset);
  tree replace_removed_params_ssa_names (tree old_name, gimple *stmt);
  bool modify_expression (tree *expr_p, bool convert);
  bool modify_assignment (gimple *stmt, gimple_seq *extra_stmts);
  bool modify_call_stmt (gcall **stmt_p, gimple *orig_stmt);
  bool modify_cfun_body ();
  void reset_debug_stmts ();
  void mark_dead_statements (tree dead_param, vec<tree> *debugstack);
  bool prepare_debug_expressions (tree dead_ssa);

  /* Declaration of the function that is being transformed.  */

  tree m_fndecl;

  /* If non-NULL, the tree-inline master data structure guiding materialization
     of the current clone.  */
  struct copy_body_data *m_id;

  /* Vector of old parameter declarations (before changing them).  */

  auto_vec<tree, 16> m_oparms;

  /* Vector of parameter declarations the function will have after
     transformation.  */

  auto_vec<tree, 16> m_new_decls;

  /* If the function type has non-NULL TYPE_ARG_TYPES, this is the vector of
     these types after transformation, otherwise an empty one.  */

  auto_vec<tree, 16> m_new_types;

  /* Vector of structures telling how to replace old parameters in the
     function body.  TODO: Even though there usually be only few, but should we
     use a hash?  */

  auto_vec<ipa_param_body_replacement, 16> m_replacements;

  /* Vector for remapping SSA_BASES from old parameter declarations that are
     being removed as a part of the transformation.  Before a new VAR_DECL is
     created, it holds the old PARM_DECL, once the variable is built it is
     stored here.  */

  auto_vec<tree> m_removed_decls;

  /* Hash to quickly lookup the item in m_removed_decls given the old decl.  */

  hash_map<tree, unsigned> m_removed_map;

  /* True iff the transformed function is a class method that is about to loose
     its this pointer and must be converted to a normal function.  */

  bool m_method2func;
};

void push_function_arg_decls (vec<tree> *args, tree fndecl);
void push_function_arg_types (vec<tree> *types, tree fntype);
void ipa_verify_edge_has_no_modifications (cgraph_edge *cs);
void ipa_edge_modifications_finalize ();


#endif	/* IPA_PARAM_MANIPULATION_H */
