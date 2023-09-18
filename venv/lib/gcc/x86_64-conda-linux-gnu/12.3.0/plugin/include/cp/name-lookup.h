/* Declarations for -*- C++ -*- name lookup routines.
   Copyright (C) 2003-2022 Free Software Foundation, Inc.
   Contributed by Gabriel Dos Reis <gdr@integrable-solutions.net>

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

#ifndef GCC_CP_NAME_LOOKUP_H
#define GCC_CP_NAME_LOOKUP_H

#include "c-family/c-common.h"


/* The datatype used to implement C++ scope.  */
struct cp_binding_level;

/* Nonzero if this binding is for a local scope, as opposed to a class
   or namespace scope.  */
#define LOCAL_BINDING_P(NODE) ((NODE)->is_local)

/* True if NODE->value is from a base class of the class which is
   currently being defined.  */
#define INHERITED_VALUE_BINDING_P(NODE) ((NODE)->value_is_inherited)

/* The IMPLICIT_TYPEDEF is hidden from ordinary name lookup (it was
   injected via a local class's friend decl). The typdef may be in the
   VALUE or the TYPE slot.  We do not get the situation where the
   value and type slots are both filled and both hidden.  */
#define HIDDEN_TYPE_BINDING_P(NODE) ((NODE)->type_is_hidden)

/* Datatype that represents binding established by a declaration between
   a name and a C++ entity.  */
struct GTY(()) cxx_binding {
  /* Link to chain together various bindings for this name.  */
  cxx_binding *previous;
  /* The non-type entity this name is bound to.  */
  tree value;
  /* The type entity this name is bound to.  */
  tree type;
  /* The scope at which this binding was made.  */
  cp_binding_level *scope;

  bool value_is_inherited : 1;
  bool is_local : 1;
  bool type_is_hidden : 1;
};

/* Datatype used to temporarily save C++ bindings (for implicit
   instantiations purposes and like).  Implemented in decl.cc.  */
struct GTY(()) cxx_saved_binding {
  /* The name of the current binding.  */
  tree identifier;
  /* The binding we're saving.  */
  cxx_binding *binding;
  tree real_type_value;
};

/* To support lazy module loading, we squirrel away a section number
   (and a couple of flags) in the binding slot of unloaded bindings.
   We rely on pointers being aligned and setting the bottom bit to
   mark a lazy value.  GTY doesn't like an array of union, so we have
   a containing struct.  */

struct GTY(()) binding_slot {
  union GTY((desc ("%1.is_lazy ()"))) binding_slot_lazy {
    tree GTY((tag ("false"))) binding;
  } u;

  operator tree & ()
  {
    gcc_checking_assert (!is_lazy ());
    return u.binding;
  }
  binding_slot &operator= (tree t)
  {
    u.binding = t;
    return *this;
  }
  bool is_lazy () const
  {
    return bool (uintptr_t (u.binding) & 1);
  }
  void set_lazy (unsigned snum)
  {
    gcc_checking_assert (!u.binding);
    u.binding = tree (uintptr_t ((snum << 1) | 1));
  }
  void or_lazy (unsigned snum)
  {
    gcc_checking_assert (is_lazy ());
    u.binding = tree (uintptr_t (u.binding) | (snum << 1));
  }
  unsigned get_lazy () const
  {
    gcc_checking_assert (is_lazy ());
    return unsigned (uintptr_t (u.binding) >> 1);
  }
};

/* Bindings for modules are held in a sparse array.  There is always a
   current TU slot, others are allocated as needed.  By construction
   of the importing mechanism we only ever need to append to the
   array.  Rather than have straight index/slot tuples, we bunch them
   up for greater packing.

   The cluster representation packs well on a 64-bit system.  */

#define BINDING_VECTOR_SLOTS_PER_CLUSTER 2
struct binding_index {
  unsigned short base;
  unsigned short span;
};

struct GTY(()) binding_cluster
{
  binding_index GTY((skip)) indices[BINDING_VECTOR_SLOTS_PER_CLUSTER];
  binding_slot slots[BINDING_VECTOR_SLOTS_PER_CLUSTER];
};

/* These two fields overlay lang flags.  So don't use those.  */
#define BINDING_VECTOR_ALLOC_CLUSTERS(NODE) \
  (BINDING_VECTOR_CHECK (NODE)->base.u.dependence_info.clique)
#define BINDING_VECTOR_NUM_CLUSTERS(NODE) \
  (BINDING_VECTOR_CHECK (NODE)->base.u.dependence_info.base)
#define BINDING_VECTOR_CLUSTER_BASE(NODE) \
  (((tree_binding_vec *)BINDING_VECTOR_CHECK (NODE))->vec)
#define BINDING_VECTOR_CLUSTER_LAST(NODE) \
  (&BINDING_VECTOR_CLUSTER (NODE, BINDING_VECTOR_NUM_CLUSTERS (NODE) - 1))
#define BINDING_VECTOR_CLUSTER(NODE,IX) \
  (((tree_binding_vec *)BINDING_VECTOR_CHECK (NODE))->vec[IX])

struct GTY(()) tree_binding_vec {
  struct tree_base base;
  tree name;
  binding_cluster GTY((length ("%h.base.u.dependence_info.base"))) vec[1];
};

/* The name of a module vector.  */
#define BINDING_VECTOR_NAME(NODE) \
  (((tree_binding_vec *)BINDING_VECTOR_CHECK (NODE))->name)

/* tree_binding_vec does uses  base.u.dependence_info.base field for
   length.  It does not have lang_flag etc available!  */

/* These two flags note if a module-vector contains deduplicated
   bindings (i.e. multiple declarations in different imports).  */
/* This binding contains duplicate references to a global module
   entity.  */
#define BINDING_VECTOR_GLOBAL_DUPS_P(NODE) \
  (BINDING_VECTOR_CHECK (NODE)->base.static_flag)
/* This binding contains duplicate references to a partioned module
   entity.  */
#define BINDING_VECTOR_PARTITION_DUPS_P(NODE) \
  (BINDING_VECTOR_CHECK (NODE)->base.volatile_flag)

/* These two flags indicate the provenence of the bindings on this
   particular vector slot.  We can of course determine this from slot
   number, but that's a relatively expensive lookup.  This avoids
   that when iterating.  */
/* This slot is part of the global module (a header unit).  */
#define MODULE_BINDING_GLOBAL_P(NODE) \
  (OVERLOAD_CHECK (NODE)->base.static_flag)
/* This slot is part of the current module (a partition or primary).  */
#define MODULE_BINDING_PARTITION_P(NODE)		\
  (OVERLOAD_CHECK (NODE)->base.volatile_flag)

extern void set_identifier_type_value (tree, tree);
extern void push_binding (tree, tree, cp_binding_level*);
extern void pop_local_binding (tree, tree);
extern void pop_bindings_and_leave_scope (void);
extern tree constructor_name (tree);
extern bool constructor_name_p (tree, tree);

/* The kinds of scopes we recognize.  */
enum scope_kind {
  sk_block = 0,      /* An ordinary block scope.  This enumerator must
			have the value zero because "cp_binding_level"
			is initialized by using "memset" to set the
			contents to zero, and the default scope kind
			is "sk_block".  */
  sk_cleanup,	     /* A scope for (pseudo-)scope for cleanup.  It is
			pseudo in that it is transparent to name lookup
			activities.  */
  sk_try,	     /* A try-block.  */
  sk_catch,	     /* A catch-block.  */
  sk_for,	     /* The scope of the variable declared in a
			init-statement.  */
  sk_cond,	     /* The scope of the variable declared in the condition
			of an if or switch statement.  */
  sk_function_parms, /* The scope containing function parameters.  */
  sk_class,	     /* The scope containing the members of a class.  */
  sk_scoped_enum,    /* The scope containing the enumerators of a C++11
                        scoped enumeration.  */
  sk_namespace,	     /* The scope containing the members of a
			namespace, including the global scope.  */
  sk_template_parms, /* A scope for template parameters.  */
  sk_template_spec,  /* Like sk_template_parms, but for an explicit
			specialization.  Since, by definition, an
			explicit specialization is introduced by
			"template <>", this scope is always empty.  */
  sk_transaction,    /* A synchronized or atomic statement.  */
  sk_omp	     /* An OpenMP structured block.  */
};

struct GTY(()) cp_class_binding {
  cxx_binding *base;
  /* The bound name.  */
  tree identifier;
};

/* For each binding contour we allocate a binding_level structure
   which records the names defined in that contour.
   Contours include:
    0) the global one
    1) one for each function definition,
       where internal declarations of the parameters appear.
    2) one for each compound statement,
       to record its declarations.

   The current meaning of a name can be found by searching the levels
   from the current one out to the global one.

   Off to the side, may be the class_binding_level.  This exists only
   to catch class-local declarations.  It is otherwise nonexistent.

   Also there may be binding levels that catch cleanups that must be
   run when exceptions occur.  Thus, to see whether a name is bound in
   the current scope, it is not enough to look in the
   CURRENT_BINDING_LEVEL.  You should use lookup_name_current_level
   instead.  */

struct GTY(()) cp_binding_level {
  /* A chain of _DECL nodes for all variables, constants, functions,
      and typedef types.  These are in the reverse of the order
      supplied.  There may be OVERLOADs on this list, too, but they
      are wrapped in TREE_LISTs; the TREE_VALUE is the OVERLOAD.  */
  tree names;

  /* Using directives.  */
  vec<tree, va_gc> *using_directives;

  /* For the binding level corresponding to a class, the entities
      declared in the class or its base classes.  */
  vec<cp_class_binding, va_gc> *class_shadowed;

  /* Similar to class_shadowed, but for IDENTIFIER_TYPE_VALUE, and
      is used for all binding levels. The TREE_PURPOSE is the name of
      the entity, the TREE_TYPE is the associated type.  In addition
      the TREE_VALUE is the IDENTIFIER_TYPE_VALUE before we entered
      the class.  */
  tree type_shadowed;

  /* For each level (except not the global one),
      a chain of BLOCK nodes for all the levels
      that were entered and exited one level down.  */
  tree blocks;

  /* The entity (namespace, class, function) the scope of which this
      binding contour corresponds to.  Otherwise NULL.  */
  tree this_entity;

  /* The binding level which this one is contained in (inherits from).  */
  cp_binding_level *level_chain;

  /* STATEMENT_LIST for statements in this binding contour.
      Only used at present for SK_CLEANUP temporary bindings.  */
  tree statement_list;

  /* Binding depth at which this level began.  */
  int binding_depth;

  /* The kind of scope that this object represents.  However, a
      SK_TEMPLATE_SPEC scope is represented with KIND set to
      SK_TEMPLATE_PARMS and EXPLICIT_SPEC_P set to true.  */
  ENUM_BITFIELD (scope_kind) kind : 4;

  /* True if this scope is an SK_TEMPLATE_SPEC scope.  This field is
      only valid if KIND == SK_TEMPLATE_PARMS.  */
  BOOL_BITFIELD explicit_spec_p : 1;

  /* true means make a BLOCK for this level regardless of all else.  */
  unsigned keep : 1;

  /* Nonzero if this level can safely have additional
      cleanup-needing variables added to it.  */
  unsigned more_cleanups_ok : 1;
  unsigned have_cleanups : 1;

  /* Transient state set if this scope is of sk_class kind
     and is in the process of defining 'this_entity'.  Reset
     on leaving the class definition to allow for the scope
     to be subsequently re-used as a non-defining scope for
     'this_entity'.  */
  unsigned defining_class_p : 1;

  /* true for SK_FUNCTION_PARMS of immediate functions.  */
  unsigned immediate_fn_ctx_p : 1;

  /* True for SK_FUNCTION_PARMS of a requires-expression.  */
  unsigned requires_expression: 1;

  /* 21 bits left to fill a 32-bit word.  */
};

/* The binding level currently in effect.  */

#define current_binding_level			\
  (*(cfun && cp_function_chain && cp_function_chain->bindings \
   ? &cp_function_chain->bindings		\
   : &scope_chain->bindings))

/* The binding level of the current class, if any.  */

#define class_binding_level scope_chain->class_bindings

/* True if SCOPE designates the global scope binding contour.  */
#define global_scope_p(SCOPE) \
  ((SCOPE) == NAMESPACE_LEVEL (global_namespace))

extern cp_binding_level *leave_scope (void);
extern bool kept_level_p (void);
extern bool global_bindings_p (void);
extern bool toplevel_bindings_p (void);
extern bool namespace_bindings_p (void);
extern bool local_bindings_p (void);
extern bool template_parm_scope_p (void);
extern scope_kind innermost_scope_kind (void);
extern cp_binding_level *begin_scope (scope_kind, tree);
extern void print_binding_stack	(void);
extern void pop_everything (void);
extern void keep_next_level (bool);
extern bool is_ancestor (tree ancestor, tree descendant);
extern bool is_nested_namespace (tree parent, tree descendant,
				 bool inline_only = false);
extern tree push_scope (tree);
extern void pop_scope (tree);
extern tree push_inner_scope (tree);
extern void pop_inner_scope (tree, tree);
extern void push_binding_level (cp_binding_level *);

extern bool handle_namespace_attrs (tree, tree);
extern void pushlevel_class (void);
extern void poplevel_class (void);

/* What kind of scopes name lookup looks in.  An enum class so we
   don't accidentally mix integers.  */
enum class LOOK_where
{
  BLOCK = 1 << 0,  /* Consider block scopes.  */ 
  CLASS = 1 << 1,  /* Consider class scopes.  */ 
  NAMESPACE = 1 << 2,  /* Consider namespace scopes.  */ 

  ALL = BLOCK | CLASS | NAMESPACE,
  BLOCK_NAMESPACE = BLOCK | NAMESPACE,
  CLASS_NAMESPACE = CLASS | NAMESPACE,
};
constexpr LOOK_where operator| (LOOK_where a, LOOK_where b)
{
  return LOOK_where (unsigned (a) | unsigned (b));
}
constexpr LOOK_where operator& (LOOK_where a, LOOK_where b)
{
  return LOOK_where (unsigned (a) & unsigned (b));
}

enum class LOOK_want
{
  NORMAL = 0,  /* Normal lookup -- non-types can hide implicit types.  */
  TYPE = 1 << 1,  /* We only want TYPE_DECLS.  */
  NAMESPACE = 1 << 2,  /* We only want NAMESPACE_DECLS.  */

  HIDDEN_FRIEND = 1 << 3, /* See hidden friends.  */
  HIDDEN_LAMBDA = 1 << 4,  /* See lambda-ignored entities.  */

  TYPE_NAMESPACE = TYPE | NAMESPACE,  /* Either NAMESPACE or TYPE.  */
};
constexpr LOOK_want operator| (LOOK_want a, LOOK_want b)
{
  return LOOK_want (unsigned (a) | unsigned (b));
}
constexpr LOOK_want operator& (LOOK_want a, LOOK_want b)
{
  return LOOK_want (unsigned (a) & unsigned (b));
}

extern tree lookup_name (tree, LOOK_where, LOOK_want = LOOK_want::NORMAL);
/* Also declared in c-family/c-common.h.  */
extern tree lookup_name (tree name);
inline tree lookup_name (tree name, LOOK_want want)
{
  return lookup_name (name, LOOK_where::ALL, want);
}

enum class TAG_how
{
  CURRENT_ONLY = 0, // Look and insert only in current scope

  GLOBAL = 1, // Unqualified lookup, innermost-non-class insertion

  INNERMOST_NON_CLASS = 2, // Look and insert only into
			   // innermost-non-class

  HIDDEN_FRIEND = 3, // As INNERMOST_NON_CLASS, but hide it
};

extern tree lookup_elaborated_type (tree, TAG_how);
extern tree get_namespace_binding (tree ns, tree id);
extern void set_global_binding (tree decl);
inline tree get_global_binding (tree id)
{
  return get_namespace_binding (NULL_TREE, id);
}
extern tree lookup_qualified_name (tree scope, tree name,
				   LOOK_want = LOOK_want::NORMAL,
				   bool = true);
extern tree lookup_qualified_name (tree scope, const char *name,
				   LOOK_want = LOOK_want::NORMAL,
				   bool = true);
extern bool pushdecl_class_level (tree);
extern tree pushdecl_namespace_level (tree, bool hiding = false);
extern bool push_class_level_binding (tree, tree);
extern tree get_local_decls ();
extern int function_parm_depth (void);
extern tree cp_namespace_decls (tree);
extern void set_decl_namespace (tree, tree, bool);
extern void push_decl_namespace (tree);
extern void pop_decl_namespace (void);
extern void do_namespace_alias (tree, tree);
extern tree do_class_using_decl (tree, tree);
extern tree lookup_arg_dependent (tree, tree, vec<tree, va_gc> *);
extern tree search_anon_aggr (tree, tree, bool = false);
extern tree get_class_binding_direct (tree, tree, bool want_type = false);
extern tree get_class_binding (tree, tree, bool want_type = false);
extern tree *find_member_slot (tree klass, tree name);
extern tree *add_member_slot (tree klass, tree name);
extern void resort_type_member_vec (void *, void *,
				    gt_pointer_operator, void *);
extern vec<tree, va_gc> *set_class_bindings (tree, int extra = 0);
extern void insert_late_enum_def_bindings (tree, tree);
extern tree innermost_non_namespace_value (tree);
extern cxx_binding *outer_binding (tree, cxx_binding *, bool);
extern void cp_emit_debug_info_for_using (tree, tree);

extern void finish_nonmember_using_decl (tree scope, tree name);
extern void finish_using_directive (tree target, tree attribs);
void push_local_extern_decl_alias (tree decl);
extern tree pushdecl (tree, bool hiding = false);
extern tree pushdecl_outermost_localscope (tree);
extern tree pushdecl_top_level (tree);
extern tree pushdecl_top_level_and_finish (tree, tree);
extern tree pushtag (tree, tree, TAG_how = TAG_how::CURRENT_ONLY);
extern int push_namespace (tree, bool make_inline = false);
extern void pop_namespace (void);
extern void push_nested_namespace (tree);
extern void pop_nested_namespace (tree);
extern void push_to_top_level (void);
extern void pop_from_top_level (void);
extern void push_using_decl_bindings (tree, tree);

/* Lower level interface for modules. */
extern tree *mergeable_namespace_slots (tree ns, tree name, bool is_global,
					tree *mvec);
extern void add_mergeable_namespace_entity (tree *slot, tree decl);
extern tree lookup_class_binding (tree ctx, tree name);
extern bool import_module_binding (tree ctx, tree name, unsigned mod,
				   unsigned snum);
extern bool set_module_binding (tree ctx, tree name, unsigned mod,
				int mod_glob_flag,
				tree value, tree type, tree visible);
extern void add_module_namespace_decl (tree ns, tree decl);

enum WMB_Flags
{
  WMB_None = 0,
  WMB_Dups = 1 << 0,
  WMB_Export = 1 << 1,
  WMB_Using = 1 << 2,
  WMB_Hidden = 1 << 3,
};

extern unsigned walk_module_binding (tree binding, bitmap partitions,
				     bool (*)(tree decl, WMB_Flags, void *data),
				     void *data);
extern tree add_imported_namespace (tree ctx, tree name, location_t,
				    unsigned module,
				    bool inline_p, bool visible_p);
extern const char *get_cxx_dialect_name (enum cxx_dialect dialect);

#endif /* GCC_CP_NAME_LOOKUP_H */
