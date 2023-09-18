/* Interprocedural analyses.
   Copyright (C) 2005-2022 Free Software Foundation, Inc.

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

#ifndef IPA_PROP_H
#define IPA_PROP_H

/* The following definitions and interfaces are used by
   interprocedural analyses or parameters.  */

#define IPA_UNDESCRIBED_USE -1

/* ipa-prop.cc stuff (ipa-cp, indirect inlining):  */

/* A jump function for a callsite represents the values passed as actual
   arguments of the callsite.  They were originally proposed in a paper called
   "Interprocedural Constant Propagation", by David Callahan, Keith D Cooper,
   Ken Kennedy, Linda Torczon in Comp86, pg 152-161.  There are three main
   types of values :

   Pass-through - the caller's formal parameter is passed as an actual
                  argument, possibly one simple operation performed on it.
   Constant     - a constant (is_gimple_ip_invariant)is passed as an actual
                  argument.
   Unknown      - neither of the above.

   IPA_JF_LOAD_AGG is a compound pass-through jump function, in which primary
   operation on formal parameter is memory dereference that loads a value from
   a part of an aggregate, which is represented or pointed to by the formal
   parameter.  Moreover, an additional unary/binary operation can be applied on
   the loaded value, and final result is passed as actual argument of callee
   (e.g. *(param_1(D) + 4) op 24 ).  It is meant to describe usage of aggregate
   parameter or by-reference parameter referenced in argument passing, commonly
   found in C++ and Fortran.

   IPA_JF_ANCESTOR is a special pass-through jump function, which means that
   the result is an address of a part of the object pointed to by the formal
   parameter to which the function refers.  It is mainly intended to represent
   getting addresses of ancestor fields in C++
   (e.g. &this_1(D)->D.1766.D.1756).  Note that if the original pointer is
   NULL, ancestor jump function must behave like a simple pass-through.

   Other pass-through functions can either simply pass on an unchanged formal
   parameter or can apply one simple binary operation to it (such jump
   functions are called polynomial).

   Jump functions are computed in ipa-prop.cc by function
   update_call_notes_after_inlining.  Some information can be lost and jump
   functions degraded accordingly when inlining, see
   update_call_notes_after_inlining in the same file.  */

enum jump_func_type
{
  IPA_JF_UNKNOWN = 0,  /* newly allocated and zeroed jump functions default */
  IPA_JF_CONST,             /* represented by field costant */
  IPA_JF_PASS_THROUGH,	    /* represented by field pass_through */
  IPA_JF_LOAD_AGG,	    /* represented by field load_agg */
  IPA_JF_ANCESTOR	    /* represented by field ancestor */
};

struct ipa_cst_ref_desc;

/* Structure holding data required to describe a constant jump function.  */
struct GTY(()) ipa_constant_data
{
  /* THe value of the constant.  */
  tree value;
  /* Pointer to the structure that describes the reference.  */
  struct ipa_cst_ref_desc GTY((skip)) *rdesc;
};

/* Structure holding data required to describe a pass-through jump function.  */

struct GTY(()) ipa_pass_through_data
{
  /* If an operation is to be performed on the original parameter, this is the
     second (constant) operand.  */
  tree operand;
  /* Number of the caller's formal parameter being passed.  */
  int formal_id;
  /* Operation that is performed on the argument before it is passed on.
     Special values which have other meaning than in normal contexts:
       - NOP_EXPR means no operation, not even type conversion.
       - ASSERT_EXPR means that only the value in operand is allowed to pass
         through (without any change), for all other values the result is
         unknown.
     Otherwise operation must be a simple binary or unary arithmetic operation
     where the caller's parameter is the first operand and (for binary
     operations) the operand field from this structure is the second one.  */
  enum tree_code operation;
  /* When the passed value is a pointer, it is set to true only when we are
     certain that no write to the object it points to has occurred since the
     caller functions started execution, except for changes noted in the
     aggregate part of the jump function (see description of
     ipa_agg_jump_function).  The flag is used only when the operation is
     NOP_EXPR.  */
  unsigned agg_preserved : 1;
  /* Set when the edge has already been used to decrement an appropriate
     reference description counter and should not be decremented again.  */
  unsigned refdesc_decremented : 1;
};

/* Structure holding data required to describe a load-value-from-aggregate
   jump function.  */

struct GTY(()) ipa_load_agg_data
{
  /* Inherit from pass through jump function, describing unary/binary
     operation on the value loaded from aggregate that is represented or
     pointed to by the formal parameter, specified by formal_id in this
     pass_through jump function data structure.  */
  struct ipa_pass_through_data pass_through;
  /* Type of the value loaded from the aggregate.  */
  tree type;
  /* Offset at which the value is located within the aggregate.  */
  HOST_WIDE_INT offset;
  /* True if loaded by reference (the aggregate is pointed to by the formal
     parameter) or false if loaded by value (the aggregate is represented
     by the formal parameter).  */
  bool by_ref;
};

/* Structure holding data required to describe an ancestor pass-through
   jump function.  */

struct GTY(()) ipa_ancestor_jf_data
{
  /* Offset of the field representing the ancestor.  */
  HOST_WIDE_INT offset;
  /* Number of the caller's formal parameter being passed.  */
  int formal_id;
  /* Flag with the same meaning like agg_preserve in ipa_pass_through_data.  */
  unsigned agg_preserved : 1;
  /* When set, the operation should not have any effect on NULL pointers.  */
  unsigned keep_null : 1;
};

/* A jump function for an aggregate part at a given offset, which describes how
   it content value is generated.  All unlisted positions are assumed to have a
   value defined in an unknown way.  */

struct GTY(()) ipa_agg_jf_item
{
  /* The offset for the aggregate part.  */
  HOST_WIDE_INT offset;

  /* Data type of the aggregate part.  */
  tree type;

  /* Jump function type.  */
  enum jump_func_type jftype;

  /* Represents a value of jump function. constant represents the actual constant
     in constant jump function content.  pass_through is used only in simple pass
     through jump function context.  load_agg is for load-value-from-aggregate
     jump function context.  */
  union jump_func_agg_value
  {
    tree GTY ((tag ("IPA_JF_CONST"))) constant;
    struct ipa_pass_through_data GTY ((tag ("IPA_JF_PASS_THROUGH"))) pass_through;
    struct ipa_load_agg_data GTY ((tag ("IPA_JF_LOAD_AGG"))) load_agg;
  } GTY ((desc ("%1.jftype"))) value;
};

/* Jump functions describing a set of aggregate contents.  */

struct GTY(()) ipa_agg_jump_function
{
  /* Description of the individual jump function item.  */
  vec<ipa_agg_jf_item, va_gc> *items;
  /* True if the data was passed by reference (as opposed to by value).  */
  bool by_ref;
};

/* An element in an aggregate part describing a known value at a given offset.
   All unlisted positions are assumed to be unknown and all listed values must
   fulfill is_gimple_ip_invariant.  */

struct ipa_agg_value
{
  /* The offset at which the known value is located within the aggregate.  */
  HOST_WIDE_INT offset;

  /* The known constant.  */
  tree value;

  /* Return true if OTHER describes same agg value.  */
  bool equal_to (const ipa_agg_value &other);
};

/* Structure describing a set of known offset/value for aggregate.  */

struct ipa_agg_value_set
{
  /* Description of the individual item.  */
  vec<ipa_agg_value> items;
  /* True if the data was passed by reference (as opposed to by value).  */
  bool by_ref;

  /* Return true if OTHER describes same agg values.  */
  bool equal_to (const ipa_agg_value_set &other)
  {
    if (by_ref != other.by_ref)
      return false;
    if (items.length () != other.items.length ())
      return false;
    for (unsigned int i = 0; i < items.length (); i++)
      if (!items[i].equal_to (other.items[i]))
	return false;
    return true;
  }

  /* Return true if there is any value for aggregate.  */
  bool is_empty () const
  {
    return items.is_empty ();
  }

  ipa_agg_value_set copy () const
  {
    ipa_agg_value_set new_copy;

    new_copy.items = items.copy ();
    new_copy.by_ref = by_ref;

    return new_copy;
  }

  void release ()
  {
    items.release ();
  }
};

/* Return copy of a vec<ipa_agg_value_set>.  */

static inline vec<ipa_agg_value_set>
ipa_copy_agg_values (const vec<ipa_agg_value_set> &aggs)
{
  vec<ipa_agg_value_set> aggs_copy = vNULL;

  if (!aggs.is_empty ())
    {
      ipa_agg_value_set *agg;
      int i;

      aggs_copy.reserve_exact (aggs.length ());

      FOR_EACH_VEC_ELT (aggs, i, agg)
	aggs_copy.quick_push (agg->copy ());
    }

  return aggs_copy;
}

/* For vec<ipa_agg_value_set>, DO NOT call release(), use below function
   instead.  Because ipa_agg_value_set contains a field of vector type, we
   should release this child vector in each element before reclaiming the
   whole vector.  */

static inline void
ipa_release_agg_values (vec<ipa_agg_value_set> &aggs,
			bool release_vector = true)
{
  ipa_agg_value_set *agg;
  int i;

  FOR_EACH_VEC_ELT (aggs, i, agg)
    agg->release ();
  if (release_vector)
    aggs.release ();
}

/* Information about zero/non-zero bits.  */
class GTY(()) ipa_bits
{
public:
  /* The propagated value.  */
  widest_int value;
  /* Mask corresponding to the value.
     Similar to ccp_lattice_t, if xth bit of mask is 0,
     implies xth bit of value is constant.  */
  widest_int mask;
};

/* Info about value ranges.  */

class GTY(()) ipa_vr
{
public:
  /* The data fields below are valid only if known is true.  */
  bool known;
  enum value_range_kind type;
  wide_int min;
  wide_int max;
  bool nonzero_p (tree) const;
};

/* A jump function for a callsite represents the values passed as actual
   arguments of the callsite. See enum jump_func_type for the various
   types of jump functions supported.  */
struct GTY (()) ipa_jump_func
{
  /* Aggregate jump function description.  See struct ipa_agg_jump_function
     and its description.  */
  struct ipa_agg_jump_function agg;

  /* Information about zero/non-zero bits.  The pointed to structure is shared
     betweed different jump functions.  Use ipa_set_jfunc_bits to set this
     field.  */
  class ipa_bits *bits;

  /* Information about value range, containing valid data only when vr_known is
     true.  The pointed to structure is shared betweed different jump
     functions.  Use ipa_set_jfunc_vr to set this field.  */
  value_range *m_vr;

  enum jump_func_type type;
  /* Represents a value of a jump function.  pass_through is used only in jump
     function context.  constant represents the actual constant in constant jump
     functions and member_cst holds constant c++ member functions.  */
  union jump_func_value
  {
    struct ipa_constant_data GTY ((tag ("IPA_JF_CONST"))) constant;
    struct ipa_pass_through_data GTY ((tag ("IPA_JF_PASS_THROUGH"))) pass_through;
    struct ipa_ancestor_jf_data GTY ((tag ("IPA_JF_ANCESTOR"))) ancestor;
  } GTY ((desc ("%1.type"))) value;
};


/* Return the constant stored in a constant jump functin JFUNC.  */

static inline tree
ipa_get_jf_constant (struct ipa_jump_func *jfunc)
{
  gcc_checking_assert (jfunc->type == IPA_JF_CONST);
  return jfunc->value.constant.value;
}

static inline struct ipa_cst_ref_desc *
ipa_get_jf_constant_rdesc (struct ipa_jump_func *jfunc)
{
  gcc_checking_assert (jfunc->type == IPA_JF_CONST);
  return jfunc->value.constant.rdesc;
}

/* Make JFUNC not participate in any further reference counting.  */

inline void
ipa_zap_jf_refdesc (ipa_jump_func *jfunc)
{
  gcc_checking_assert (jfunc->type == IPA_JF_CONST);
  jfunc->value.constant.rdesc = NULL;
}

/* Return the operand of a pass through jmp function JFUNC.  */

static inline tree
ipa_get_jf_pass_through_operand (struct ipa_jump_func *jfunc)
{
  gcc_checking_assert (jfunc->type == IPA_JF_PASS_THROUGH);
  return jfunc->value.pass_through.operand;
}

/* Return the number of the caller's formal parameter that a pass through jump
   function JFUNC refers to.  */

static inline int
ipa_get_jf_pass_through_formal_id (struct ipa_jump_func *jfunc)
{
  gcc_checking_assert (jfunc->type == IPA_JF_PASS_THROUGH);
  return jfunc->value.pass_through.formal_id;
}

/* Return operation of a pass through jump function JFUNC.  */

static inline enum tree_code
ipa_get_jf_pass_through_operation (struct ipa_jump_func *jfunc)
{
  gcc_checking_assert (jfunc->type == IPA_JF_PASS_THROUGH);
  return jfunc->value.pass_through.operation;
}

/* Return the agg_preserved flag of a pass through jump function JFUNC.  */

static inline bool
ipa_get_jf_pass_through_agg_preserved (struct ipa_jump_func *jfunc)
{
  gcc_checking_assert (jfunc->type == IPA_JF_PASS_THROUGH);
  return jfunc->value.pass_through.agg_preserved;
}

/* Return the refdesc_decremented flag of a pass through jump function
   JFUNC.  */

inline bool
ipa_get_jf_pass_through_refdesc_decremented (struct ipa_jump_func *jfunc)
{
  gcc_checking_assert (jfunc->type == IPA_JF_PASS_THROUGH);
  return jfunc->value.pass_through.refdesc_decremented;
}

/* Set the refdesc_decremented flag of a pass through jump function JFUNC to
   VALUE.  */

inline void
ipa_set_jf_pass_through_refdesc_decremented (ipa_jump_func *jfunc, bool value)
{
  gcc_checking_assert (jfunc->type == IPA_JF_PASS_THROUGH);
  jfunc->value.pass_through.refdesc_decremented = value;
}

/* Return true if pass through jump function JFUNC preserves type
   information.  */

static inline bool
ipa_get_jf_pass_through_type_preserved (struct ipa_jump_func *jfunc)
{
  gcc_checking_assert (jfunc->type == IPA_JF_PASS_THROUGH);
  return jfunc->value.pass_through.agg_preserved;
}

/* Return the offset of an ancestor jump function JFUNC.  */

static inline HOST_WIDE_INT
ipa_get_jf_ancestor_offset (struct ipa_jump_func *jfunc)
{
  gcc_checking_assert (jfunc->type == IPA_JF_ANCESTOR);
  return jfunc->value.ancestor.offset;
}

/* Return the number of the caller's formal parameter that an ancestor jump
   function JFUNC refers to.  */

static inline int
ipa_get_jf_ancestor_formal_id (struct ipa_jump_func *jfunc)
{
  gcc_checking_assert (jfunc->type == IPA_JF_ANCESTOR);
  return jfunc->value.ancestor.formal_id;
}

/* Return the agg_preserved flag of an ancestor jump function JFUNC.  */

static inline bool
ipa_get_jf_ancestor_agg_preserved (struct ipa_jump_func *jfunc)
{
  gcc_checking_assert (jfunc->type == IPA_JF_ANCESTOR);
  return jfunc->value.ancestor.agg_preserved;
}

/* Return true if ancestor jump function JFUNC presrves type information.  */

static inline bool
ipa_get_jf_ancestor_type_preserved (struct ipa_jump_func *jfunc)
{
  gcc_checking_assert (jfunc->type == IPA_JF_ANCESTOR);
  return jfunc->value.ancestor.agg_preserved;
}

/* Return if jfunc represents an operation whether we first check the formal
   parameter for non-NULLness unless it does not matter because the offset is
   zero anyway.  */

static inline bool
ipa_get_jf_ancestor_keep_null (struct ipa_jump_func *jfunc)
{
  gcc_checking_assert (jfunc->type == IPA_JF_ANCESTOR);
  return jfunc->value.ancestor.keep_null;
}

/* Class for allocating a bundle of various potentially known properties about
   actual arguments of a particular call on stack for the usual case and on
   heap only if there are unusually many arguments.  The data is deallocated
   when the instance of this class goes out of scope or is otherwise
   destructed.  */

class ipa_auto_call_arg_values
{
public:
  ~ipa_auto_call_arg_values ();

  /* If m_known_vals (vector of known "scalar" values) is sufficiantly long,
     return its element at INDEX, otherwise return NULL.  */
  tree safe_sval_at (int index)
  {
    /* TODO: Assert non-negative index here and test.  */
    if ((unsigned) index < m_known_vals.length ())
      return m_known_vals[index];
    return NULL;
  }

  /* If m_known_aggs is sufficiantly long, return the pointer rto its element
     at INDEX, otherwise return NULL.  */
  ipa_agg_value_set *safe_aggval_at (int index)
  {
    /* TODO: Assert non-negative index here and test.  */
    if ((unsigned) index < m_known_aggs.length ())
      return &m_known_aggs[index];
    return NULL;
  }

  /* Vector describing known values of parameters.  */
  auto_vec<tree, 32> m_known_vals;

  /* Vector describing known polymorphic call contexts.  */
  auto_vec<ipa_polymorphic_call_context, 32> m_known_contexts;

  /* Vector describing known aggregate values.  */
  auto_vec<ipa_agg_value_set, 32> m_known_aggs;

  /* Vector describing known value ranges of arguments.  */
  auto_vec<value_range, 32> m_known_value_ranges;
};

/* Class bundling the various potentially known properties about actual
   arguments of a particular call.  This variant does not deallocate the
   bundled data in any way.  */

class ipa_call_arg_values
{
public:
  /* Default constructor, setting the vectors to empty ones.  */
  ipa_call_arg_values ()
  {}

  /* Construct this general variant of the bundle from the variant which uses
     auto_vecs to hold the vectors.  This means that vectors of objects
     constructed with this constructor should not be changed because if they
     get reallocated, the member vectors and the underlying auto_vecs would get
     out of sync.  */
  ipa_call_arg_values (ipa_auto_call_arg_values *aavals)
    : m_known_vals (aavals->m_known_vals.to_vec_legacy ()),
      m_known_contexts (aavals->m_known_contexts.to_vec_legacy ()),
      m_known_aggs (aavals->m_known_aggs.to_vec_legacy ()),
      m_known_value_ranges (aavals->m_known_value_ranges.to_vec_legacy ())
  {}

  /* If m_known_vals (vector of known "scalar" values) is sufficiantly long,
     return its element at INDEX, otherwise return NULL.  */
  tree safe_sval_at (int index)
  {
    /* TODO: Assert non-negative index here and test.  */
    if ((unsigned) index < m_known_vals.length ())
      return m_known_vals[index];
    return NULL;
  }

  /* If m_known_aggs is sufficiantly long, return the pointer rto its element
     at INDEX, otherwise return NULL.  */
  ipa_agg_value_set *safe_aggval_at (int index)
  {
    /* TODO: Assert non-negative index here and test.  */
    if ((unsigned) index < m_known_aggs.length ())
      return &m_known_aggs[index];
    return NULL;
  }

  /* Vector describing known values of parameters.  */
  vec<tree> m_known_vals = vNULL;

  /* Vector describing known polymorphic call contexts.  */
  vec<ipa_polymorphic_call_context> m_known_contexts = vNULL;

  /* Vector describing known aggregate values.  */
  vec<ipa_agg_value_set> m_known_aggs = vNULL;

  /* Vector describing known value ranges of arguments.  */
  vec<value_range> m_known_value_ranges = vNULL;
};


/* Summary describing a single formal parameter.  */

struct GTY(()) ipa_param_descriptor
{
  /* In analysis and modification phase, this is the PARAM_DECL of this
     parameter, in IPA LTO phase, this is the type of the described
     parameter or NULL if not known.  Do not read this field directly but
     through ipa_get_param and ipa_get_type as appropriate.  */
  tree decl_or_type;
  /* If all uses of the parameter are described by ipa-prop structures, this
     says how many there are.  If any use could not be described by means of
     ipa-prop structures (which include flag dereferenced below), this is
     IPA_UNDESCRIBED_USE.  */
  int controlled_uses;
  unsigned int move_cost : 27;
  /* The parameter is used.  */
  unsigned used : 1;
  unsigned used_by_ipa_predicates : 1;
  unsigned used_by_indirect_call : 1;
  unsigned used_by_polymorphic_call : 1;
  /* Set to true when in addition to being used in call statements, the
     parameter has also been used for loads (but not for writes, does not
     escape, etc.).  This allows us to identify parameters p which are only
     used as *p, and so when we propagate a constant to them, we can generate a
     LOAD and not ADDR reference to them.  */
  unsigned load_dereferenced : 1;
};

/* ipa_node_params stores information related to formal parameters of functions
   and some other information for interprocedural passes that operate on
   parameters (such as ipa-cp).  */

class GTY((for_user)) ipa_node_params
{
public:
  /* Default constructor.  */
  ipa_node_params ();

  /* Default destructor.  */
  ~ipa_node_params ();

  /* Information about individual formal parameters that are gathered when
     summaries are generated. */
  vec<ipa_param_descriptor, va_gc> *descriptors;
  /* Pointer to an array of structures describing individual formal
     parameters.  */
  class ipcp_param_lattices * GTY((skip)) lattices;
  /* Only for versioned nodes this field would not be NULL,
     it points to the node that IPA cp cloned from.  */
  struct cgraph_node * GTY((skip)) ipcp_orig_node;
  /* If this node is an ipa-cp clone, these are the known constants that
     describe what it has been specialized for.  */
  vec<tree> GTY((skip)) known_csts;
  /* If this node is an ipa-cp clone, these are the known polymorphic contexts
     that describe what it has been specialized for.  */
  vec<ipa_polymorphic_call_context> GTY((skip)) known_contexts;
  /* Whether the param uses analysis and jump function computation has already
     been performed.  */
  unsigned analysis_done : 1;
  /* Whether the function is enqueued in ipa-cp propagation stack.  */
  unsigned node_enqueued : 1;
  /* Whether we should create a specialized version based on values that are
     known to be constant in all contexts.  */
  unsigned do_clone_for_all_contexts : 1;
  /* Set if this is an IPA-CP clone for all contexts.  */
  unsigned is_all_contexts_clone : 1;
  /* Node has been completely replaced by clones and will be removed after
     ipa-cp is finished.  */
  unsigned node_dead : 1;
  /* Node is involved in a recursion, potentionally indirect.  */
  unsigned node_within_scc : 1;
  /* Node contains only direct recursion.  */
  unsigned node_is_self_scc : 1;
  /* Node is calling a private function called only once.  */
  unsigned node_calling_single_call : 1;
  /* False when there is something makes versioning impossible.  */
  unsigned versionable : 1;
};

inline
ipa_node_params::ipa_node_params ()
: descriptors (NULL), lattices (NULL), ipcp_orig_node (NULL),
  known_csts (vNULL), known_contexts (vNULL), analysis_done (0),
  node_enqueued (0), do_clone_for_all_contexts (0), is_all_contexts_clone (0),
  node_dead (0), node_within_scc (0), node_is_self_scc (0),
  node_calling_single_call (0), versionable (0)
{
}

inline
ipa_node_params::~ipa_node_params ()
{
  free (lattices);
  vec_free (descriptors);
  known_csts.release ();
  known_contexts.release ();
}

/* Intermediate information that we get from alias analysis about a particular
   parameter in a particular basic_block.  When a parameter or the memory it
   references is marked modified, we use that information in all dominated
   blocks without consulting alias analysis oracle.  */

struct ipa_param_aa_status
{
  /* Set when this structure contains meaningful information.  If not, the
     structure describing a dominating BB should be used instead.  */
  bool valid;

  /* Whether we have seen something which might have modified the data in
     question.  PARM is for the parameter itself, REF is for data it points to
     but using the alias type of individual accesses and PT is the same thing
     but for computing aggregate pass-through functions using a very inclusive
     ao_ref.  */
  bool parm_modified, ref_modified, pt_modified;
};

/* Information related to a given BB that used only when looking at function
   body.  */

struct ipa_bb_info
{
  /* Call graph edges going out of this BB.  */
  vec<cgraph_edge *> cg_edges;
  /* Alias analysis statuses of each formal parameter at this bb.  */
  vec<ipa_param_aa_status> param_aa_statuses;
};

/* Structure with global information that is only used when looking at function
   body. */

struct ipa_func_body_info
{
  /* The node that is being analyzed.  */
  cgraph_node *node;

  /* Its info.  */
  class ipa_node_params *info;

  /* Information about individual BBs. */
  vec<ipa_bb_info> bb_infos;

  /* Number of parameters.  */
  int param_count;

  /* Number of statements we are still allowed to walked by when analyzing this
     function.  */
  unsigned int aa_walk_budget;
};

/* ipa_node_params access functions.  Please use these to access fields that
   are or will be shared among various passes.  */

/* Return the number of formal parameters. */

static inline int
ipa_get_param_count (class ipa_node_params *info)
{
  return vec_safe_length (info->descriptors);
}

/* Return the parameter declaration in DESCRIPTORS at index I and assert it is
   indeed a PARM_DECL.  */

static inline tree
ipa_get_param (const vec<ipa_param_descriptor, va_gc> &descriptors, int i)
{
  tree t = descriptors[i].decl_or_type;
  gcc_checking_assert (TREE_CODE (t) == PARM_DECL);
  return t;
}

/* Return the declaration of Ith formal parameter of the function corresponding
   to INFO.  Note there is no setter function as this array is built just once
   using ipa_initialize_node_params.  This function should not be called in
   WPA.  */

static inline tree
ipa_get_param (class ipa_node_params *info, int i)
{
  gcc_checking_assert (info->descriptors);
  return ipa_get_param (*info->descriptors, i);
}

/* Return the type of Ith formal parameter of the function corresponding
   to INFO if it is known or NULL if not.  */

static inline tree
ipa_get_type (class ipa_node_params *info, int i)
{
  if (vec_safe_length (info->descriptors) <= (unsigned) i)
    return NULL;
  tree t = (*info->descriptors)[i].decl_or_type;
  if (!t)
    return NULL;
  if (TYPE_P (t))
    return t;
  gcc_checking_assert (TREE_CODE (t) == PARM_DECL);
  return TREE_TYPE (t);
}

/* Return the move cost of Ith formal parameter of the function corresponding
   to INFO.  */

static inline int
ipa_get_param_move_cost (class ipa_node_params *info, int i)
{
  gcc_checking_assert (info->descriptors);
  return (*info->descriptors)[i].move_cost;
}

/* Set the used flag corresponding to the Ith formal parameter of the function
   associated with INFO to VAL.  */

static inline void
ipa_set_param_used (class ipa_node_params *info, int i, bool val)
{
  gcc_checking_assert (info->descriptors);
  (*info->descriptors)[i].used = val;
}

/* Set the used_by_ipa_predicates flag corresponding to the Ith formal
   parameter of the function associated with INFO to VAL.  */

static inline void
ipa_set_param_used_by_ipa_predicates (class ipa_node_params *info, int i, bool val)
{
  gcc_checking_assert (info->descriptors);
  (*info->descriptors)[i].used_by_ipa_predicates = val;
}

/* Set the used_by_indirect_call flag corresponding to the Ith formal
   parameter of the function associated with INFO to VAL.  */

static inline void
ipa_set_param_used_by_indirect_call (class ipa_node_params *info, int i, bool val)
{
  gcc_checking_assert (info->descriptors);
  (*info->descriptors)[i].used_by_indirect_call = val;
}

/* Set the .used_by_polymorphic_call flag corresponding to the Ith formal
   parameter of the function associated with INFO to VAL.  */

static inline void
ipa_set_param_used_by_polymorphic_call (class ipa_node_params *info, int i, bool val)
{
  gcc_checking_assert (info->descriptors);
  (*info->descriptors)[i].used_by_polymorphic_call = val;
}

/* Return how many uses described by ipa-prop a parameter has or
   IPA_UNDESCRIBED_USE if there is a use that is not described by these
   structures.  */
static inline int
ipa_get_controlled_uses (class ipa_node_params *info, int i)
{
  /* FIXME: introducing speculation causes out of bounds access here.  */
  if (vec_safe_length (info->descriptors) > (unsigned)i)
    return (*info->descriptors)[i].controlled_uses;
  return IPA_UNDESCRIBED_USE;
}

/* Set the controlled counter of a given parameter.  */

static inline void
ipa_set_controlled_uses (class ipa_node_params *info, int i, int val)
{
  gcc_checking_assert (info->descriptors);
  (*info->descriptors)[i].controlled_uses = val;
}

/* Assuming a parameter does not have IPA_UNDESCRIBED_USE controlled uses,
   return flag which indicates it has been dereferenced but only in a load.  */
static inline int
ipa_get_param_load_dereferenced (class ipa_node_params *info, int i)
{
  gcc_assert (ipa_get_controlled_uses (info, i) != IPA_UNDESCRIBED_USE);
  return (*info->descriptors)[i].load_dereferenced;
}

/* Set the load_dereferenced flag of a given parameter.  */

static inline void
ipa_set_param_load_dereferenced (class ipa_node_params *info, int i, bool val)
{
  gcc_checking_assert (info->descriptors);
  (*info->descriptors)[i].load_dereferenced = val;
}

/* Return the used flag corresponding to the Ith formal parameter of the
   function associated with INFO.  */

static inline bool
ipa_is_param_used (class ipa_node_params *info, int i)
{
  gcc_checking_assert (info->descriptors);
  return (*info->descriptors)[i].used;
}

/* Return the used_by_ipa_predicates flag corresponding to the Ith formal
   parameter of the function associated with INFO.  */

static inline bool
ipa_is_param_used_by_ipa_predicates (class ipa_node_params *info, int i)
{
  gcc_checking_assert (info->descriptors);
  return (*info->descriptors)[i].used_by_ipa_predicates;
}

/* Return the used_by_indirect_call flag corresponding to the Ith formal
   parameter of the function associated with INFO.  */

static inline bool
ipa_is_param_used_by_indirect_call (class ipa_node_params *info, int i)
{
  gcc_checking_assert (info->descriptors);
  return (*info->descriptors)[i].used_by_indirect_call;
}

/* Return the used_by_polymorphic_call flag corresponding to the Ith formal
   parameter of the function associated with INFO.  */

static inline bool
ipa_is_param_used_by_polymorphic_call (class ipa_node_params *info, int i)
{
  gcc_checking_assert (info->descriptors);
  return (*info->descriptors)[i].used_by_polymorphic_call;
}

/* Information about replacements done in aggregates for a given node (each
   node has its linked list).  */
struct GTY(()) ipa_agg_replacement_value
{
  /* Next item in the linked list.  */
  struct ipa_agg_replacement_value *next;
  /* Offset within the aggregate.  */
  HOST_WIDE_INT offset;
  /* The constant value.  */
  tree value;
  /* The parameter index.  */
  int index;
  /* Whether the value was passed by reference.  */
  bool by_ref;
};

/* Structure holding information for the transformation phase of IPA-CP.  */

struct GTY(()) ipcp_transformation
{
  /* Linked list of known aggregate values.  */
  ipa_agg_replacement_value *agg_values;
  /* Known bits information.  */
  vec<ipa_bits *, va_gc> *bits;
  /* Value range information.  */
  vec<ipa_vr, va_gc> *m_vr;

  /* Default constructor.  */
  ipcp_transformation ()
  : agg_values (NULL), bits (NULL), m_vr (NULL)
  { }

  /* Default destructor.  */
  ~ipcp_transformation ()
  {
    ipa_agg_replacement_value *agg = agg_values;
    while (agg)
      {
	ipa_agg_replacement_value *next = agg->next;
	ggc_free (agg);
	agg = next;
      }
    vec_free (bits);
    vec_free (m_vr);
  }
};

void ipa_set_node_agg_value_chain (struct cgraph_node *node,
				   struct ipa_agg_replacement_value *aggvals);
void ipcp_transformation_initialize (void);
void ipcp_free_transformation_sum (void);

/* ipa_edge_args stores information related to a callsite and particularly its
   arguments.  It can be accessed by the IPA_EDGE_REF macro.  */

class GTY((for_user)) ipa_edge_args
{
 public:

  /* Default constructor.  */
  ipa_edge_args () : jump_functions (NULL), polymorphic_call_contexts (NULL)
    {}

  /* Destructor.  */
  ~ipa_edge_args ()
    {
      unsigned int i;
      ipa_jump_func *jf;
      FOR_EACH_VEC_SAFE_ELT (jump_functions, i, jf)
	vec_free (jf->agg.items);
      vec_free (jump_functions);
      vec_free (polymorphic_call_contexts);
    }

  /* Vectors of the callsite's jump function and polymorphic context
     information of each parameter.  */
  vec<ipa_jump_func, va_gc> *jump_functions;
  vec<ipa_polymorphic_call_context, va_gc> *polymorphic_call_contexts;
};

/* ipa_edge_args access functions.  Please use these to access fields that
   are or will be shared among various passes.  */

/* Return the number of actual arguments. */

static inline int
ipa_get_cs_argument_count (class ipa_edge_args *args)
{
  return vec_safe_length (args->jump_functions);
}

/* Returns a pointer to the jump function for the ith argument.  Please note
   there is no setter function as jump functions are all set up in
   ipa_compute_jump_functions. */

static inline struct ipa_jump_func *
ipa_get_ith_jump_func (class ipa_edge_args *args, int i)
{
  return &(*args->jump_functions)[i];
}

/* Returns a pointer to the polymorphic call context for the ith argument.
   NULL if contexts are not computed.  */
static inline class ipa_polymorphic_call_context *
ipa_get_ith_polymorhic_call_context (class ipa_edge_args *args, int i)
{
  if (!args->polymorphic_call_contexts)
    return NULL;
  return &(*args->polymorphic_call_contexts)[i];
}

/* Function summary for ipa_node_params.  */
class GTY((user)) ipa_node_params_t: public function_summary <ipa_node_params *>
{
public:
  ipa_node_params_t (symbol_table *table, bool ggc):
    function_summary<ipa_node_params *> (table, ggc)
  {
    disable_insertion_hook ();
  }

  /* Hook that is called by summary when a node is duplicated.  */
  virtual void duplicate (cgraph_node *node,
			  cgraph_node *node2,
			  ipa_node_params *data,
			  ipa_node_params *data2);
};

/* Summary to manange ipa_edge_args structures.  */

class GTY((user)) ipa_edge_args_sum_t : public call_summary <ipa_edge_args *>
{
 public:
  ipa_edge_args_sum_t (symbol_table *table, bool ggc)
    : call_summary<ipa_edge_args *> (table, ggc) { }

  void remove (cgraph_edge *edge)
  {
    call_summary <ipa_edge_args *>::remove (edge);
  }

  /* Hook that is called by summary when an edge is removed.  */
  virtual void remove (cgraph_edge *cs, ipa_edge_args *args);
  /* Hook that is called by summary when an edge is duplicated.  */
  virtual void duplicate (cgraph_edge *src,
			  cgraph_edge *dst,
			  ipa_edge_args *old_args,
			  ipa_edge_args *new_args);
};

/* Function summary where the parameter infos are actually stored. */
extern GTY(()) ipa_node_params_t * ipa_node_params_sum;
/* Call summary to store information about edges such as jump functions.  */
extern GTY(()) ipa_edge_args_sum_t *ipa_edge_args_sum;

/* Function summary for IPA-CP transformation.  */
class ipcp_transformation_t
: public function_summary<ipcp_transformation *>
{
public:
  ipcp_transformation_t (symbol_table *table, bool ggc):
    function_summary<ipcp_transformation *> (table, ggc) {}

  ~ipcp_transformation_t () {}

  static ipcp_transformation_t *create_ggc (symbol_table *symtab)
  {
    ipcp_transformation_t *summary
      = new (ggc_alloc_no_dtor <ipcp_transformation_t> ())
      ipcp_transformation_t (symtab, true);
    return summary;
  }
  /* Hook that is called by summary when a node is duplicated.  */
  virtual void duplicate (cgraph_node *node,
			  cgraph_node *node2,
			  ipcp_transformation *data,
			  ipcp_transformation *data2);
};

/* Function summary where the IPA CP transformations are actually stored.  */
extern GTY(()) function_summary <ipcp_transformation *> *ipcp_transformation_sum;

/* Creating and freeing ipa_node_params and ipa_edge_args.  */
void ipa_create_all_node_params (void);
void ipa_create_all_edge_args (void);
void ipa_check_create_edge_args (void);
void ipa_free_all_node_params (void);
void ipa_free_all_edge_args (void);
void ipa_free_all_structures_after_ipa_cp (void);
void ipa_free_all_structures_after_iinln (void);

void ipa_register_cgraph_hooks (void);
int count_formal_params (tree fndecl);

/* This function ensures the array of node param infos is big enough to
   accommodate a structure for all nodes and reallocates it if not.  */

static inline void
ipa_check_create_node_params (void)
{
  if (!ipa_node_params_sum)
    ipa_node_params_sum
      = (new (ggc_alloc_no_dtor <ipa_node_params_t> ())
	 ipa_node_params_t (symtab, true));
}

/* Returns true if edge summary contains a record for EDGE.  The main purpose
   of this function is that debug dumping function can check info availability
   without causing allocations.  */

static inline bool
ipa_edge_args_info_available_for_edge_p (struct cgraph_edge *edge)
{
  return ipa_edge_args_sum->exists (edge);
}

static inline ipcp_transformation *
ipcp_get_transformation_summary (cgraph_node *node)
{
  if (ipcp_transformation_sum == NULL)
    return NULL;

  return ipcp_transformation_sum->get (node);
}

/* Return the aggregate replacements for NODE, if there are any.  */

static inline struct ipa_agg_replacement_value *
ipa_get_agg_replacements_for_node (cgraph_node *node)
{
  ipcp_transformation *ts = ipcp_get_transformation_summary (node);
  return ts ? ts->agg_values : NULL;
}

/* Function formal parameters related computations.  */
void ipa_initialize_node_params (struct cgraph_node *node);
bool ipa_propagate_indirect_call_infos (struct cgraph_edge *cs,
					vec<cgraph_edge *> *new_edges);

/* Indirect edge processing and target discovery.  */
tree ipa_get_indirect_edge_target (struct cgraph_edge *ie,
				   ipa_call_arg_values *avals,
				   bool *speculative);
tree ipa_get_indirect_edge_target (struct cgraph_edge *ie,
				   ipa_auto_call_arg_values *avals,
				   bool *speculative);
struct cgraph_edge *ipa_make_edge_direct_to_target (struct cgraph_edge *, tree,
						    bool speculative = false);
tree ipa_impossible_devirt_target (struct cgraph_edge *, tree);
ipa_bits *ipa_get_ipa_bits_for_value (const widest_int &value,
				      const widest_int &mask);


/* Functions related to both.  */
void ipa_analyze_node (struct cgraph_node *);

/* Aggregate jump function related functions.  */
tree ipa_find_agg_cst_for_param (const ipa_agg_value_set *agg, tree scalar,
				 HOST_WIDE_INT offset, bool by_ref,
				 bool *from_global_constant = NULL);
bool ipa_load_from_parm_agg (struct ipa_func_body_info *fbi,
			     vec<ipa_param_descriptor, va_gc> *descriptors,
			     gimple *stmt, tree op, int *index_p,
			     HOST_WIDE_INT *offset_p, poly_int64 *size_p,
			     bool *by_ref, bool *guaranteed_unmodified = NULL);

/* Debugging interface.  */
void ipa_print_node_params (FILE *, struct cgraph_node *node);
void ipa_print_all_params (FILE *);
void ipa_print_node_jump_functions (FILE *f, struct cgraph_node *node);
void ipa_print_all_jump_functions (FILE * f);
void ipcp_verify_propagated_values (void);

template <typename value>
class ipcp_value;

extern object_allocator<ipcp_value<tree> > ipcp_cst_values_pool;
extern object_allocator<ipcp_value<ipa_polymorphic_call_context> >
  ipcp_poly_ctx_values_pool;

template <typename valtype>
struct ipcp_value_source;

extern object_allocator<ipcp_value_source<tree> > ipcp_sources_pool;

struct ipcp_agg_lattice;

extern object_allocator<ipcp_agg_lattice> ipcp_agg_lattice_pool;

void ipa_dump_agg_replacement_values (FILE *f,
				      struct ipa_agg_replacement_value *av);
void ipa_prop_write_jump_functions (void);
void ipa_prop_read_jump_functions (void);
void ipcp_write_transformation_summaries (void);
void ipcp_read_transformation_summaries (void);
int ipa_get_param_decl_index (class ipa_node_params *, tree);
tree ipa_value_from_jfunc (class ipa_node_params *info,
			   struct ipa_jump_func *jfunc, tree type);
unsigned int ipcp_transform_function (struct cgraph_node *node);
ipa_polymorphic_call_context ipa_context_from_jfunc (ipa_node_params *,
						     cgraph_edge *,
						     int,
						     ipa_jump_func *);
value_range ipa_value_range_from_jfunc (ipa_node_params *, cgraph_edge *,
					ipa_jump_func *, tree);
ipa_agg_value_set ipa_agg_value_set_from_jfunc (ipa_node_params *,
						cgraph_node *,
						ipa_agg_jump_function *);
void ipa_dump_param (FILE *, class ipa_node_params *info, int i);
void ipa_release_body_info (struct ipa_func_body_info *);
tree ipa_get_callee_param_type (struct cgraph_edge *e, int i);
bool ipcp_get_parm_bits (tree, tree *, widest_int *);
bool unadjusted_ptr_and_unit_offset (tree op, tree *ret,
				     poly_int64 *offset_ret);

/* From tree-sra.cc:  */
tree build_ref_for_offset (location_t, tree, poly_int64, bool, tree,
			   gimple_stmt_iterator *, bool);

/* In ipa-cp.cc  */
void ipa_cp_cc_finalize (void);

#endif /* IPA_PROP_H */
