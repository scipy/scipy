/* Vectorizer
   Copyright (C) 2003-2022 Free Software Foundation, Inc.
   Contributed by Dorit Naishlos <dorit@il.ibm.com>

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

#ifndef GCC_TREE_VECTORIZER_H
#define GCC_TREE_VECTORIZER_H

typedef class _stmt_vec_info *stmt_vec_info;
typedef struct _slp_tree *slp_tree;

#include "tree-data-ref.h"
#include "tree-hash-traits.h"
#include "target.h"
#include "internal-fn.h"
#include "tree-ssa-operands.h"
#include "gimple-match.h"

/* Used for naming of new temporaries.  */
enum vect_var_kind {
  vect_simple_var,
  vect_pointer_var,
  vect_scalar_var,
  vect_mask_var
};

/* Defines type of operation.  */
enum operation_type {
  unary_op = 1,
  binary_op,
  ternary_op
};

/* Define type of available alignment support.  */
enum dr_alignment_support {
  dr_unaligned_unsupported,
  dr_unaligned_supported,
  dr_explicit_realign,
  dr_explicit_realign_optimized,
  dr_aligned
};

/* Define type of def-use cross-iteration cycle.  */
enum vect_def_type {
  vect_uninitialized_def = 0,
  vect_constant_def = 1,
  vect_external_def,
  vect_internal_def,
  vect_induction_def,
  vect_reduction_def,
  vect_double_reduction_def,
  vect_nested_cycle,
  vect_unknown_def_type
};

/* Define type of reduction.  */
enum vect_reduction_type {
  TREE_CODE_REDUCTION,
  COND_REDUCTION,
  INTEGER_INDUC_COND_REDUCTION,
  CONST_COND_REDUCTION,

  /* Retain a scalar phi and use a FOLD_EXTRACT_LAST within the loop
     to implement:

       for (int i = 0; i < VF; ++i)
         res = cond[i] ? val[i] : res;  */
  EXTRACT_LAST_REDUCTION,

  /* Use a folding reduction within the loop to implement:

       for (int i = 0; i < VF; ++i)
	 res = res OP val[i];

     (with no reassocation).  */
  FOLD_LEFT_REDUCTION
};

#define VECTORIZABLE_CYCLE_DEF(D) (((D) == vect_reduction_def)           \
                                   || ((D) == vect_double_reduction_def) \
                                   || ((D) == vect_nested_cycle))

/* Structure to encapsulate information about a group of like
   instructions to be presented to the target cost model.  */
struct stmt_info_for_cost {
  int count;
  enum vect_cost_for_stmt kind;
  enum vect_cost_model_location where;
  stmt_vec_info stmt_info;
  slp_tree node;
  tree vectype;
  int misalign;
};

typedef vec<stmt_info_for_cost> stmt_vector_for_cost;

/* Maps base addresses to an innermost_loop_behavior and the stmt it was
   derived from that gives the maximum known alignment for that base.  */
typedef hash_map<tree_operand_hash,
		 std::pair<stmt_vec_info, innermost_loop_behavior *> >
	  vec_base_alignments;

/* Represents elements [START, START + LENGTH) of cyclical array OPS*
   (i.e. OPS repeated to give at least START + LENGTH elements)  */
struct vect_scalar_ops_slice
{
  tree op (unsigned int i) const;
  bool all_same_p () const;

  vec<tree> *ops;
  unsigned int start;
  unsigned int length;
};

/* Return element I of the slice.  */
inline tree
vect_scalar_ops_slice::op (unsigned int i) const
{
  return (*ops)[(i + start) % ops->length ()];
}

/* Hash traits for vect_scalar_ops_slice.  */
struct vect_scalar_ops_slice_hash : typed_noop_remove<vect_scalar_ops_slice>
{
  typedef vect_scalar_ops_slice value_type;
  typedef vect_scalar_ops_slice compare_type;

  static const bool empty_zero_p = true;

  static void mark_deleted (value_type &s) { s.length = ~0U; }
  static void mark_empty (value_type &s) { s.length = 0; }
  static bool is_deleted (const value_type &s) { return s.length == ~0U; }
  static bool is_empty (const value_type &s) { return s.length == 0; }
  static hashval_t hash (const value_type &);
  static bool equal (const value_type &, const compare_type &);
};

/************************************************************************
  SLP
 ************************************************************************/
typedef vec<std::pair<unsigned, unsigned> > lane_permutation_t;
typedef vec<unsigned> load_permutation_t;

/* A computation tree of an SLP instance.  Each node corresponds to a group of
   stmts to be packed in a SIMD stmt.  */
struct _slp_tree {
  _slp_tree ();
  ~_slp_tree ();

  /* Nodes that contain def-stmts of this node statements operands.  */
  vec<slp_tree> children;

  /* A group of scalar stmts to be vectorized together.  */
  vec<stmt_vec_info> stmts;
  /* A group of scalar operands to be vectorized together.  */
  vec<tree> ops;
  /* The representative that should be used for analysis and
     code generation.  */
  stmt_vec_info representative;

  /* Load permutation relative to the stores, NULL if there is no
     permutation.  */
  load_permutation_t load_permutation;
  /* Lane permutation of the operands scalar lanes encoded as pairs
     of { operand number, lane number }.  The number of elements
     denotes the number of output lanes.  */
  lane_permutation_t lane_permutation;

  tree vectype;
  /* Vectorized stmt/s.  */
  vec<gimple *> vec_stmts;
  vec<tree> vec_defs;
  /* Number of vector stmts that are created to replace the group of scalar
     stmts. It is calculated during the transformation phase as the number of
     scalar elements in one scalar iteration (GROUP_SIZE) multiplied by VF
     divided by vector size.  */
  unsigned int vec_stmts_size;

  /* Reference count in the SLP graph.  */
  unsigned int refcnt;
  /* The maximum number of vector elements for the subtree rooted
     at this node.  */
  poly_uint64 max_nunits;
  /* The DEF type of this node.  */
  enum vect_def_type def_type;
  /* The number of scalar lanes produced by this node.  */
  unsigned int lanes;
  /* The operation of this node.  */
  enum tree_code code;

  int vertex;

  /* If not NULL this is a cached failed SLP discovery attempt with
     the lanes that failed during SLP discovery as 'false'.  This is
     a copy of the matches array.  */
  bool *failed;

  /* Allocate from slp_tree_pool.  */
  static void *operator new (size_t);

  /* Return memory to slp_tree_pool.  */
  static void operator delete (void *, size_t);

  /* Linked list of nodes to release when we free the slp_tree_pool.  */
  slp_tree next_node;
  slp_tree prev_node;
};

/* The enum describes the type of operations that an SLP instance
   can perform. */

enum slp_instance_kind {
    slp_inst_kind_store,
    slp_inst_kind_reduc_group,
    slp_inst_kind_reduc_chain,
    slp_inst_kind_bb_reduc,
    slp_inst_kind_ctor
};

/* SLP instance is a sequence of stmts in a loop that can be packed into
   SIMD stmts.  */
typedef class _slp_instance {
public:
  /* The root of SLP tree.  */
  slp_tree root;

  /* For vector constructors, the constructor stmt that the SLP tree is built
     from, NULL otherwise.  */
  vec<stmt_vec_info> root_stmts;

  /* The unrolling factor required to vectorized this SLP instance.  */
  poly_uint64 unrolling_factor;

  /* The group of nodes that contain loads of this SLP instance.  */
  vec<slp_tree> loads;

  /* The SLP node containing the reduction PHIs.  */
  slp_tree reduc_phis;

  /* Vector cost of this entry to the SLP graph.  */
  stmt_vector_for_cost cost_vec;

  /* If this instance is the main entry of a subgraph the set of
     entries into the same subgraph, including itself.  */
  vec<_slp_instance *> subgraph_entries;

  /* The type of operation the SLP instance is performing.  */
  slp_instance_kind kind;

  dump_user_location_t location () const;
} *slp_instance;


/* Access Functions.  */
#define SLP_INSTANCE_TREE(S)                     (S)->root
#define SLP_INSTANCE_UNROLLING_FACTOR(S)         (S)->unrolling_factor
#define SLP_INSTANCE_LOADS(S)                    (S)->loads
#define SLP_INSTANCE_ROOT_STMTS(S)               (S)->root_stmts
#define SLP_INSTANCE_KIND(S)                     (S)->kind

#define SLP_TREE_CHILDREN(S)                     (S)->children
#define SLP_TREE_SCALAR_STMTS(S)                 (S)->stmts
#define SLP_TREE_SCALAR_OPS(S)                   (S)->ops
#define SLP_TREE_REF_COUNT(S)                    (S)->refcnt
#define SLP_TREE_VEC_STMTS(S)                    (S)->vec_stmts
#define SLP_TREE_VEC_DEFS(S)                     (S)->vec_defs
#define SLP_TREE_NUMBER_OF_VEC_STMTS(S)          (S)->vec_stmts_size
#define SLP_TREE_LOAD_PERMUTATION(S)             (S)->load_permutation
#define SLP_TREE_LANE_PERMUTATION(S)             (S)->lane_permutation
#define SLP_TREE_DEF_TYPE(S)			 (S)->def_type
#define SLP_TREE_VECTYPE(S)			 (S)->vectype
#define SLP_TREE_REPRESENTATIVE(S)		 (S)->representative
#define SLP_TREE_LANES(S)			 (S)->lanes
#define SLP_TREE_CODE(S)			 (S)->code

/* Key for map that records association between
   scalar conditions and corresponding loop mask, and
   is populated by vect_record_loop_mask.  */

struct scalar_cond_masked_key
{
  scalar_cond_masked_key (tree t, unsigned ncopies_)
    : ncopies (ncopies_)
  {
    get_cond_ops_from_tree (t);
  }

  void get_cond_ops_from_tree (tree);

  unsigned ncopies;
  bool inverted_p;
  tree_code code;
  tree op0;
  tree op1;
};

template<>
struct default_hash_traits<scalar_cond_masked_key>
{
  typedef scalar_cond_masked_key compare_type;
  typedef scalar_cond_masked_key value_type;

  static inline hashval_t
  hash (value_type v)
  {
    inchash::hash h;
    h.add_int (v.code);
    inchash::add_expr (v.op0, h, 0);
    inchash::add_expr (v.op1, h, 0);
    h.add_int (v.ncopies);
    h.add_flag (v.inverted_p);
    return h.end ();
  }

  static inline bool
  equal (value_type existing, value_type candidate)
  {
    return (existing.ncopies == candidate.ncopies
	    && existing.code == candidate.code
	    && existing.inverted_p == candidate.inverted_p
	    && operand_equal_p (existing.op0, candidate.op0, 0)
	    && operand_equal_p (existing.op1, candidate.op1, 0));
  }

  static const bool empty_zero_p = true;

  static inline void
  mark_empty (value_type &v)
  {
    v.ncopies = 0;
    v.inverted_p = false;
  }

  static inline bool
  is_empty (value_type v)
  {
    return v.ncopies == 0;
  }

  static inline void mark_deleted (value_type &) {}

  static inline bool is_deleted (const value_type &)
  {
    return false;
  }

  static inline void remove (value_type &) {}
};

typedef hash_set<scalar_cond_masked_key> scalar_cond_masked_set_type;

/* Key and map that records association between vector conditions and
   corresponding loop mask, and is populated by prepare_vec_mask.  */

typedef pair_hash<tree_operand_hash, tree_operand_hash> tree_cond_mask_hash;
typedef hash_set<tree_cond_mask_hash> vec_cond_masked_set_type;

/* Describes two objects whose addresses must be unequal for the vectorized
   loop to be valid.  */
typedef std::pair<tree, tree> vec_object_pair;

/* Records that vectorization is only possible if abs (EXPR) >= MIN_VALUE.
   UNSIGNED_P is true if we can assume that abs (EXPR) == EXPR.  */
class vec_lower_bound {
public:
  vec_lower_bound () {}
  vec_lower_bound (tree e, bool u, poly_uint64 m)
    : expr (e), unsigned_p (u), min_value (m) {}

  tree expr;
  bool unsigned_p;
  poly_uint64 min_value;
};

/* Vectorizer state shared between different analyses like vector sizes
   of the same CFG region.  */
class vec_info_shared {
public:
  vec_info_shared();
  ~vec_info_shared();

  void save_datarefs();
  void check_datarefs();

  /* The number of scalar stmts.  */
  unsigned n_stmts;

  /* All data references.  Freed by free_data_refs, so not an auto_vec.  */
  vec<data_reference_p> datarefs;
  vec<data_reference> datarefs_copy;

  /* The loop nest in which the data dependences are computed.  */
  auto_vec<loop_p> loop_nest;

  /* All data dependences.  Freed by free_dependence_relations, so not
     an auto_vec.  */
  vec<ddr_p> ddrs;
};

/* Vectorizer state common between loop and basic-block vectorization.  */
class vec_info {
public:
  typedef hash_set<int_hash<machine_mode, E_VOIDmode, E_BLKmode> > mode_set;
  enum vec_kind { bb, loop };

  vec_info (vec_kind, vec_info_shared *);
  ~vec_info ();

  stmt_vec_info add_stmt (gimple *);
  stmt_vec_info add_pattern_stmt (gimple *, stmt_vec_info);
  stmt_vec_info lookup_stmt (gimple *);
  stmt_vec_info lookup_def (tree);
  stmt_vec_info lookup_single_use (tree);
  class dr_vec_info *lookup_dr (data_reference *);
  void move_dr (stmt_vec_info, stmt_vec_info);
  void remove_stmt (stmt_vec_info);
  void replace_stmt (gimple_stmt_iterator *, stmt_vec_info, gimple *);
  void insert_on_entry (stmt_vec_info, gimple *);
  void insert_seq_on_entry (stmt_vec_info, gimple_seq);

  /* The type of vectorization.  */
  vec_kind kind;

  /* Shared vectorizer state.  */
  vec_info_shared *shared;

  /* The mapping of GIMPLE UID to stmt_vec_info.  */
  vec<stmt_vec_info> stmt_vec_infos;
  /* Whether the above mapping is complete.  */
  bool stmt_vec_info_ro;

  /* The SLP graph.  */
  auto_vec<slp_instance> slp_instances;

  /* Maps base addresses to an innermost_loop_behavior that gives the maximum
     known alignment for that base.  */
  vec_base_alignments base_alignments;

  /* All interleaving chains of stores, represented by the first
     stmt in the chain.  */
  auto_vec<stmt_vec_info> grouped_stores;

  /* The set of vector modes used in the vectorized region.  */
  mode_set used_vector_modes;

  /* The argument we should pass to related_vector_mode when looking up
     the vector mode for a scalar mode, or VOIDmode if we haven't yet
     made any decisions about which vector modes to use.  */
  machine_mode vector_mode;

private:
  stmt_vec_info new_stmt_vec_info (gimple *stmt);
  void set_vinfo_for_stmt (gimple *, stmt_vec_info, bool = true);
  void free_stmt_vec_infos ();
  void free_stmt_vec_info (stmt_vec_info);
};

class _loop_vec_info;
class _bb_vec_info;

template<>
template<>
inline bool
is_a_helper <_loop_vec_info *>::test (vec_info *i)
{
  return i->kind == vec_info::loop;
}

template<>
template<>
inline bool
is_a_helper <_bb_vec_info *>::test (vec_info *i)
{
  return i->kind == vec_info::bb;
}

/* In general, we can divide the vector statements in a vectorized loop
   into related groups ("rgroups") and say that for each rgroup there is
   some nS such that the rgroup operates on nS values from one scalar
   iteration followed by nS values from the next.  That is, if VF is the
   vectorization factor of the loop, the rgroup operates on a sequence:

     (1,1) (1,2) ... (1,nS) (2,1) ... (2,nS) ... (VF,1) ... (VF,nS)

   where (i,j) represents a scalar value with index j in a scalar
   iteration with index i.

   [ We use the term "rgroup" to emphasise that this grouping isn't
     necessarily the same as the grouping of statements used elsewhere.
     For example, if we implement a group of scalar loads using gather
     loads, we'll use a separate gather load for each scalar load, and
     thus each gather load will belong to its own rgroup. ]

   In general this sequence will occupy nV vectors concatenated
   together.  If these vectors have nL lanes each, the total number
   of scalar values N is given by:

       N = nS * VF = nV * nL

   None of nS, VF, nV and nL are required to be a power of 2.  nS and nV
   are compile-time constants but VF and nL can be variable (if the target
   supports variable-length vectors).

   In classical vectorization, each iteration of the vector loop would
   handle exactly VF iterations of the original scalar loop.  However,
   in vector loops that are able to operate on partial vectors, a
   particular iteration of the vector loop might handle fewer than VF
   iterations of the scalar loop.  The vector lanes that correspond to
   iterations of the scalar loop are said to be "active" and the other
   lanes are said to be "inactive".

   In such vector loops, many rgroups need to be controlled to ensure
   that they have no effect for the inactive lanes.  Conceptually, each
   such rgroup needs a sequence of booleans in the same order as above,
   but with each (i,j) replaced by a boolean that indicates whether
   iteration i is active.  This sequence occupies nV vector controls
   that again have nL lanes each.  Thus the control sequence as a whole
   consists of VF independent booleans that are each repeated nS times.

   Taking mask-based approach as a partially-populated vectors example.
   We make the simplifying assumption that if a sequence of nV masks is
   suitable for one (nS,nL) pair, we can reuse it for (nS/2,nL/2) by
   VIEW_CONVERTing it.  This holds for all current targets that support
   fully-masked loops.  For example, suppose the scalar loop is:

     float *f;
     double *d;
     for (int i = 0; i < n; ++i)
       {
	 f[i * 2 + 0] += 1.0f;
	 f[i * 2 + 1] += 2.0f;
	 d[i] += 3.0;
       }

   and suppose that vectors have 256 bits.  The vectorized f accesses
   will belong to one rgroup and the vectorized d access to another:

     f rgroup: nS = 2, nV = 1, nL = 8
     d rgroup: nS = 1, nV = 1, nL = 4
	       VF = 4

     [ In this simple example the rgroups do correspond to the normal
       SLP grouping scheme. ]

   If only the first three lanes are active, the masks we need are:

     f rgroup: 1 1 | 1 1 | 1 1 | 0 0
     d rgroup:  1  |  1  |  1  |  0

   Here we can use a mask calculated for f's rgroup for d's, but not
   vice versa.

   Thus for each value of nV, it is enough to provide nV masks, with the
   mask being calculated based on the highest nL (or, equivalently, based
   on the highest nS) required by any rgroup with that nV.  We therefore
   represent the entire collection of masks as a two-level table, with the
   first level being indexed by nV - 1 (since nV == 0 doesn't exist) and
   the second being indexed by the mask index 0 <= i < nV.  */

/* The controls (like masks or lengths) needed by rgroups with nV vectors,
   according to the description above.  */
struct rgroup_controls {
  /* The largest nS for all rgroups that use these controls.  */
  unsigned int max_nscalars_per_iter;

  /* For the largest nS recorded above, the loop controls divide each scalar
     into FACTOR equal-sized pieces.  This is useful if we need to split
     element-based accesses into byte-based accesses.  */
  unsigned int factor;

  /* This is a vector type with MAX_NSCALARS_PER_ITER * VF / nV elements.
     For mask-based controls, it is the type of the masks in CONTROLS.
     For length-based controls, it can be any vector type that has the
     specified number of elements; the type of the elements doesn't matter.  */
  tree type;

  /* A vector of nV controls, in iteration order.  */
  vec<tree> controls;

  /* In case of len_load and len_store with a bias there is only one
     rgroup.  This holds the adjusted loop length for the this rgroup.  */
  tree bias_adjusted_ctrl;
};

typedef auto_vec<rgroup_controls> vec_loop_masks;

typedef auto_vec<rgroup_controls> vec_loop_lens;

typedef auto_vec<std::pair<data_reference*, tree> > drs_init_vec;

/* Information about a reduction accumulator from the main loop that could
   conceivably be reused as the input to a reduction in an epilogue loop.  */
struct vect_reusable_accumulator {
  /* The final value of the accumulator, which forms the input to the
     reduction operation.  */
  tree reduc_input;

  /* The stmt_vec_info that describes the reduction (i.e. the one for
     which is_reduc_info is true).  */
  stmt_vec_info reduc_info;
};

/*-----------------------------------------------------------------*/
/* Info on vectorized loops.                                       */
/*-----------------------------------------------------------------*/
typedef class _loop_vec_info : public vec_info {
public:
  _loop_vec_info (class loop *, vec_info_shared *);
  ~_loop_vec_info ();

  /* The loop to which this info struct refers to.  */
  class loop *loop;

  /* The loop basic blocks.  */
  basic_block *bbs;

  /* Number of latch executions.  */
  tree num_itersm1;
  /* Number of iterations.  */
  tree num_iters;
  /* Number of iterations of the original loop.  */
  tree num_iters_unchanged;
  /* Condition under which this loop is analyzed and versioned.  */
  tree num_iters_assumptions;

  /* The cost of the vector code.  */
  class vector_costs *vector_costs;

  /* The cost of the scalar code.  */
  class vector_costs *scalar_costs;

  /* Threshold of number of iterations below which vectorization will not be
     performed. It is calculated from MIN_PROFITABLE_ITERS and
     param_min_vect_loop_bound.  */
  unsigned int th;

  /* When applying loop versioning, the vector form should only be used
     if the number of scalar iterations is >= this value, on top of all
     the other requirements.  Ignored when loop versioning is not being
     used.  */
  poly_uint64 versioning_threshold;

  /* Unrolling factor  */
  poly_uint64 vectorization_factor;

  /* If this loop is an epilogue loop whose main loop can be skipped,
     MAIN_LOOP_EDGE is the edge from the main loop to this loop's
     preheader.  SKIP_MAIN_LOOP_EDGE is then the edge that skips the
     main loop and goes straight to this loop's preheader.

     Both fields are null otherwise.  */
  edge main_loop_edge;
  edge skip_main_loop_edge;

  /* If this loop is an epilogue loop that might be skipped after executing
     the main loop, this edge is the one that skips the epilogue.  */
  edge skip_this_loop_edge;

  /* The vectorized form of a standard reduction replaces the original
     scalar code's final result (a loop-closed SSA PHI) with the result
     of a vector-to-scalar reduction operation.  After vectorization,
     this variable maps these vector-to-scalar results to information
     about the reductions that generated them.  */
  hash_map<tree, vect_reusable_accumulator> reusable_accumulators;

  /* The number of times that the target suggested we unroll the vector loop
     in order to promote more ILP.  This value will be used to re-analyze the
     loop for vectorization and if successful the value will be folded into
     vectorization_factor (and therefore exactly divides
     vectorization_factor).  */
  unsigned int suggested_unroll_factor;

  /* Maximum runtime vectorization factor, or MAX_VECTORIZATION_FACTOR
     if there is no particular limit.  */
  unsigned HOST_WIDE_INT max_vectorization_factor;

  /* The masks that a fully-masked loop should use to avoid operating
     on inactive scalars.  */
  vec_loop_masks masks;

  /* The lengths that a loop with length should use to avoid operating
     on inactive scalars.  */
  vec_loop_lens lens;

  /* Set of scalar conditions that have loop mask applied.  */
  scalar_cond_masked_set_type scalar_cond_masked_set;

  /* Set of vector conditions that have loop mask applied.  */
  vec_cond_masked_set_type vec_cond_masked_set;

  /* If we are using a loop mask to align memory addresses, this variable
     contains the number of vector elements that we should skip in the
     first iteration of the vector loop (i.e. the number of leading
     elements that should be false in the first mask).  */
  tree mask_skip_niters;

  /* The type that the loop control IV should be converted to before
     testing which of the VF scalars are active and inactive.
     Only meaningful if LOOP_VINFO_USING_PARTIAL_VECTORS_P.  */
  tree rgroup_compare_type;

  /* For #pragma omp simd if (x) loops the x expression.  If constant 0,
     the loop should not be vectorized, if constant non-zero, simd_if_cond
     shouldn't be set and loop vectorized normally, if SSA_NAME, the loop
     should be versioned on that condition, using scalar loop if the condition
     is false and vectorized loop otherwise.  */
  tree simd_if_cond;

  /* The type that the vector loop control IV should have when
     LOOP_VINFO_USING_PARTIAL_VECTORS_P is true.  */
  tree rgroup_iv_type;

  /* Unknown DRs according to which loop was peeled.  */
  class dr_vec_info *unaligned_dr;

  /* peeling_for_alignment indicates whether peeling for alignment will take
     place, and what the peeling factor should be:
     peeling_for_alignment = X means:
        If X=0: Peeling for alignment will not be applied.
        If X>0: Peel first X iterations.
        If X=-1: Generate a runtime test to calculate the number of iterations
                 to be peeled, using the dataref recorded in the field
                 unaligned_dr.  */
  int peeling_for_alignment;

  /* The mask used to check the alignment of pointers or arrays.  */
  int ptr_mask;

  /* Data Dependence Relations defining address ranges that are candidates
     for a run-time aliasing check.  */
  auto_vec<ddr_p> may_alias_ddrs;

  /* Data Dependence Relations defining address ranges together with segment
     lengths from which the run-time aliasing check is built.  */
  auto_vec<dr_with_seg_len_pair_t> comp_alias_ddrs;

  /* Check that the addresses of each pair of objects is unequal.  */
  auto_vec<vec_object_pair> check_unequal_addrs;

  /* List of values that are required to be nonzero.  This is used to check
     whether things like "x[i * n] += 1;" are safe and eventually gets added
     to the checks for lower bounds below.  */
  auto_vec<tree> check_nonzero;

  /* List of values that need to be checked for a minimum value.  */
  auto_vec<vec_lower_bound> lower_bounds;

  /* Statements in the loop that have data references that are candidates for a
     runtime (loop versioning) misalignment check.  */
  auto_vec<stmt_vec_info> may_misalign_stmts;

  /* Reduction cycles detected in the loop. Used in loop-aware SLP.  */
  auto_vec<stmt_vec_info> reductions;

  /* All reduction chains in the loop, represented by the first
     stmt in the chain.  */
  auto_vec<stmt_vec_info> reduction_chains;

  /* Cost vector for a single scalar iteration.  */
  auto_vec<stmt_info_for_cost> scalar_cost_vec;

  /* Map of IV base/step expressions to inserted name in the preheader.  */
  hash_map<tree_operand_hash, tree> *ivexpr_map;

  /* Map of OpenMP "omp simd array" scan variables to corresponding
     rhs of the store of the initializer.  */
  hash_map<tree, tree> *scan_map;

  /* The unrolling factor needed to SLP the loop. In case of that pure SLP is
     applied to the loop, i.e., no unrolling is needed, this is 1.  */
  poly_uint64 slp_unrolling_factor;

  /* The factor used to over weight those statements in an inner loop
     relative to the loop being vectorized.  */
  unsigned int inner_loop_cost_factor;

  /* Is the loop vectorizable? */
  bool vectorizable;

  /* Records whether we still have the option of vectorizing this loop
     using partially-populated vectors; in other words, whether it is
     still possible for one iteration of the vector loop to handle
     fewer than VF scalars.  */
  bool can_use_partial_vectors_p;

  /* True if we've decided to use partially-populated vectors, so that
     the vector loop can handle fewer than VF scalars.  */
  bool using_partial_vectors_p;

  /* True if we've decided to use partially-populated vectors for the
     epilogue of loop.  */
  bool epil_using_partial_vectors_p;

  /* The bias for len_load and len_store.  For now, only 0 and -1 are
     supported.  -1 must be used when a backend does not support
     len_load/len_store with a length of zero.  */
  signed char partial_load_store_bias;

  /* When we have grouped data accesses with gaps, we may introduce invalid
     memory accesses.  We peel the last iteration of the loop to prevent
     this.  */
  bool peeling_for_gaps;

  /* When the number of iterations is not a multiple of the vector size
     we need to peel off iterations at the end to form an epilogue loop.  */
  bool peeling_for_niter;

  /* True if there are no loop carried data dependencies in the loop.
     If loop->safelen <= 1, then this is always true, either the loop
     didn't have any loop carried data dependencies, or the loop is being
     vectorized guarded with some runtime alias checks, or couldn't
     be vectorized at all, but then this field shouldn't be used.
     For loop->safelen >= 2, the user has asserted that there are no
     backward dependencies, but there still could be loop carried forward
     dependencies in such loops.  This flag will be false if normal
     vectorizer data dependency analysis would fail or require versioning
     for alias, but because of loop->safelen >= 2 it has been vectorized
     even without versioning for alias.  E.g. in:
     #pragma omp simd
     for (int i = 0; i < m; i++)
       a[i] = a[i + k] * c;
     (or #pragma simd or #pragma ivdep) we can vectorize this and it will
     DTRT even for k > 0 && k < m, but without safelen we would not
     vectorize this, so this field would be false.  */
  bool no_data_dependencies;

  /* Mark loops having masked stores.  */
  bool has_mask_store;

  /* Queued scaling factor for the scalar loop.  */
  profile_probability scalar_loop_scaling;

  /* If if-conversion versioned this loop before conversion, this is the
     loop version without if-conversion.  */
  class loop *scalar_loop;

  /* For loops being epilogues of already vectorized loops
     this points to the original vectorized loop.  Otherwise NULL.  */
  _loop_vec_info *orig_loop_info;

  /* Used to store loop_vec_infos of epilogues of this loop during
     analysis.  */
  vec<_loop_vec_info *> epilogue_vinfos;

} *loop_vec_info;

/* Access Functions.  */
#define LOOP_VINFO_LOOP(L)                 (L)->loop
#define LOOP_VINFO_BBS(L)                  (L)->bbs
#define LOOP_VINFO_NITERSM1(L)             (L)->num_itersm1
#define LOOP_VINFO_NITERS(L)               (L)->num_iters
/* Since LOOP_VINFO_NITERS and LOOP_VINFO_NITERSM1 can change after
   prologue peeling retain total unchanged scalar loop iterations for
   cost model.  */
#define LOOP_VINFO_NITERS_UNCHANGED(L)     (L)->num_iters_unchanged
#define LOOP_VINFO_NITERS_ASSUMPTIONS(L)   (L)->num_iters_assumptions
#define LOOP_VINFO_COST_MODEL_THRESHOLD(L) (L)->th
#define LOOP_VINFO_VERSIONING_THRESHOLD(L) (L)->versioning_threshold
#define LOOP_VINFO_VECTORIZABLE_P(L)       (L)->vectorizable
#define LOOP_VINFO_CAN_USE_PARTIAL_VECTORS_P(L) (L)->can_use_partial_vectors_p
#define LOOP_VINFO_USING_PARTIAL_VECTORS_P(L) (L)->using_partial_vectors_p
#define LOOP_VINFO_EPIL_USING_PARTIAL_VECTORS_P(L)                             \
  (L)->epil_using_partial_vectors_p
#define LOOP_VINFO_PARTIAL_LOAD_STORE_BIAS(L) (L)->partial_load_store_bias
#define LOOP_VINFO_VECT_FACTOR(L)          (L)->vectorization_factor
#define LOOP_VINFO_MAX_VECT_FACTOR(L)      (L)->max_vectorization_factor
#define LOOP_VINFO_MASKS(L)                (L)->masks
#define LOOP_VINFO_LENS(L)                 (L)->lens
#define LOOP_VINFO_MASK_SKIP_NITERS(L)     (L)->mask_skip_niters
#define LOOP_VINFO_RGROUP_COMPARE_TYPE(L)  (L)->rgroup_compare_type
#define LOOP_VINFO_RGROUP_IV_TYPE(L)       (L)->rgroup_iv_type
#define LOOP_VINFO_PTR_MASK(L)             (L)->ptr_mask
#define LOOP_VINFO_N_STMTS(L)		   (L)->shared->n_stmts
#define LOOP_VINFO_LOOP_NEST(L)            (L)->shared->loop_nest
#define LOOP_VINFO_DATAREFS(L)             (L)->shared->datarefs
#define LOOP_VINFO_DDRS(L)                 (L)->shared->ddrs
#define LOOP_VINFO_INT_NITERS(L)           (TREE_INT_CST_LOW ((L)->num_iters))
#define LOOP_VINFO_PEELING_FOR_ALIGNMENT(L) (L)->peeling_for_alignment
#define LOOP_VINFO_UNALIGNED_DR(L)         (L)->unaligned_dr
#define LOOP_VINFO_MAY_MISALIGN_STMTS(L)   (L)->may_misalign_stmts
#define LOOP_VINFO_MAY_ALIAS_DDRS(L)       (L)->may_alias_ddrs
#define LOOP_VINFO_COMP_ALIAS_DDRS(L)      (L)->comp_alias_ddrs
#define LOOP_VINFO_CHECK_UNEQUAL_ADDRS(L)  (L)->check_unequal_addrs
#define LOOP_VINFO_CHECK_NONZERO(L)        (L)->check_nonzero
#define LOOP_VINFO_LOWER_BOUNDS(L)         (L)->lower_bounds
#define LOOP_VINFO_GROUPED_STORES(L)       (L)->grouped_stores
#define LOOP_VINFO_SLP_INSTANCES(L)        (L)->slp_instances
#define LOOP_VINFO_SLP_UNROLLING_FACTOR(L) (L)->slp_unrolling_factor
#define LOOP_VINFO_REDUCTIONS(L)           (L)->reductions
#define LOOP_VINFO_REDUCTION_CHAINS(L)     (L)->reduction_chains
#define LOOP_VINFO_PEELING_FOR_GAPS(L)     (L)->peeling_for_gaps
#define LOOP_VINFO_PEELING_FOR_NITER(L)    (L)->peeling_for_niter
#define LOOP_VINFO_NO_DATA_DEPENDENCIES(L) (L)->no_data_dependencies
#define LOOP_VINFO_SCALAR_LOOP(L)	   (L)->scalar_loop
#define LOOP_VINFO_SCALAR_LOOP_SCALING(L)  (L)->scalar_loop_scaling
#define LOOP_VINFO_HAS_MASK_STORE(L)       (L)->has_mask_store
#define LOOP_VINFO_SCALAR_ITERATION_COST(L) (L)->scalar_cost_vec
#define LOOP_VINFO_ORIG_LOOP_INFO(L)       (L)->orig_loop_info
#define LOOP_VINFO_SIMD_IF_COND(L)         (L)->simd_if_cond
#define LOOP_VINFO_INNER_LOOP_COST_FACTOR(L) (L)->inner_loop_cost_factor

#define LOOP_VINFO_FULLY_MASKED_P(L)		\
  (LOOP_VINFO_USING_PARTIAL_VECTORS_P (L)	\
   && !LOOP_VINFO_MASKS (L).is_empty ())

#define LOOP_VINFO_FULLY_WITH_LENGTH_P(L)	\
  (LOOP_VINFO_USING_PARTIAL_VECTORS_P (L)	\
   && !LOOP_VINFO_LENS (L).is_empty ())

#define LOOP_REQUIRES_VERSIONING_FOR_ALIGNMENT(L)	\
  ((L)->may_misalign_stmts.length () > 0)
#define LOOP_REQUIRES_VERSIONING_FOR_ALIAS(L)		\
  ((L)->comp_alias_ddrs.length () > 0 \
   || (L)->check_unequal_addrs.length () > 0 \
   || (L)->lower_bounds.length () > 0)
#define LOOP_REQUIRES_VERSIONING_FOR_NITERS(L)		\
  (LOOP_VINFO_NITERS_ASSUMPTIONS (L))
#define LOOP_REQUIRES_VERSIONING_FOR_SIMD_IF_COND(L)	\
  (LOOP_VINFO_SIMD_IF_COND (L))
#define LOOP_REQUIRES_VERSIONING(L)			\
  (LOOP_REQUIRES_VERSIONING_FOR_ALIGNMENT (L)		\
   || LOOP_REQUIRES_VERSIONING_FOR_ALIAS (L)		\
   || LOOP_REQUIRES_VERSIONING_FOR_NITERS (L)		\
   || LOOP_REQUIRES_VERSIONING_FOR_SIMD_IF_COND (L))

#define LOOP_VINFO_NITERS_KNOWN_P(L)          \
  (tree_fits_shwi_p ((L)->num_iters) && tree_to_shwi ((L)->num_iters) > 0)

#define LOOP_VINFO_EPILOGUE_P(L) \
  (LOOP_VINFO_ORIG_LOOP_INFO (L) != NULL)

#define LOOP_VINFO_ORIG_MAX_VECT_FACTOR(L) \
  (LOOP_VINFO_MAX_VECT_FACTOR (LOOP_VINFO_ORIG_LOOP_INFO (L)))

/* Wrapper for loop_vec_info, for tracking success/failure, where a non-NULL
   value signifies success, and a NULL value signifies failure, supporting
   propagating an opt_problem * describing the failure back up the call
   stack.  */
typedef opt_pointer_wrapper <loop_vec_info> opt_loop_vec_info;

static inline loop_vec_info
loop_vec_info_for_loop (class loop *loop)
{
  return (loop_vec_info) loop->aux;
}

struct slp_root
{
  slp_root (slp_instance_kind kind_, vec<stmt_vec_info> stmts_,
	    vec<stmt_vec_info> roots_)
    : kind(kind_), stmts(stmts_), roots(roots_) {}
  slp_instance_kind kind;
  vec<stmt_vec_info> stmts;
  vec<stmt_vec_info> roots;
};

typedef class _bb_vec_info : public vec_info
{
public:
  _bb_vec_info (vec<basic_block> bbs, vec_info_shared *);
  ~_bb_vec_info ();

  /* The region we are operating on.  bbs[0] is the entry, excluding
     its PHI nodes.  In the future we might want to track an explicit
     entry edge to cover bbs[0] PHI nodes and have a region entry
     insert location.  */
  vec<basic_block> bbs;

  vec<slp_root> roots;
} *bb_vec_info;

#define BB_VINFO_BB(B)               (B)->bb
#define BB_VINFO_GROUPED_STORES(B)   (B)->grouped_stores
#define BB_VINFO_SLP_INSTANCES(B)    (B)->slp_instances
#define BB_VINFO_DATAREFS(B)         (B)->shared->datarefs
#define BB_VINFO_DDRS(B)             (B)->shared->ddrs

/*-----------------------------------------------------------------*/
/* Info on vectorized defs.                                        */
/*-----------------------------------------------------------------*/
enum stmt_vec_info_type {
  undef_vec_info_type = 0,
  load_vec_info_type,
  store_vec_info_type,
  shift_vec_info_type,
  op_vec_info_type,
  call_vec_info_type,
  call_simd_clone_vec_info_type,
  assignment_vec_info_type,
  condition_vec_info_type,
  comparison_vec_info_type,
  reduc_vec_info_type,
  induc_vec_info_type,
  type_promotion_vec_info_type,
  type_demotion_vec_info_type,
  type_conversion_vec_info_type,
  cycle_phi_info_type,
  lc_phi_info_type,
  phi_info_type,
  loop_exit_ctrl_vec_info_type
};

/* Indicates whether/how a variable is used in the scope of loop/basic
   block.  */
enum vect_relevant {
  vect_unused_in_scope = 0,

  /* The def is only used outside the loop.  */
  vect_used_only_live,
  /* The def is in the inner loop, and the use is in the outer loop, and the
     use is a reduction stmt.  */
  vect_used_in_outer_by_reduction,
  /* The def is in the inner loop, and the use is in the outer loop (and is
     not part of reduction).  */
  vect_used_in_outer,

  /* defs that feed computations that end up (only) in a reduction. These
     defs may be used by non-reduction stmts, but eventually, any
     computations/values that are affected by these defs are used to compute
     a reduction (i.e. don't get stored to memory, for example). We use this
     to identify computations that we can change the order in which they are
     computed.  */
  vect_used_by_reduction,

  vect_used_in_scope
};

/* The type of vectorization that can be applied to the stmt: regular loop-based
   vectorization; pure SLP - the stmt is a part of SLP instances and does not
   have uses outside SLP instances; or hybrid SLP and loop-based - the stmt is
   a part of SLP instance and also must be loop-based vectorized, since it has
   uses outside SLP sequences.

   In the loop context the meanings of pure and hybrid SLP are slightly
   different. By saying that pure SLP is applied to the loop, we mean that we
   exploit only intra-iteration parallelism in the loop; i.e., the loop can be
   vectorized without doing any conceptual unrolling, cause we don't pack
   together stmts from different iterations, only within a single iteration.
   Loop hybrid SLP means that we exploit both intra-iteration and
   inter-iteration parallelism (e.g., number of elements in the vector is 4
   and the slp-group-size is 2, in which case we don't have enough parallelism
   within an iteration, so we obtain the rest of the parallelism from subsequent
   iterations by unrolling the loop by 2).  */
enum slp_vect_type {
  loop_vect = 0,
  pure_slp,
  hybrid
};

/* Says whether a statement is a load, a store of a vectorized statement
   result, or a store of an invariant value.  */
enum vec_load_store_type {
  VLS_LOAD,
  VLS_STORE,
  VLS_STORE_INVARIANT
};

/* Describes how we're going to vectorize an individual load or store,
   or a group of loads or stores.  */
enum vect_memory_access_type {
  /* An access to an invariant address.  This is used only for loads.  */
  VMAT_INVARIANT,

  /* A simple contiguous access.  */
  VMAT_CONTIGUOUS,

  /* A contiguous access that goes down in memory rather than up,
     with no additional permutation.  This is used only for stores
     of invariants.  */
  VMAT_CONTIGUOUS_DOWN,

  /* A simple contiguous access in which the elements need to be permuted
     after loading or before storing.  Only used for loop vectorization;
     SLP uses separate permutes.  */
  VMAT_CONTIGUOUS_PERMUTE,

  /* A simple contiguous access in which the elements need to be reversed
     after loading or before storing.  */
  VMAT_CONTIGUOUS_REVERSE,

  /* An access that uses IFN_LOAD_LANES or IFN_STORE_LANES.  */
  VMAT_LOAD_STORE_LANES,

  /* An access in which each scalar element is loaded or stored
     individually.  */
  VMAT_ELEMENTWISE,

  /* A hybrid of VMAT_CONTIGUOUS and VMAT_ELEMENTWISE, used for grouped
     SLP accesses.  Each unrolled iteration uses a contiguous load
     or store for the whole group, but the groups from separate iterations
     are combined in the same way as for VMAT_ELEMENTWISE.  */
  VMAT_STRIDED_SLP,

  /* The access uses gather loads or scatter stores.  */
  VMAT_GATHER_SCATTER
};

class dr_vec_info {
public:
  /* The data reference itself.  */
  data_reference *dr;
  /* The statement that contains the data reference.  */
  stmt_vec_info stmt;
  /* The analysis group this DR belongs to when doing BB vectorization.
     DRs of the same group belong to the same conditional execution context.  */
  unsigned group;
  /* The misalignment in bytes of the reference, or -1 if not known.  */
  int misalignment;
  /* The byte alignment that we'd ideally like the reference to have,
     and the value that misalignment is measured against.  */
  poly_uint64 target_alignment;
  /* If true the alignment of base_decl needs to be increased.  */
  bool base_misaligned;
  tree base_decl;

  /* Stores current vectorized loop's offset.  To be added to the DR's
     offset to calculate current offset of data reference.  */
  tree offset;
};

typedef struct data_reference *dr_p;

class _stmt_vec_info {
public:

  enum stmt_vec_info_type type;

  /* Indicates whether this stmts is part of a computation whose result is
     used outside the loop.  */
  bool live;

  /* Stmt is part of some pattern (computation idiom)  */
  bool in_pattern_p;

  /* True if the statement was created during pattern recognition as
     part of the replacement for RELATED_STMT.  This implies that the
     statement isn't part of any basic block, although for convenience
     its gimple_bb is the same as for RELATED_STMT.  */
  bool pattern_stmt_p;

  /* Is this statement vectorizable or should it be skipped in (partial)
     vectorization.  */
  bool vectorizable;

  /* The stmt to which this info struct refers to.  */
  gimple *stmt;

  /* The vector type to be used for the LHS of this statement.  */
  tree vectype;

  /* The vectorized stmts.  */
  vec<gimple *> vec_stmts;

  /* The following is relevant only for stmts that contain a non-scalar
     data-ref (array/pointer/struct access). A GIMPLE stmt is expected to have
     at most one such data-ref.  */

  dr_vec_info dr_aux;

  /* Information about the data-ref relative to this loop
     nest (the loop that is being considered for vectorization).  */
  innermost_loop_behavior dr_wrt_vec_loop;

  /* For loop PHI nodes, the base and evolution part of it.  This makes sure
     this information is still available in vect_update_ivs_after_vectorizer
     where we may not be able to re-analyze the PHI nodes evolution as
     peeling for the prologue loop can make it unanalyzable.  The evolution
     part is still correct after peeling, but the base may have changed from
     the version here.  */
  tree loop_phi_evolution_base_unchanged;
  tree loop_phi_evolution_part;

  /* Used for various bookkeeping purposes, generally holding a pointer to
     some other stmt S that is in some way "related" to this stmt.
     Current use of this field is:
        If this stmt is part of a pattern (i.e. the field 'in_pattern_p' is
        true): S is the "pattern stmt" that represents (and replaces) the
        sequence of stmts that constitutes the pattern.  Similarly, the
        related_stmt of the "pattern stmt" points back to this stmt (which is
        the last stmt in the original sequence of stmts that constitutes the
        pattern).  */
  stmt_vec_info related_stmt;

  /* Used to keep a sequence of def stmts of a pattern stmt if such exists.
     The sequence is attached to the original statement rather than the
     pattern statement.  */
  gimple_seq pattern_def_seq;

  /* Selected SIMD clone's function info.  First vector element
     is SIMD clone's function decl, followed by a pair of trees (base + step)
     for linear arguments (pair of NULLs for other arguments).  */
  vec<tree> simd_clone_info;

  /* Classify the def of this stmt.  */
  enum vect_def_type def_type;

  /*  Whether the stmt is SLPed, loop-based vectorized, or both.  */
  enum slp_vect_type slp_type;

  /* Interleaving and reduction chains info.  */
  /* First element in the group.  */
  stmt_vec_info first_element;
  /* Pointer to the next element in the group.  */
  stmt_vec_info next_element;
  /* The size of the group.  */
  unsigned int size;
  /* For stores, number of stores from this group seen. We vectorize the last
     one.  */
  unsigned int store_count;
  /* For loads only, the gap from the previous load. For consecutive loads, GAP
     is 1.  */
  unsigned int gap;

  /* The minimum negative dependence distance this stmt participates in
     or zero if none.  */
  unsigned int min_neg_dist;

  /* Not all stmts in the loop need to be vectorized. e.g, the increment
     of the loop induction variable and computation of array indexes. relevant
     indicates whether the stmt needs to be vectorized.  */
  enum vect_relevant relevant;

  /* For loads if this is a gather, for stores if this is a scatter.  */
  bool gather_scatter_p;

  /* True if this is an access with loop-invariant stride.  */
  bool strided_p;

  /* For both loads and stores.  */
  unsigned simd_lane_access_p : 3;

  /* Classifies how the load or store is going to be implemented
     for loop vectorization.  */
  vect_memory_access_type memory_access_type;

  /* For INTEGER_INDUC_COND_REDUCTION, the initial value to be used.  */
  tree induc_cond_initial_val;

  /* If not NULL the value to be added to compute final reduction value.  */
  tree reduc_epilogue_adjustment;

  /* On a reduction PHI the reduction type as detected by
     vect_is_simple_reduction and vectorizable_reduction.  */
  enum vect_reduction_type reduc_type;

  /* The original reduction code, to be used in the epilogue.  */
  code_helper reduc_code;
  /* An internal function we should use in the epilogue.  */
  internal_fn reduc_fn;

  /* On a stmt participating in the reduction the index of the operand
     on the reduction SSA cycle.  */
  int reduc_idx;

  /* On a reduction PHI the def returned by vect_force_simple_reduction.
     On the def returned by vect_force_simple_reduction the
     corresponding PHI.  */
  stmt_vec_info reduc_def;

  /* The vector input type relevant for reduction vectorization.  */
  tree reduc_vectype_in;

  /* The vector type for performing the actual reduction.  */
  tree reduc_vectype;

  /* If IS_REDUC_INFO is true and if the vector code is performing
     N scalar reductions in parallel, this variable gives the initial
     scalar values of those N reductions.  */
  vec<tree> reduc_initial_values;

  /* If IS_REDUC_INFO is true and if the vector code is performing
     N scalar reductions in parallel, this variable gives the vectorized code's
     final (scalar) result for each of those N reductions.  In other words,
     REDUC_SCALAR_RESULTS[I] replaces the original scalar code's loop-closed
     SSA PHI for reduction number I.  */
  vec<tree> reduc_scalar_results;

  /* Only meaningful if IS_REDUC_INFO.  If non-null, the reduction is
     being performed by an epilogue loop and we have decided to reuse
     this accumulator from the main loop.  */
  vect_reusable_accumulator *reused_accumulator;

  /* Whether we force a single cycle PHI during reduction vectorization.  */
  bool force_single_cycle;

  /* Whether on this stmt reduction meta is recorded.  */
  bool is_reduc_info;

  /* If nonzero, the lhs of the statement could be truncated to this
     many bits without affecting any users of the result.  */
  unsigned int min_output_precision;

  /* If nonzero, all non-boolean input operands have the same precision,
     and they could each be truncated to this many bits without changing
     the result.  */
  unsigned int min_input_precision;

  /* If OPERATION_BITS is nonzero, the statement could be performed on
     an integer with the sign and number of bits given by OPERATION_SIGN
     and OPERATION_BITS without changing the result.  */
  unsigned int operation_precision;
  signop operation_sign;

  /* If the statement produces a boolean result, this value describes
     how we should choose the associated vector type.  The possible
     values are:

     - an integer precision N if we should use the vector mask type
       associated with N-bit integers.  This is only used if all relevant
       input booleans also want the vector mask type for N-bit integers,
       or if we can convert them into that form by pattern-matching.

     - ~0U if we considered choosing a vector mask type but decided
       to treat the boolean as a normal integer type instead.

     - 0 otherwise.  This means either that the operation isn't one that
       could have a vector mask type (and so should have a normal vector
       type instead) or that we simply haven't made a choice either way.  */
  unsigned int mask_precision;

  /* True if this is only suitable for SLP vectorization.  */
  bool slp_vect_only_p;

  /* True if this is a pattern that can only be handled by SLP
     vectorization.  */
  bool slp_vect_pattern_only_p;
};

/* Information about a gather/scatter call.  */
struct gather_scatter_info {
  /* The internal function to use for the gather/scatter operation,
     or IFN_LAST if a built-in function should be used instead.  */
  internal_fn ifn;

  /* The FUNCTION_DECL for the built-in gather/scatter function,
     or null if an internal function should be used instead.  */
  tree decl;

  /* The loop-invariant base value.  */
  tree base;

  /* The original scalar offset, which is a non-loop-invariant SSA_NAME.  */
  tree offset;

  /* Each offset element should be multiplied by this amount before
     being added to the base.  */
  int scale;

  /* The definition type for the vectorized offset.  */
  enum vect_def_type offset_dt;

  /* The type of the vectorized offset.  */
  tree offset_vectype;

  /* The type of the scalar elements after loading or before storing.  */
  tree element_type;

  /* The type of the scalar elements being loaded or stored.  */
  tree memory_type;
};

/* Access Functions.  */
#define STMT_VINFO_TYPE(S)                 (S)->type
#define STMT_VINFO_STMT(S)                 (S)->stmt
#define STMT_VINFO_RELEVANT(S)             (S)->relevant
#define STMT_VINFO_LIVE_P(S)               (S)->live
#define STMT_VINFO_VECTYPE(S)              (S)->vectype
#define STMT_VINFO_VEC_STMTS(S)            (S)->vec_stmts
#define STMT_VINFO_VECTORIZABLE(S)         (S)->vectorizable
#define STMT_VINFO_DATA_REF(S)             ((S)->dr_aux.dr + 0)
#define STMT_VINFO_GATHER_SCATTER_P(S)	   (S)->gather_scatter_p
#define STMT_VINFO_STRIDED_P(S)	   	   (S)->strided_p
#define STMT_VINFO_MEMORY_ACCESS_TYPE(S)   (S)->memory_access_type
#define STMT_VINFO_SIMD_LANE_ACCESS_P(S)   (S)->simd_lane_access_p
#define STMT_VINFO_VEC_INDUC_COND_INITIAL_VAL(S) (S)->induc_cond_initial_val
#define STMT_VINFO_REDUC_EPILOGUE_ADJUSTMENT(S) (S)->reduc_epilogue_adjustment
#define STMT_VINFO_REDUC_IDX(S)		   (S)->reduc_idx
#define STMT_VINFO_FORCE_SINGLE_CYCLE(S)   (S)->force_single_cycle

#define STMT_VINFO_DR_WRT_VEC_LOOP(S)      (S)->dr_wrt_vec_loop
#define STMT_VINFO_DR_BASE_ADDRESS(S)      (S)->dr_wrt_vec_loop.base_address
#define STMT_VINFO_DR_INIT(S)              (S)->dr_wrt_vec_loop.init
#define STMT_VINFO_DR_OFFSET(S)            (S)->dr_wrt_vec_loop.offset
#define STMT_VINFO_DR_STEP(S)              (S)->dr_wrt_vec_loop.step
#define STMT_VINFO_DR_BASE_ALIGNMENT(S)    (S)->dr_wrt_vec_loop.base_alignment
#define STMT_VINFO_DR_BASE_MISALIGNMENT(S) \
  (S)->dr_wrt_vec_loop.base_misalignment
#define STMT_VINFO_DR_OFFSET_ALIGNMENT(S) \
  (S)->dr_wrt_vec_loop.offset_alignment
#define STMT_VINFO_DR_STEP_ALIGNMENT(S) \
  (S)->dr_wrt_vec_loop.step_alignment

#define STMT_VINFO_DR_INFO(S) \
  (gcc_checking_assert ((S)->dr_aux.stmt == (S)), &(S)->dr_aux)

#define STMT_VINFO_IN_PATTERN_P(S)         (S)->in_pattern_p
#define STMT_VINFO_RELATED_STMT(S)         (S)->related_stmt
#define STMT_VINFO_PATTERN_DEF_SEQ(S)      (S)->pattern_def_seq
#define STMT_VINFO_SIMD_CLONE_INFO(S)	   (S)->simd_clone_info
#define STMT_VINFO_DEF_TYPE(S)             (S)->def_type
#define STMT_VINFO_GROUPED_ACCESS(S) \
  ((S)->dr_aux.dr && DR_GROUP_FIRST_ELEMENT(S))
#define STMT_VINFO_LOOP_PHI_EVOLUTION_BASE_UNCHANGED(S) (S)->loop_phi_evolution_base_unchanged
#define STMT_VINFO_LOOP_PHI_EVOLUTION_PART(S) (S)->loop_phi_evolution_part
#define STMT_VINFO_MIN_NEG_DIST(S)	(S)->min_neg_dist
#define STMT_VINFO_REDUC_TYPE(S)	(S)->reduc_type
#define STMT_VINFO_REDUC_CODE(S)	(S)->reduc_code
#define STMT_VINFO_REDUC_FN(S)		(S)->reduc_fn
#define STMT_VINFO_REDUC_DEF(S)		(S)->reduc_def
#define STMT_VINFO_REDUC_VECTYPE(S)     (S)->reduc_vectype
#define STMT_VINFO_REDUC_VECTYPE_IN(S)  (S)->reduc_vectype_in
#define STMT_VINFO_SLP_VECT_ONLY(S)     (S)->slp_vect_only_p
#define STMT_VINFO_SLP_VECT_ONLY_PATTERN(S) (S)->slp_vect_pattern_only_p

#define DR_GROUP_FIRST_ELEMENT(S) \
  (gcc_checking_assert ((S)->dr_aux.dr), (S)->first_element)
#define DR_GROUP_NEXT_ELEMENT(S) \
  (gcc_checking_assert ((S)->dr_aux.dr), (S)->next_element)
#define DR_GROUP_SIZE(S) \
  (gcc_checking_assert ((S)->dr_aux.dr), (S)->size)
#define DR_GROUP_STORE_COUNT(S) \
  (gcc_checking_assert ((S)->dr_aux.dr), (S)->store_count)
#define DR_GROUP_GAP(S) \
  (gcc_checking_assert ((S)->dr_aux.dr), (S)->gap)

#define REDUC_GROUP_FIRST_ELEMENT(S) \
  (gcc_checking_assert (!(S)->dr_aux.dr), (S)->first_element)
#define REDUC_GROUP_NEXT_ELEMENT(S) \
  (gcc_checking_assert (!(S)->dr_aux.dr), (S)->next_element)
#define REDUC_GROUP_SIZE(S) \
  (gcc_checking_assert (!(S)->dr_aux.dr), (S)->size)

#define STMT_VINFO_RELEVANT_P(S)          ((S)->relevant != vect_unused_in_scope)

#define HYBRID_SLP_STMT(S)                ((S)->slp_type == hybrid)
#define PURE_SLP_STMT(S)                  ((S)->slp_type == pure_slp)
#define STMT_SLP_TYPE(S)                   (S)->slp_type

/* Contains the scalar or vector costs for a vec_info.  */
class vector_costs
{
public:
  vector_costs (vec_info *, bool);
  virtual ~vector_costs () {}

  /* Update the costs in response to adding COUNT copies of a statement.

     - WHERE specifies whether the cost occurs in the loop prologue,
       the loop body, or the loop epilogue.
     - KIND is the kind of statement, which is always meaningful.
     - STMT_INFO or NODE, if nonnull, describe the statement that will be
       vectorized.
     - VECTYPE, if nonnull, is the vector type that the vectorized
       statement will operate on.  Note that this should be used in
       preference to STMT_VINFO_VECTYPE (STMT_INFO) since the latter
       is not correct for SLP.
     - for unaligned_load and unaligned_store statements, MISALIGN is
       the byte misalignment of the load or store relative to the target's
       preferred alignment for VECTYPE, or DR_MISALIGNMENT_UNKNOWN
       if the misalignment is not known.

     Return the calculated cost as well as recording it.  The return
     value is used for dumping purposes.  */
  virtual unsigned int add_stmt_cost (int count, vect_cost_for_stmt kind,
				      stmt_vec_info stmt_info,
				      slp_tree node,
				      tree vectype, int misalign,
				      vect_cost_model_location where);

  /* Finish calculating the cost of the code.  The results can be
     read back using the functions below.

     If the costs describe vector code, SCALAR_COSTS gives the costs
     of the corresponding scalar code, otherwise it is null.  */
  virtual void finish_cost (const vector_costs *scalar_costs);

  /* The costs in THIS and OTHER both describe ways of vectorizing
     a main loop.  Return true if the costs described by THIS are
     cheaper than the costs described by OTHER.  Return false if any
     of the following are true:

     - THIS and OTHER are of equal cost
     - OTHER is better than THIS
     - we can't be sure about the relative costs of THIS and OTHER.  */
  virtual bool better_main_loop_than_p (const vector_costs *other) const;

  /* Likewise, but the costs in THIS and OTHER both describe ways of
     vectorizing an epilogue loop of MAIN_LOOP.  */
  virtual bool better_epilogue_loop_than_p (const vector_costs *other,
					    loop_vec_info main_loop) const;

  unsigned int prologue_cost () const;
  unsigned int body_cost () const;
  unsigned int epilogue_cost () const;
  unsigned int outside_cost () const;
  unsigned int total_cost () const;
  unsigned int suggested_unroll_factor () const;

protected:
  unsigned int record_stmt_cost (stmt_vec_info, vect_cost_model_location,
				 unsigned int);
  unsigned int adjust_cost_for_freq (stmt_vec_info, vect_cost_model_location,
				     unsigned int);
  int compare_inside_loop_cost (const vector_costs *) const;
  int compare_outside_loop_cost (const vector_costs *) const;

  /* The region of code that we're considering vectorizing.  */
  vec_info *m_vinfo;

  /* True if we're costing the scalar code, false if we're costing
     the vector code.  */
  bool m_costing_for_scalar;

  /* The costs of the three regions, indexed by vect_cost_model_location.  */
  unsigned int m_costs[3];

  /* The suggested unrolling factor determined at finish_cost.  */
  unsigned int m_suggested_unroll_factor;

  /* True if finish_cost has been called.  */
  bool m_finished;
};

/* Create costs for VINFO.  COSTING_FOR_SCALAR is true if the costs
   are for scalar code, false if they are for vector code.  */

inline
vector_costs::vector_costs (vec_info *vinfo, bool costing_for_scalar)
  : m_vinfo (vinfo),
    m_costing_for_scalar (costing_for_scalar),
    m_costs (),
    m_suggested_unroll_factor(1),
    m_finished (false)
{
}

/* Return the cost of the prologue code (in abstract units).  */

inline unsigned int
vector_costs::prologue_cost () const
{
  gcc_checking_assert (m_finished);
  return m_costs[vect_prologue];
}

/* Return the cost of the body code (in abstract units).  */

inline unsigned int
vector_costs::body_cost () const
{
  gcc_checking_assert (m_finished);
  return m_costs[vect_body];
}

/* Return the cost of the epilogue code (in abstract units).  */

inline unsigned int
vector_costs::epilogue_cost () const
{
  gcc_checking_assert (m_finished);
  return m_costs[vect_epilogue];
}

/* Return the cost of the prologue and epilogue code (in abstract units).  */

inline unsigned int
vector_costs::outside_cost () const
{
  return prologue_cost () + epilogue_cost ();
}

/* Return the cost of the prologue, body and epilogue code
   (in abstract units).  */

inline unsigned int
vector_costs::total_cost () const
{
  return body_cost () + outside_cost ();
}

/* Return the suggested unroll factor.  */

inline unsigned int
vector_costs::suggested_unroll_factor () const
{
  gcc_checking_assert (m_finished);
  return m_suggested_unroll_factor;
}

#define VECT_MAX_COST 1000

/* The maximum number of intermediate steps required in multi-step type
   conversion.  */
#define MAX_INTERM_CVT_STEPS         3

#define MAX_VECTORIZATION_FACTOR INT_MAX

/* Nonzero if TYPE represents a (scalar) boolean type or type
   in the middle-end compatible with it (unsigned precision 1 integral
   types).  Used to determine which types should be vectorized as
   VECTOR_BOOLEAN_TYPE_P.  */

#define VECT_SCALAR_BOOLEAN_TYPE_P(TYPE) \
  (TREE_CODE (TYPE) == BOOLEAN_TYPE		\
   || ((TREE_CODE (TYPE) == INTEGER_TYPE	\
	|| TREE_CODE (TYPE) == ENUMERAL_TYPE)	\
       && TYPE_PRECISION (TYPE) == 1		\
       && TYPE_UNSIGNED (TYPE)))

static inline bool
nested_in_vect_loop_p (class loop *loop, stmt_vec_info stmt_info)
{
  return (loop->inner
	  && (loop->inner == (gimple_bb (stmt_info->stmt))->loop_father));
}

/* PHI is either a scalar reduction phi or a scalar induction phi.
   Return the initial value of the variable on entry to the containing
   loop.  */

static inline tree
vect_phi_initial_value (gphi *phi)
{
  basic_block bb = gimple_bb (phi);
  edge pe = loop_preheader_edge (bb->loop_father);
  gcc_assert (pe->dest == bb);
  return PHI_ARG_DEF_FROM_EDGE (phi, pe);
}

/* Return true if STMT_INFO should produce a vector mask type rather than
   a normal nonmask type.  */

static inline bool
vect_use_mask_type_p (stmt_vec_info stmt_info)
{
  return stmt_info->mask_precision && stmt_info->mask_precision != ~0U;
}

/* Return TRUE if a statement represented by STMT_INFO is a part of a
   pattern.  */

static inline bool
is_pattern_stmt_p (stmt_vec_info stmt_info)
{
  return stmt_info->pattern_stmt_p;
}

/* If STMT_INFO is a pattern statement, return the statement that it
   replaces, otherwise return STMT_INFO itself.  */

inline stmt_vec_info
vect_orig_stmt (stmt_vec_info stmt_info)
{
  if (is_pattern_stmt_p (stmt_info))
    return STMT_VINFO_RELATED_STMT (stmt_info);
  return stmt_info;
}

/* Return the later statement between STMT1_INFO and STMT2_INFO.  */

static inline stmt_vec_info
get_later_stmt (stmt_vec_info stmt1_info, stmt_vec_info stmt2_info)
{
  if (gimple_uid (vect_orig_stmt (stmt1_info)->stmt)
      > gimple_uid (vect_orig_stmt (stmt2_info)->stmt))
    return stmt1_info;
  else
    return stmt2_info;
}

/* If STMT_INFO has been replaced by a pattern statement, return the
   replacement statement, otherwise return STMT_INFO itself.  */

inline stmt_vec_info
vect_stmt_to_vectorize (stmt_vec_info stmt_info)
{
  if (STMT_VINFO_IN_PATTERN_P (stmt_info))
    return STMT_VINFO_RELATED_STMT (stmt_info);
  return stmt_info;
}

/* Return true if BB is a loop header.  */

static inline bool
is_loop_header_bb_p (basic_block bb)
{
  if (bb == (bb->loop_father)->header)
    return true;
  gcc_checking_assert (EDGE_COUNT (bb->preds) == 1);
  return false;
}

/* Return pow2 (X).  */

static inline int
vect_pow2 (int x)
{
  int i, res = 1;

  for (i = 0; i < x; i++)
    res *= 2;

  return res;
}

/* Alias targetm.vectorize.builtin_vectorization_cost.  */

static inline int
builtin_vectorization_cost (enum vect_cost_for_stmt type_of_cost,
			    tree vectype, int misalign)
{
  return targetm.vectorize.builtin_vectorization_cost (type_of_cost,
						       vectype, misalign);
}

/* Get cost by calling cost target builtin.  */

static inline
int vect_get_stmt_cost (enum vect_cost_for_stmt type_of_cost)
{
  return builtin_vectorization_cost (type_of_cost, NULL, 0);
}

/* Alias targetm.vectorize.init_cost.  */

static inline vector_costs *
init_cost (vec_info *vinfo, bool costing_for_scalar)
{
  return targetm.vectorize.create_costs (vinfo, costing_for_scalar);
}

extern void dump_stmt_cost (FILE *, int, enum vect_cost_for_stmt,
			    stmt_vec_info, slp_tree, tree, int, unsigned,
			    enum vect_cost_model_location);

/* Alias targetm.vectorize.add_stmt_cost.  */

static inline unsigned
add_stmt_cost (vector_costs *costs, int count,
	       enum vect_cost_for_stmt kind,
	       stmt_vec_info stmt_info, slp_tree node,
	       tree vectype, int misalign,
	       enum vect_cost_model_location where)
{
  unsigned cost = costs->add_stmt_cost (count, kind, stmt_info, node, vectype,
					misalign, where);
  if (dump_file && (dump_flags & TDF_DETAILS))
    dump_stmt_cost (dump_file, count, kind, stmt_info, node, vectype, misalign,
		    cost, where);
  return cost;
}

static inline unsigned
add_stmt_cost (vector_costs *costs, int count, enum vect_cost_for_stmt kind,
	       enum vect_cost_model_location where)
{
  gcc_assert (kind == cond_branch_taken || kind == cond_branch_not_taken
	      || kind == scalar_stmt);
  return add_stmt_cost (costs, count, kind, NULL, NULL, NULL_TREE, 0, where);
}

/* Alias targetm.vectorize.add_stmt_cost.  */

static inline unsigned
add_stmt_cost (vector_costs *costs, stmt_info_for_cost *i)
{
  return add_stmt_cost (costs, i->count, i->kind, i->stmt_info, i->node,
			i->vectype, i->misalign, i->where);
}

/* Alias targetm.vectorize.finish_cost.  */

static inline void
finish_cost (vector_costs *costs, const vector_costs *scalar_costs,
	     unsigned *prologue_cost, unsigned *body_cost,
	     unsigned *epilogue_cost, unsigned *suggested_unroll_factor = NULL)
{
  costs->finish_cost (scalar_costs);
  *prologue_cost = costs->prologue_cost ();
  *body_cost = costs->body_cost ();
  *epilogue_cost = costs->epilogue_cost ();
  if (suggested_unroll_factor)
    *suggested_unroll_factor = costs->suggested_unroll_factor ();
}

inline void
add_stmt_costs (vector_costs *costs, stmt_vector_for_cost *cost_vec)
{
  stmt_info_for_cost *cost;
  unsigned i;
  FOR_EACH_VEC_ELT (*cost_vec, i, cost)
    add_stmt_cost (costs, cost->count, cost->kind, cost->stmt_info,
		   cost->node, cost->vectype, cost->misalign, cost->where);
}

/*-----------------------------------------------------------------*/
/* Info on data references alignment.                              */
/*-----------------------------------------------------------------*/
#define DR_MISALIGNMENT_UNKNOWN (-1)
#define DR_MISALIGNMENT_UNINITIALIZED (-2)

inline void
set_dr_misalignment (dr_vec_info *dr_info, int val)
{
  dr_info->misalignment = val;
}

extern int dr_misalignment (dr_vec_info *dr_info, tree vectype,
			    poly_int64 offset = 0);

#define SET_DR_MISALIGNMENT(DR, VAL) set_dr_misalignment (DR, VAL)

/* Only defined once DR_MISALIGNMENT is defined.  */
static inline const poly_uint64
dr_target_alignment (dr_vec_info *dr_info)
{
  if (STMT_VINFO_GROUPED_ACCESS (dr_info->stmt))
    dr_info = STMT_VINFO_DR_INFO (DR_GROUP_FIRST_ELEMENT (dr_info->stmt));
  return dr_info->target_alignment;
}
#define DR_TARGET_ALIGNMENT(DR) dr_target_alignment (DR)

static inline void
set_dr_target_alignment (dr_vec_info *dr_info, poly_uint64 val)
{
  dr_info->target_alignment = val;
}
#define SET_DR_TARGET_ALIGNMENT(DR, VAL) set_dr_target_alignment (DR, VAL)

/* Return true if data access DR_INFO is aligned to the targets
   preferred alignment for VECTYPE (which may be less than a full vector).  */

static inline bool
aligned_access_p (dr_vec_info *dr_info, tree vectype)
{
  return (dr_misalignment (dr_info, vectype) == 0);
}

/* Return TRUE if the (mis-)alignment of the data access is known with
   respect to the targets preferred alignment for VECTYPE, and FALSE
   otherwise.  */

static inline bool
known_alignment_for_access_p (dr_vec_info *dr_info, tree vectype)
{
  return (dr_misalignment (dr_info, vectype) != DR_MISALIGNMENT_UNKNOWN);
}

/* Return the minimum alignment in bytes that the vectorized version
   of DR_INFO is guaranteed to have.  */

static inline unsigned int
vect_known_alignment_in_bytes (dr_vec_info *dr_info, tree vectype)
{
  int misalignment = dr_misalignment (dr_info, vectype);
  if (misalignment == DR_MISALIGNMENT_UNKNOWN)
    return TYPE_ALIGN_UNIT (TREE_TYPE (DR_REF (dr_info->dr)));
  else if (misalignment == 0)
    return known_alignment (DR_TARGET_ALIGNMENT (dr_info));
  return misalignment & -misalignment;
}

/* Return the behavior of DR_INFO with respect to the vectorization context
   (which for outer loop vectorization might not be the behavior recorded
   in DR_INFO itself).  */

static inline innermost_loop_behavior *
vect_dr_behavior (vec_info *vinfo, dr_vec_info *dr_info)
{
  stmt_vec_info stmt_info = dr_info->stmt;
  loop_vec_info loop_vinfo = dyn_cast<loop_vec_info> (vinfo);
  if (loop_vinfo == NULL
      || !nested_in_vect_loop_p (LOOP_VINFO_LOOP (loop_vinfo), stmt_info))
    return &DR_INNERMOST (dr_info->dr);
  else
    return &STMT_VINFO_DR_WRT_VEC_LOOP (stmt_info);
}

/* Return the offset calculated by adding the offset of this DR_INFO to the
   corresponding data_reference's offset.  If CHECK_OUTER then use
   vect_dr_behavior to select the appropriate data_reference to use.  */

inline tree
get_dr_vinfo_offset (vec_info *vinfo,
		     dr_vec_info *dr_info, bool check_outer = false)
{
  innermost_loop_behavior *base;
  if (check_outer)
    base = vect_dr_behavior (vinfo, dr_info);
  else
    base = &dr_info->dr->innermost;

  tree offset = base->offset;

  if (!dr_info->offset)
    return offset;

  offset = fold_convert (sizetype, offset);
  return fold_build2 (PLUS_EXPR, TREE_TYPE (dr_info->offset), offset,
		      dr_info->offset);
}


/* Return the vect cost model for LOOP.  */
static inline enum vect_cost_model
loop_cost_model (loop_p loop)
{
  if (loop != NULL
      && loop->force_vectorize
      && flag_simd_cost_model != VECT_COST_MODEL_DEFAULT)
    return flag_simd_cost_model;
  return flag_vect_cost_model;
}

/* Return true if the vect cost model is unlimited.  */
static inline bool
unlimited_cost_model (loop_p loop)
{
  return loop_cost_model (loop) == VECT_COST_MODEL_UNLIMITED;
}

/* Return true if the loop described by LOOP_VINFO is fully-masked and
   if the first iteration should use a partial mask in order to achieve
   alignment.  */

static inline bool
vect_use_loop_mask_for_alignment_p (loop_vec_info loop_vinfo)
{
  return (LOOP_VINFO_FULLY_MASKED_P (loop_vinfo)
	  && LOOP_VINFO_PEELING_FOR_ALIGNMENT (loop_vinfo));
}

/* Return the number of vectors of type VECTYPE that are needed to get
   NUNITS elements.  NUNITS should be based on the vectorization factor,
   so it is always a known multiple of the number of elements in VECTYPE.  */

static inline unsigned int
vect_get_num_vectors (poly_uint64 nunits, tree vectype)
{
  return exact_div (nunits, TYPE_VECTOR_SUBPARTS (vectype)).to_constant ();
}

/* Return the number of copies needed for loop vectorization when
   a statement operates on vectors of type VECTYPE.  This is the
   vectorization factor divided by the number of elements in
   VECTYPE and is always known at compile time.  */

static inline unsigned int
vect_get_num_copies (loop_vec_info loop_vinfo, tree vectype)
{
  return vect_get_num_vectors (LOOP_VINFO_VECT_FACTOR (loop_vinfo), vectype);
}

/* Update maximum unit count *MAX_NUNITS so that it accounts for
   NUNITS.  *MAX_NUNITS can be 1 if we haven't yet recorded anything.  */

static inline void
vect_update_max_nunits (poly_uint64 *max_nunits, poly_uint64 nunits)
{
  /* All unit counts have the form vec_info::vector_size * X for some
     rational X, so two unit sizes must have a common multiple.
     Everything is a multiple of the initial value of 1.  */
  *max_nunits = force_common_multiple (*max_nunits, nunits);
}

/* Update maximum unit count *MAX_NUNITS so that it accounts for
   the number of units in vector type VECTYPE.  *MAX_NUNITS can be 1
   if we haven't yet recorded any vector types.  */

static inline void
vect_update_max_nunits (poly_uint64 *max_nunits, tree vectype)
{
  vect_update_max_nunits (max_nunits, TYPE_VECTOR_SUBPARTS (vectype));
}

/* Return the vectorization factor that should be used for costing
   purposes while vectorizing the loop described by LOOP_VINFO.
   Pick a reasonable estimate if the vectorization factor isn't
   known at compile time.  */

static inline unsigned int
vect_vf_for_cost (loop_vec_info loop_vinfo)
{
  return estimated_poly_value (LOOP_VINFO_VECT_FACTOR (loop_vinfo));
}

/* Estimate the number of elements in VEC_TYPE for costing purposes.
   Pick a reasonable estimate if the exact number isn't known at
   compile time.  */

static inline unsigned int
vect_nunits_for_cost (tree vec_type)
{
  return estimated_poly_value (TYPE_VECTOR_SUBPARTS (vec_type));
}

/* Return the maximum possible vectorization factor for LOOP_VINFO.  */

static inline unsigned HOST_WIDE_INT
vect_max_vf (loop_vec_info loop_vinfo)
{
  unsigned HOST_WIDE_INT vf;
  if (LOOP_VINFO_VECT_FACTOR (loop_vinfo).is_constant (&vf))
    return vf;
  return MAX_VECTORIZATION_FACTOR;
}

/* Return the size of the value accessed by unvectorized data reference
   DR_INFO.  This is only valid once STMT_VINFO_VECTYPE has been calculated
   for the associated gimple statement, since that guarantees that DR_INFO
   accesses either a scalar or a scalar equivalent.  ("Scalar equivalent"
   here includes things like V1SI, which can be vectorized in the same way
   as a plain SI.)  */

inline unsigned int
vect_get_scalar_dr_size (dr_vec_info *dr_info)
{
  return tree_to_uhwi (TYPE_SIZE_UNIT (TREE_TYPE (DR_REF (dr_info->dr))));
}

/* Return true if LOOP_VINFO requires a runtime check for whether the
   vector loop is profitable.  */

inline bool
vect_apply_runtime_profitability_check_p (loop_vec_info loop_vinfo)
{
  unsigned int th = LOOP_VINFO_COST_MODEL_THRESHOLD (loop_vinfo);
  return (!LOOP_VINFO_NITERS_KNOWN_P (loop_vinfo)
	  && th >= vect_vf_for_cost (loop_vinfo));
}

/* Source location + hotness information. */
extern dump_user_location_t vect_location;

/* A macro for calling:
     dump_begin_scope (MSG, vect_location);
   via an RAII object, thus printing "=== MSG ===\n" to the dumpfile etc,
   and then calling
     dump_end_scope ();
   once the object goes out of scope, thus capturing the nesting of
   the scopes.

   These scopes affect dump messages within them: dump messages at the
   top level implicitly default to MSG_PRIORITY_USER_FACING, whereas those
   in a nested scope implicitly default to MSG_PRIORITY_INTERNALS.  */

#define DUMP_VECT_SCOPE(MSG) \
  AUTO_DUMP_SCOPE (MSG, vect_location)

/* A sentinel class for ensuring that the "vect_location" global gets
   reset at the end of a scope.

   The "vect_location" global is used during dumping and contains a
   location_t, which could contain references to a tree block via the
   ad-hoc data.  This data is used for tracking inlining information,
   but it's not a GC root; it's simply assumed that such locations never
   get accessed if the blocks are optimized away.

   Hence we need to ensure that such locations are purged at the end
   of any operations using them (e.g. via this class).  */

class auto_purge_vect_location
{
 public:
  ~auto_purge_vect_location ();
};

/*-----------------------------------------------------------------*/
/* Function prototypes.                                            */
/*-----------------------------------------------------------------*/

/* Simple loop peeling and versioning utilities for vectorizer's purposes -
   in tree-vect-loop-manip.cc.  */
extern void vect_set_loop_condition (class loop *, loop_vec_info,
				     tree, tree, tree, bool);
extern bool slpeel_can_duplicate_loop_p (const class loop *, const_edge);
class loop *slpeel_tree_duplicate_loop_to_edge_cfg (class loop *,
						     class loop *, edge);
class loop *vect_loop_versioning (loop_vec_info, gimple *);
extern class loop *vect_do_peeling (loop_vec_info, tree, tree,
				    tree *, tree *, tree *, int, bool, bool,
				    tree *);
extern tree vect_get_main_loop_result (loop_vec_info, tree, tree);
extern void vect_prepare_for_masked_peels (loop_vec_info);
extern dump_user_location_t find_loop_location (class loop *);
extern bool vect_can_advance_ivs_p (loop_vec_info);
extern void vect_update_inits_of_drs (loop_vec_info, tree, tree_code);

/* In tree-vect-stmts.cc.  */
extern tree get_related_vectype_for_scalar_type (machine_mode, tree,
						 poly_uint64 = 0);
extern tree get_vectype_for_scalar_type (vec_info *, tree, unsigned int = 0);
extern tree get_vectype_for_scalar_type (vec_info *, tree, slp_tree);
extern tree get_mask_type_for_scalar_type (vec_info *, tree, unsigned int = 0);
extern tree get_same_sized_vectype (tree, tree);
extern bool vect_chooses_same_modes_p (vec_info *, machine_mode);
extern bool vect_get_loop_mask_type (loop_vec_info);
extern bool vect_is_simple_use (tree, vec_info *, enum vect_def_type *,
				stmt_vec_info * = NULL, gimple ** = NULL);
extern bool vect_is_simple_use (tree, vec_info *, enum vect_def_type *,
				tree *, stmt_vec_info * = NULL,
				gimple ** = NULL);
extern bool vect_is_simple_use (vec_info *, stmt_vec_info, slp_tree,
				unsigned, tree *, slp_tree *,
				enum vect_def_type *,
				tree *, stmt_vec_info * = NULL);
extern bool vect_maybe_update_slp_op_vectype (slp_tree, tree);
extern bool supportable_widening_operation (vec_info *,
					    enum tree_code, stmt_vec_info,
					    tree, tree, enum tree_code *,
					    enum tree_code *, int *,
					    vec<tree> *);
extern bool supportable_narrowing_operation (enum tree_code, tree, tree,
					     enum tree_code *, int *,
					     vec<tree> *);

extern unsigned record_stmt_cost (stmt_vector_for_cost *, int,
				  enum vect_cost_for_stmt, stmt_vec_info,
				  tree, int, enum vect_cost_model_location);
extern unsigned record_stmt_cost (stmt_vector_for_cost *, int,
				  enum vect_cost_for_stmt, slp_tree,
				  tree, int, enum vect_cost_model_location);
extern unsigned record_stmt_cost (stmt_vector_for_cost *, int,
				  enum vect_cost_for_stmt,
				  enum vect_cost_model_location);

/* Overload of record_stmt_cost with VECTYPE derived from STMT_INFO.  */

static inline unsigned
record_stmt_cost (stmt_vector_for_cost *body_cost_vec, int count,
		  enum vect_cost_for_stmt kind, stmt_vec_info stmt_info,
		  int misalign, enum vect_cost_model_location where)
{
  return record_stmt_cost (body_cost_vec, count, kind, stmt_info,
			   STMT_VINFO_VECTYPE (stmt_info), misalign, where);
}

extern void vect_finish_replace_stmt (vec_info *, stmt_vec_info, gimple *);
extern void vect_finish_stmt_generation (vec_info *, stmt_vec_info, gimple *,
					 gimple_stmt_iterator *);
extern opt_result vect_mark_stmts_to_be_vectorized (loop_vec_info, bool *);
extern tree vect_get_store_rhs (stmt_vec_info);
void vect_get_vec_defs_for_operand (vec_info *vinfo, stmt_vec_info, unsigned,
				    tree op, vec<tree> *, tree = NULL);
void vect_get_vec_defs (vec_info *, stmt_vec_info, slp_tree, unsigned,
			tree, vec<tree> *,
			tree = NULL, vec<tree> * = NULL,
			tree = NULL, vec<tree> * = NULL,
			tree = NULL, vec<tree> * = NULL);
void vect_get_vec_defs (vec_info *, stmt_vec_info, slp_tree, unsigned,
			tree, vec<tree> *, tree,
			tree = NULL, vec<tree> * = NULL, tree = NULL,
			tree = NULL, vec<tree> * = NULL, tree = NULL,
			tree = NULL, vec<tree> * = NULL, tree = NULL);
extern tree vect_init_vector (vec_info *, stmt_vec_info, tree, tree,
                              gimple_stmt_iterator *);
extern tree vect_get_slp_vect_def (slp_tree, unsigned);
extern bool vect_transform_stmt (vec_info *, stmt_vec_info,
				 gimple_stmt_iterator *,
				 slp_tree, slp_instance);
extern void vect_remove_stores (vec_info *, stmt_vec_info);
extern bool vect_nop_conversion_p (stmt_vec_info);
extern opt_result vect_analyze_stmt (vec_info *, stmt_vec_info, bool *,
				     slp_tree,
				     slp_instance, stmt_vector_for_cost *);
extern void vect_get_load_cost (vec_info *, stmt_vec_info, int,
				dr_alignment_support, int, bool,
				unsigned int *, unsigned int *,
				stmt_vector_for_cost *,
				stmt_vector_for_cost *, bool);
extern void vect_get_store_cost (vec_info *, stmt_vec_info, int,
				 dr_alignment_support, int,
				 unsigned int *, stmt_vector_for_cost *);
extern bool vect_supportable_shift (vec_info *, enum tree_code, tree);
extern tree vect_gen_perm_mask_any (tree, const vec_perm_indices &);
extern tree vect_gen_perm_mask_checked (tree, const vec_perm_indices &);
extern void optimize_mask_stores (class loop*);
extern tree vect_gen_while (gimple_seq *, tree, tree, tree,
			    const char * = nullptr);
extern tree vect_gen_while_not (gimple_seq *, tree, tree, tree);
extern opt_result vect_get_vector_types_for_stmt (vec_info *,
						  stmt_vec_info, tree *,
						  tree *, unsigned int = 0);
extern opt_tree vect_get_mask_type_for_stmt (stmt_vec_info, unsigned int = 0);

/* In tree-vect-data-refs.cc.  */
extern bool vect_can_force_dr_alignment_p (const_tree, poly_uint64);
extern enum dr_alignment_support vect_supportable_dr_alignment
				   (vec_info *, dr_vec_info *, tree, int);
extern tree vect_get_smallest_scalar_type (stmt_vec_info, tree);
extern opt_result vect_analyze_data_ref_dependences (loop_vec_info, unsigned int *);
extern bool vect_slp_analyze_instance_dependence (vec_info *, slp_instance);
extern opt_result vect_enhance_data_refs_alignment (loop_vec_info);
extern opt_result vect_analyze_data_refs_alignment (loop_vec_info);
extern bool vect_slp_analyze_instance_alignment (vec_info *, slp_instance);
extern opt_result vect_analyze_data_ref_accesses (vec_info *, vec<int> *);
extern opt_result vect_prune_runtime_alias_test_list (loop_vec_info);
extern bool vect_gather_scatter_fn_p (vec_info *, bool, bool, tree, tree,
				      tree, int, internal_fn *, tree *);
extern bool vect_check_gather_scatter (stmt_vec_info, loop_vec_info,
				       gather_scatter_info *);
extern opt_result vect_find_stmt_data_reference (loop_p, gimple *,
						 vec<data_reference_p> *,
						 vec<int> *, int);
extern opt_result vect_analyze_data_refs (vec_info *, poly_uint64 *, bool *);
extern void vect_record_base_alignments (vec_info *);
extern tree vect_create_data_ref_ptr (vec_info *,
				      stmt_vec_info, tree, class loop *, tree,
				      tree *, gimple_stmt_iterator *,
				      gimple **, bool,
				      tree = NULL_TREE);
extern tree bump_vector_ptr (vec_info *, tree, gimple *, gimple_stmt_iterator *,
			     stmt_vec_info, tree);
extern void vect_copy_ref_info (tree, tree);
extern tree vect_create_destination_var (tree, tree);
extern bool vect_grouped_store_supported (tree, unsigned HOST_WIDE_INT);
extern bool vect_store_lanes_supported (tree, unsigned HOST_WIDE_INT, bool);
extern bool vect_grouped_load_supported (tree, bool, unsigned HOST_WIDE_INT);
extern bool vect_load_lanes_supported (tree, unsigned HOST_WIDE_INT, bool);
extern void vect_permute_store_chain (vec_info *, vec<tree> &,
				      unsigned int, stmt_vec_info,
				      gimple_stmt_iterator *, vec<tree> *);
extern tree vect_setup_realignment (vec_info *,
				    stmt_vec_info, gimple_stmt_iterator *,
				    tree *, enum dr_alignment_support, tree,
	                            class loop **);
extern void vect_transform_grouped_load (vec_info *, stmt_vec_info, vec<tree>,
					 int, gimple_stmt_iterator *);
extern void vect_record_grouped_load_vectors (vec_info *,
					      stmt_vec_info, vec<tree>);
extern tree vect_get_new_vect_var (tree, enum vect_var_kind, const char *);
extern tree vect_get_new_ssa_name (tree, enum vect_var_kind,
				   const char * = NULL);
extern tree vect_create_addr_base_for_vector_ref (vec_info *,
						  stmt_vec_info, gimple_seq *,
						  tree);

/* In tree-vect-loop.cc.  */
extern tree neutral_op_for_reduction (tree, code_helper, tree);
extern widest_int vect_iv_limit_for_partial_vectors (loop_vec_info loop_vinfo);
bool vect_rgroup_iv_might_wrap_p (loop_vec_info, rgroup_controls *);
/* Used in tree-vect-loop-manip.cc */
extern opt_result vect_determine_partial_vectors_and_peeling (loop_vec_info,
							      bool);
/* Used in gimple-loop-interchange.c and tree-parloops.cc.  */
extern bool check_reduction_path (dump_user_location_t, loop_p, gphi *, tree,
				  enum tree_code);
extern bool needs_fold_left_reduction_p (tree, code_helper);
/* Drive for loop analysis stage.  */
extern opt_loop_vec_info vect_analyze_loop (class loop *, vec_info_shared *);
extern tree vect_build_loop_niters (loop_vec_info, bool * = NULL);
extern void vect_gen_vector_loop_niters (loop_vec_info, tree, tree *,
					 tree *, bool);
extern tree vect_halve_mask_nunits (tree, machine_mode);
extern tree vect_double_mask_nunits (tree, machine_mode);
extern void vect_record_loop_mask (loop_vec_info, vec_loop_masks *,
				   unsigned int, tree, tree);
extern tree vect_get_loop_mask (gimple_stmt_iterator *, vec_loop_masks *,
				unsigned int, tree, unsigned int);
extern void vect_record_loop_len (loop_vec_info, vec_loop_lens *, unsigned int,
				  tree, unsigned int);
extern tree vect_get_loop_len (loop_vec_info, vec_loop_lens *, unsigned int,
			       unsigned int);
extern gimple_seq vect_gen_len (tree, tree, tree, tree);
extern stmt_vec_info info_for_reduction (vec_info *, stmt_vec_info);
extern bool reduction_fn_for_scalar_code (code_helper, internal_fn *);

/* Drive for loop transformation stage.  */
extern class loop *vect_transform_loop (loop_vec_info, gimple *);
struct vect_loop_form_info
{
  tree number_of_iterations;
  tree number_of_iterationsm1;
  tree assumptions;
  gcond *loop_cond;
  gcond *inner_loop_cond;
};
extern opt_result vect_analyze_loop_form (class loop *, vect_loop_form_info *);
extern loop_vec_info vect_create_loop_vinfo (class loop *, vec_info_shared *,
					     const vect_loop_form_info *,
					     loop_vec_info = nullptr);
extern bool vectorizable_live_operation (vec_info *,
					 stmt_vec_info, gimple_stmt_iterator *,
					 slp_tree, slp_instance, int,
					 bool, stmt_vector_for_cost *);
extern bool vectorizable_reduction (loop_vec_info, stmt_vec_info,
				    slp_tree, slp_instance,
				    stmt_vector_for_cost *);
extern bool vectorizable_induction (loop_vec_info, stmt_vec_info,
				    gimple **, slp_tree,
				    stmt_vector_for_cost *);
extern bool vect_transform_reduction (loop_vec_info, stmt_vec_info,
				      gimple_stmt_iterator *,
				      gimple **, slp_tree);
extern bool vect_transform_cycle_phi (loop_vec_info, stmt_vec_info,
				      gimple **,
				      slp_tree, slp_instance);
extern bool vectorizable_lc_phi (loop_vec_info, stmt_vec_info,
				 gimple **, slp_tree);
extern bool vectorizable_phi (vec_info *, stmt_vec_info, gimple **, slp_tree,
			      stmt_vector_for_cost *);
extern bool vect_emulated_vector_p (tree);
extern bool vect_can_vectorize_without_simd_p (tree_code);
extern bool vect_can_vectorize_without_simd_p (code_helper);
extern int vect_get_known_peeling_cost (loop_vec_info, int, int *,
					stmt_vector_for_cost *,
					stmt_vector_for_cost *,
					stmt_vector_for_cost *);
extern tree cse_and_gimplify_to_preheader (loop_vec_info, tree);

/* In tree-vect-slp.cc.  */
extern void vect_slp_init (void);
extern void vect_slp_fini (void);
extern void vect_free_slp_instance (slp_instance);
extern bool vect_transform_slp_perm_load (vec_info *, slp_tree, const vec<tree> &,
					  gimple_stmt_iterator *, poly_uint64,
					  bool, unsigned *,
					  unsigned * = nullptr, bool = false);
extern bool vect_slp_analyze_operations (vec_info *);
extern void vect_schedule_slp (vec_info *, const vec<slp_instance> &);
extern opt_result vect_analyze_slp (vec_info *, unsigned);
extern bool vect_make_slp_decision (loop_vec_info);
extern void vect_detect_hybrid_slp (loop_vec_info);
extern void vect_optimize_slp (vec_info *);
extern void vect_gather_slp_loads (vec_info *);
extern void vect_get_slp_defs (slp_tree, vec<tree> *);
extern void vect_get_slp_defs (vec_info *, slp_tree, vec<vec<tree> > *,
			       unsigned n = -1U);
extern bool vect_slp_if_converted_bb (basic_block bb, loop_p orig_loop);
extern bool vect_slp_function (function *);
extern stmt_vec_info vect_find_last_scalar_stmt_in_slp (slp_tree);
extern stmt_vec_info vect_find_first_scalar_stmt_in_slp (slp_tree);
extern bool is_simple_and_all_uses_invariant (stmt_vec_info, loop_vec_info);
extern bool can_duplicate_and_interleave_p (vec_info *, unsigned int, tree,
					    unsigned int * = NULL,
					    tree * = NULL, tree * = NULL);
extern void duplicate_and_interleave (vec_info *, gimple_seq *, tree,
				      const vec<tree> &, unsigned int, vec<tree> &);
extern int vect_get_place_in_interleaving_chain (stmt_vec_info, stmt_vec_info);
extern slp_tree vect_create_new_slp_node (unsigned, tree_code);
extern void vect_free_slp_tree (slp_tree);
extern bool compatible_calls_p (gcall *, gcall *);

/* In tree-vect-patterns.cc.  */
extern void
vect_mark_pattern_stmts (vec_info *, stmt_vec_info, gimple *, tree);

/* Pattern recognition functions.
   Additional pattern recognition functions can (and will) be added
   in the future.  */
void vect_pattern_recog (vec_info *);

/* In tree-vectorizer.cc.  */
unsigned vectorize_loops (void);
void vect_free_loop_info_assumptions (class loop *);
gimple *vect_loop_vectorized_call (class loop *, gcond **cond = NULL);
bool vect_stmt_dominates_stmt_p (gimple *, gimple *);

/* SLP Pattern matcher types, tree-vect-slp-patterns.cc.  */

/* Forward declaration of possible two operands operation that can be matched
   by the complex numbers pattern matchers.  */
enum _complex_operation : unsigned;

/* All possible load permute values that could result from the partial data-flow
   analysis.  */
typedef enum _complex_perm_kinds {
   PERM_UNKNOWN,
   PERM_EVENODD,
   PERM_ODDEVEN,
   PERM_ODDODD,
   PERM_EVENEVEN,
   /* Can be combined with any other PERM values.  */
   PERM_TOP
} complex_perm_kinds_t;

/* Cache from nodes to the load permutation they represent.  */
typedef hash_map <slp_tree, complex_perm_kinds_t>
  slp_tree_to_load_perm_map_t;

/* Cache from nodes pair to being compatible or not.  */
typedef pair_hash <nofree_ptr_hash <_slp_tree>,
		   nofree_ptr_hash <_slp_tree>> slp_node_hash;
typedef hash_map <slp_node_hash, bool> slp_compat_nodes_map_t;


/* Vector pattern matcher base class.  All SLP pattern matchers must inherit
   from this type.  */

class vect_pattern
{
  protected:
    /* The number of arguments that the IFN requires.  */
    unsigned m_num_args;

    /* The internal function that will be used when a pattern is created.  */
    internal_fn m_ifn;

    /* The current node being inspected.  */
    slp_tree *m_node;

    /* The list of operands to be the children for the node produced when the
       internal function is created.  */
    vec<slp_tree> m_ops;

    /* Default constructor where NODE is the root of the tree to inspect.  */
    vect_pattern (slp_tree *node, vec<slp_tree> *m_ops, internal_fn ifn)
    {
      this->m_ifn = ifn;
      this->m_node = node;
      this->m_ops.create (0);
      if (m_ops)
	this->m_ops.safe_splice (*m_ops);
    }

  public:

    /* Create a new instance of the pattern matcher class of the given type.  */
    static vect_pattern* recognize (slp_tree_to_load_perm_map_t *,
				    slp_compat_nodes_map_t *, slp_tree *);

    /* Build the pattern from the data collected so far.  */
    virtual void build (vec_info *) = 0;

    /* Default destructor.  */
    virtual ~vect_pattern ()
    {
	this->m_ops.release ();
    }
};

/* Function pointer to create a new pattern matcher from a generic type.  */
typedef vect_pattern* (*vect_pattern_decl_t) (slp_tree_to_load_perm_map_t *,
					      slp_compat_nodes_map_t *,
					      slp_tree *);

/* List of supported pattern matchers.  */
extern vect_pattern_decl_t slp_patterns[];

/* Number of supported pattern matchers.  */
extern size_t num__slp_patterns;

/* ----------------------------------------------------------------------
   Target support routines
   -----------------------------------------------------------------------
   The following routines are provided to simplify costing decisions in
   target code.  Please add more as needed.  */

/* Return true if an operaton of kind KIND for STMT_INFO represents
   the extraction of an element from a vector in preparation for
   storing the element to memory.  */
inline bool
vect_is_store_elt_extraction (vect_cost_for_stmt kind, stmt_vec_info stmt_info)
{
  return (kind == vec_to_scalar
	  && STMT_VINFO_DATA_REF (stmt_info)
	  && DR_IS_WRITE (STMT_VINFO_DATA_REF (stmt_info)));
}

/* Return true if STMT_INFO represents part of a reduction.  */
inline bool
vect_is_reduction (stmt_vec_info stmt_info)
{
  return STMT_VINFO_REDUC_IDX (stmt_info) >= 0;
}

/* If STMT_INFO describes a reduction, return the vect_reduction_type
   of the reduction it describes, otherwise return -1.  */
inline int
vect_reduc_type (vec_info *vinfo, stmt_vec_info stmt_info)
{
  if (loop_vec_info loop_vinfo = dyn_cast<loop_vec_info> (vinfo))
    if (STMT_VINFO_REDUC_DEF (stmt_info))
      {
	stmt_vec_info reduc_info = info_for_reduction (loop_vinfo, stmt_info);
	return int (STMT_VINFO_REDUC_TYPE (reduc_info));
      }
  return -1;
}

/* If STMT_INFO is a COND_EXPR that includes an embedded comparison, return the
   scalar type of the values being compared.  Return null otherwise.  */
inline tree
vect_embedded_comparison_type (stmt_vec_info stmt_info)
{
  if (auto *assign = dyn_cast<gassign *> (stmt_info->stmt))
    if (gimple_assign_rhs_code (assign) == COND_EXPR)
      {
	tree cond = gimple_assign_rhs1 (assign);
	if (COMPARISON_CLASS_P (cond))
	  return TREE_TYPE (TREE_OPERAND (cond, 0));
      }
  return NULL_TREE;
}

/* If STMT_INFO is a comparison or contains an embedded comparison, return the
   scalar type of the values being compared.  Return null otherwise.  */
inline tree
vect_comparison_type (stmt_vec_info stmt_info)
{
  if (auto *assign = dyn_cast<gassign *> (stmt_info->stmt))
    if (TREE_CODE_CLASS (gimple_assign_rhs_code (assign)) == tcc_comparison)
      return TREE_TYPE (gimple_assign_rhs1 (assign));
  return vect_embedded_comparison_type (stmt_info);
}

/* Return true if STMT_INFO extends the result of a load.  */
inline bool
vect_is_extending_load (class vec_info *vinfo, stmt_vec_info stmt_info)
{
  /* Although this is quite large for an inline function, this part
     at least should be inline.  */
  gassign *assign = dyn_cast <gassign *> (stmt_info->stmt);
  if (!assign || !CONVERT_EXPR_CODE_P (gimple_assign_rhs_code (assign)))
    return false;

  tree rhs = gimple_assign_rhs1 (stmt_info->stmt);
  tree lhs_type = TREE_TYPE (gimple_assign_lhs (assign));
  tree rhs_type = TREE_TYPE (rhs);
  if (!INTEGRAL_TYPE_P (lhs_type)
      || !INTEGRAL_TYPE_P (rhs_type)
      || TYPE_PRECISION (lhs_type) <= TYPE_PRECISION (rhs_type))
    return false;

  stmt_vec_info def_stmt_info = vinfo->lookup_def (rhs);
  return (def_stmt_info
	  && STMT_VINFO_DATA_REF (def_stmt_info)
	  && DR_IS_READ (STMT_VINFO_DATA_REF (def_stmt_info)));
}

/* Return true if STMT_INFO is an integer truncation.  */
inline bool
vect_is_integer_truncation (stmt_vec_info stmt_info)
{
  gassign *assign = dyn_cast <gassign *> (stmt_info->stmt);
  if (!assign || !CONVERT_EXPR_CODE_P (gimple_assign_rhs_code (assign)))
    return false;

  tree lhs_type = TREE_TYPE (gimple_assign_lhs (assign));
  tree rhs_type = TREE_TYPE (gimple_assign_rhs1 (assign));
  return (INTEGRAL_TYPE_P (lhs_type)
	  && INTEGRAL_TYPE_P (rhs_type)
	  && TYPE_PRECISION (lhs_type) < TYPE_PRECISION (rhs_type));
}

#endif  /* GCC_TREE_VECTORIZER_H  */
