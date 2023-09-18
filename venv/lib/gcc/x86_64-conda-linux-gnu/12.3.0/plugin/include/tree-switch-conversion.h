/* Tree switch conversion for GNU compiler.
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
<http://www.gnu.org/licenses/>.  */

#ifndef TREE_SWITCH_CONVERSION_H
#define TREE_SWITCH_CONVERSION_H

namespace tree_switch_conversion {

/* Type of cluster.  */

enum cluster_type
{
  SIMPLE_CASE,
  JUMP_TABLE,
  BIT_TEST
};

#define PRINT_CASE(f,c) print_generic_expr (f, c)

/* Abstract base class for representing a cluster of cases.

   Here is the inheritance hierarachy, and the enum_cluster_type
   values for the concrete subclasses:

   cluster
   |-simple_cluster (SIMPLE_CASE)
   `-group_cluster
     |-jump_table_cluster (JUMP_TABLE)
     `-bit_test_cluster   (BIT_TEST).  */

class cluster
{
public:
  /* Constructor.  */
  inline cluster (tree case_label_expr, basic_block case_bb,
		  profile_probability prob, profile_probability subtree_prob);

  /* Destructor.  */
  virtual ~cluster ()
  {}

  /* Return type.  */
  virtual cluster_type get_type () = 0;

  /* Get low value covered by a cluster.  */
  virtual tree get_low () = 0;

  /* Get high value covered by a cluster.  */
  virtual tree get_high () = 0;

  /* Debug content of a cluster.  */
  virtual void debug () = 0;

  /* Dump content of a cluster.  */
  virtual void dump (FILE *f, bool details = false) = 0;

  /* Emit GIMPLE code to handle the cluster.  */
  virtual void emit (tree, tree, tree, basic_block, location_t) = 0;

  /* Return true if a cluster handles only a single case value and the
     value is not a range.  */
  virtual bool is_single_value_p ()
  {
    return false;
  }

  /* Return range of a cluster.  If value would overflow in type of LOW,
     then return 0.  */
  static unsigned HOST_WIDE_INT get_range (tree low, tree high)
  {
    wide_int w = wi::to_wide (high) - wi::to_wide (low);
    if (wi::neg_p (w, TYPE_SIGN (TREE_TYPE (low))) || !wi::fits_uhwi_p (w))
      return 0;
    return w.to_uhwi () + 1;
  }

  /* Case label.  */
  tree m_case_label_expr;

  /* Basic block of the case.  */
  basic_block m_case_bb;

  /* Probability of taking this cluster.  */
  profile_probability m_prob;

  /* Probability of reaching subtree rooted at this node.  */
  profile_probability m_subtree_prob;

protected:
  /* Default constructor.  */
  cluster () {}
};

cluster::cluster (tree case_label_expr, basic_block case_bb,
		  profile_probability prob, profile_probability subtree_prob):
  m_case_label_expr (case_label_expr), m_case_bb (case_bb), m_prob (prob),
  m_subtree_prob (subtree_prob)
{
}

/* Subclass of cluster representing a simple contiguous range
   from [low..high].  */

class simple_cluster: public cluster
{
public:
  /* Constructor.  */
  inline simple_cluster (tree low, tree high, tree case_label_expr,
			 basic_block case_bb, profile_probability prob,
			 bool has_forward_bb = false);

  /* Destructor.  */
  ~simple_cluster ()
  {}

  cluster_type
  get_type ()
  {
    return SIMPLE_CASE;
  }

  tree
  get_low ()
  {
    return m_low;
  }

  tree
  get_high ()
  {
    return m_high;
  }

  void set_high (tree high)
  {
    m_high = high;
  }

  void
  debug ()
  {
    dump (stderr);
  }

  void
  dump (FILE *f, bool details ATTRIBUTE_UNUSED = false)
  {
    PRINT_CASE (f, get_low ());
    if (get_low () != get_high ())
      {
	fprintf (f, "-");
	PRINT_CASE (f, get_high ());
      }
    fprintf (f, " ");
  }

  void emit (tree, tree, tree, basic_block, location_t)
  {
    gcc_unreachable ();
  }

  bool is_single_value_p ()
  {
    return tree_int_cst_equal (get_low (), get_high ());
  }

  /* Return number of comparisons needed for the case.  */
  unsigned
  get_comparison_count ()
  {
    return m_range_p ? 2 : 1;
  }

  /* Low value of the case.  */
  tree m_low;

  /* High value of the case.  */
  tree m_high;

  /* True if case is a range.  */
  bool m_range_p;

  /* True if the case will use a forwarder BB.  */
  bool m_has_forward_bb;
};

simple_cluster::simple_cluster (tree low, tree high, tree case_label_expr,
				basic_block case_bb, profile_probability prob,
				bool has_forward_bb):
  cluster (case_label_expr, case_bb, prob, prob),
  m_low (low), m_high (high), m_has_forward_bb (has_forward_bb)
{
  m_range_p = m_high != NULL;
  if (m_high == NULL)
    m_high = m_low;
}

/* Abstract subclass of jump table and bit test cluster,
   handling a collection of simple_cluster instances.  */

class group_cluster: public cluster
{
public:
  /* Constructor.  */
  group_cluster (vec<cluster *> &clusters, unsigned start, unsigned end);

  /* Destructor.  */
  ~group_cluster ();

  tree
  get_low ()
  {
    return m_cases[0]->get_low ();
  }

  tree
  get_high ()
  {
    return m_cases[m_cases.length () - 1]->get_high ();
  }

  void
  debug ()
  {
    dump (stderr);
  }

  void dump (FILE *f, bool details = false);

  /* List of simple clusters handled by the group.  */
  vec<simple_cluster *> m_cases;
};

/* Concrete subclass of group_cluster representing a collection
   of cases to be implemented as a jump table.
   The "emit" vfunc gernerates a nested switch statement which
   is later lowered to a jump table.  */

class jump_table_cluster: public group_cluster
{
public:
  /* Constructor.  */
  jump_table_cluster (vec<cluster *> &clusters, unsigned start, unsigned end)
  : group_cluster (clusters, start, end)
  {}

  cluster_type
  get_type ()
  {
    return JUMP_TABLE;
  }

  void emit (tree index_expr, tree index_type,
	     tree default_label_expr, basic_block default_bb, location_t loc);

  /* Find jump tables of given CLUSTERS, where all members of the vector
     are of type simple_cluster.  New clusters are returned.  */
  static vec<cluster *> find_jump_tables (vec<cluster *> &clusters);

  /* Return true when cluster starting at START and ending at END (inclusive)
     can build a jump-table.  COMPARISON_COUNT is number of comparison
     operations needed if the clusters are expanded as decision tree.
     MAX_RATIO tells about the maximum code growth (in percent).  */
  static bool can_be_handled (const vec<cluster *> &clusters, unsigned start,
			      unsigned end, unsigned HOST_WIDE_INT max_ratio,
			      unsigned HOST_WIDE_INT comparison_count);

  /* Return true if cluster starting at START and ending at END (inclusive)
     is profitable transformation.  */
  static bool is_beneficial (const vec<cluster *> &clusters, unsigned start,
			     unsigned end);

  /* Return the smallest number of different values for which it is best
     to use a jump-table instead of a tree of conditional branches.  */
  static inline unsigned int case_values_threshold (void);

  /* Return whether jump table expansion is allowed.  */
  static inline bool is_enabled (void);
};

/* A GIMPLE switch statement can be expanded to a short sequence of bit-wise
comparisons.  "switch(x)" is converted into "if ((1 << (x-MINVAL)) & CST)"
where CST and MINVAL are integer constants.  This is better than a series
of compare-and-banch insns in some cases,  e.g. we can implement:

	if ((x==4) || (x==6) || (x==9) || (x==11))

as a single bit test:

	if ((1<<x) & ((1<<4)|(1<<6)|(1<<9)|(1<<11)))

This transformation is only applied if the number of case targets is small,
if CST constains at least 3 bits, and "1 << x" is cheap.  The bit tests are
performed in "word_mode".

The following example shows the code the transformation generates:

	int bar(int x)
	{
		switch (x)
		{
		case '0':  case '1':  case '2':  case '3':  case '4':
		case '5':  case '6':  case '7':  case '8':  case '9':
		case 'A':  case 'B':  case 'C':  case 'D':  case 'E':
		case 'F':
			return 1;
		}
		return 0;
	}

==>

	bar (int x)
	{
		tmp1 = x - 48;
		if (tmp1 > (70 - 48)) goto L2;
		tmp2 = 1 << tmp1;
		tmp3 = 0b11111100000001111111111;
		if ((tmp2 & tmp3) != 0) goto L1 ; else goto L2;
	L1:
		return 1;
	L2:
		return 0;
	}

TODO: There are still some improvements to this transformation that could
be implemented:

* A narrower mode than word_mode could be used if that is cheaper, e.g.
  for x86_64 where a narrower-mode shift may result in smaller code.

* The compounded constant could be shifted rather than the one.  The
  test would be either on the sign bit or on the least significant bit,
  depending on the direction of the shift.  On some machines, the test
  for the branch would be free if the bit to test is already set by the
  shift operation.

This transformation was contributed by Roger Sayle, see this e-mail:
   http://gcc.gnu.org/ml/gcc-patches/2003-01/msg01950.html
*/

class bit_test_cluster: public group_cluster
{
public:
  /* Constructor.  */
  bit_test_cluster (vec<cluster *> &clusters, unsigned start, unsigned end,
		    bool handles_entire_switch)
  :group_cluster (clusters, start, end),
  m_handles_entire_switch (handles_entire_switch)
  {}

  cluster_type
  get_type ()
  {
    return BIT_TEST;
  }

/*  Expand a switch statement by a short sequence of bit-wise
    comparisons.  "switch(x)" is effectively converted into
    "if ((1 << (x-MINVAL)) & CST)" where CST and MINVAL are
    integer constants.

    INDEX_EXPR is the value being switched on.

    MINVAL is the lowest case value of in the case nodes,
    and RANGE is highest value minus MINVAL.  MINVAL and RANGE
    are not guaranteed to be of the same type as INDEX_EXPR
    (the gimplifier doesn't change the type of case label values,
    and MINVAL and RANGE are derived from those values).
    MAXVAL is MINVAL + RANGE.

    There *MUST* be max_case_bit_tests or less unique case
    node targets.  */
  void emit (tree index_expr, tree index_type,
	     tree default_label_expr, basic_block default_bb, location_t loc);

  /* Find bit tests of given CLUSTERS, where all members of the vector
     are of type simple_cluster.  New clusters are returned.  */
  static vec<cluster *> find_bit_tests (vec<cluster *> &clusters);

  /* Return true when RANGE of case values with UNIQ labels
     can build a bit test.  */
  static bool can_be_handled (unsigned HOST_WIDE_INT range, unsigned uniq);

  /* Return true when cluster starting at START and ending at END (inclusive)
     can build a bit test.  */
  static bool can_be_handled (const vec<cluster *> &clusters, unsigned start,
			      unsigned end);

  /* Return true when COUNT of cases of UNIQ labels is beneficial for bit test
     transformation.  */
  static bool is_beneficial (unsigned count, unsigned uniq);

  /* Return true if cluster starting at START and ending at END (inclusive)
     is profitable transformation.  */
  static bool is_beneficial (const vec<cluster *> &clusters, unsigned start,
			     unsigned end);

/* Split the basic block at the statement pointed to by GSIP, and insert
   a branch to the target basic block of E_TRUE conditional on tree
   expression COND.

   It is assumed that there is already an edge from the to-be-split
   basic block to E_TRUE->dest block.  This edge is removed, and the
   profile information on the edge is re-used for the new conditional
   jump.

   The CFG is updated.  The dominator tree will not be valid after
   this transformation, but the immediate dominators are updated if
   UPDATE_DOMINATORS is true.

   Returns the newly created basic block.  */
  static basic_block hoist_edge_and_branch_if_true (gimple_stmt_iterator *gsip,
						    tree cond,
						    basic_block case_bb,
						    profile_probability prob,
						    location_t);

  /* Return whether bit test expansion is allowed.  */
  static inline bool is_enabled (void)
  {
    return flag_bit_tests;
  }

  /* True when the jump table handles an entire switch statement.  */
  bool m_handles_entire_switch;

  /* Maximum number of different basic blocks that can be handled by
     a bit test.  */
  static const int m_max_case_bit_tests = 3;
};

/* Helper struct to find minimal clusters.  */

class min_cluster_item
{
public:
  /* Constructor.  */
  min_cluster_item (unsigned count, unsigned start, unsigned non_jt_cases):
    m_count (count), m_start (start), m_non_jt_cases (non_jt_cases)
  {}

  /* Count of clusters.  */
  unsigned m_count;

  /* Index where is cluster boundary.  */
  unsigned m_start;

  /* Total number of cases that will not be in a jump table.  */
  unsigned m_non_jt_cases;
};

/* Helper struct to represent switch decision tree.  */

class case_tree_node
{
public:
  /* Empty Constructor.  */
  case_tree_node ();

  /* Return true when it has a child.  */
  bool has_child ()
  {
    return m_left != NULL || m_right != NULL;
  }

  /* Left son in binary tree.  */
  case_tree_node *m_left;

  /* Right son in binary tree; also node chain.  */
  case_tree_node *m_right;

  /* Parent of node in binary tree.  */
  case_tree_node *m_parent;

  /* Cluster represented by this tree node.  */
  cluster *m_c;
};

inline
case_tree_node::case_tree_node ():
  m_left (NULL), m_right (NULL), m_parent (NULL), m_c (NULL)
{
}

unsigned int
jump_table_cluster::case_values_threshold (void)
{
  unsigned int threshold = param_case_values_threshold;

  if (threshold == 0)
    threshold = targetm.case_values_threshold ();

  return threshold;
}

/* Return whether jump table expansion is allowed.  */
bool jump_table_cluster::is_enabled (void)
{
  /* If neither casesi or tablejump is available, or flag_jump_tables
     over-ruled us, we really have no choice.  */
  if (!targetm.have_casesi () && !targetm.have_tablejump ())
    return false;
  if (!flag_jump_tables)
    return false;
#ifndef ASM_OUTPUT_ADDR_DIFF_ELT
  if (flag_pic)
    return false;
#endif

  return true;
}

/* A case_bit_test represents a set of case nodes that may be
   selected from using a bit-wise comparison.  HI and LO hold
   the integer to be tested against, TARGET_EDGE contains the
   edge to the basic block to jump to upon success and BITS
   counts the number of case nodes handled by this test,
   typically the number of bits set in HI:LO.  The LABEL field
   is used to quickly identify all cases in this set without
   looking at label_to_block for every case label.  */

class case_bit_test
{
public:
  wide_int mask;
  basic_block target_bb;
  tree label;
  int bits;

  /* Comparison function for qsort to order bit tests by decreasing
     probability of execution.  */
  static int cmp (const void *p1, const void *p2);
};

class switch_decision_tree
{
public:
  /* Constructor.  */
  switch_decision_tree (gswitch *swtch): m_switch (swtch), m_phi_mapping (),
    m_case_bbs (), m_case_node_pool ("struct case_node pool"),
    m_case_list (NULL)
  {
  }

  /* Analyze switch statement and return true when the statement is expanded
     as decision tree.  */
  bool analyze_switch_statement ();

  /* Attempt to expand CLUSTERS as a decision tree.  Return true when
     expanded.  */
  bool try_switch_expansion (vec<cluster *> &clusters);
  /* Compute the number of case labels that correspond to each outgoing edge of
     switch statement.  Record this information in the aux field of the edge.
     */
  void compute_cases_per_edge ();

  /* Before switch transformation, record all SSA_NAMEs defined in switch BB
     and used in a label basic block.  */
  void record_phi_operand_mapping ();

  /* Append new operands to PHI statements that were introduced due to
     addition of new edges to case labels.  */
  void fix_phi_operands_for_edges ();

  /* Generate a decision tree, switching on INDEX_EXPR and jumping to
     one of the labels in CASE_LIST or to the DEFAULT_LABEL.

     We generate a binary decision tree to select the appropriate target
     code.  */
  void emit (basic_block bb, tree index_expr,
	     profile_probability default_prob, tree index_type);

  /* Emit step-by-step code to select a case for the value of INDEX.
     The thus generated decision tree follows the form of the
     case-node binary tree NODE, whose nodes represent test conditions.
     DEFAULT_PROB is probability of cases leading to default BB.
     INDEX_TYPE is the type of the index of the switch.  */
  basic_block emit_case_nodes (basic_block bb, tree index,
			       case_tree_node *node,
			       profile_probability default_prob,
			       tree index_type, location_t);

  /* Take an ordered list of case nodes
     and transform them into a near optimal binary tree,
     on the assumption that any target code selection value is as
     likely as any other.

     The transformation is performed by splitting the ordered
     list into two equal sections plus a pivot.  The parts are
     then attached to the pivot as left and right branches.  Each
     branch is then transformed recursively.  */
  static void balance_case_nodes (case_tree_node **head,
				  case_tree_node *parent);

  /* Dump ROOT, a list or tree of case nodes, to file F.  */
  static void dump_case_nodes (FILE *f, case_tree_node *root, int indent_step,
			       int indent_level);

  /* Add an unconditional jump to CASE_BB that happens in basic block BB.  */
  static void emit_jump (basic_block bb, basic_block case_bb);

  /* Generate code to compare OP0 with OP1 so that the condition codes are
     set and to jump to LABEL_BB if the condition is true.
     COMPARISON is the GIMPLE comparison (EQ, NE, GT, etc.).
     PROB is the probability of jumping to LABEL_BB.  */
  static basic_block emit_cmp_and_jump_insns (basic_block bb, tree op0,
					      tree op1, tree_code comparison,
					      basic_block label_bb,
					      profile_probability prob,
					      location_t);

  /* Generate code to jump to LABEL if OP0 and OP1 are equal in mode MODE.
     PROB is the probability of jumping to LABEL_BB.  */
  static basic_block do_jump_if_equal (basic_block bb, tree op0, tree op1,
				       basic_block label_bb,
				       profile_probability prob,
				       location_t);

  /* Reset the aux field of all outgoing edges of switch basic block.  */
  static inline void reset_out_edges_aux (gswitch *swtch);

  /* Switch statement.  */
  gswitch *m_switch;

  /* Map of PHI nodes that have to be fixed after expansion.  */
  hash_map<tree, tree> m_phi_mapping;

  /* List of basic blocks that belong to labels of the switch.  */
  auto_vec<basic_block> m_case_bbs;

  /* Basic block with default label.  */
  basic_block m_default_bb;

  /* A pool for case nodes.  */
  object_allocator<case_tree_node> m_case_node_pool;

  /* Balanced tree of case nodes.  */
  case_tree_node *m_case_list;
};

/*
     Switch initialization conversion

The following pass changes simple initializations of scalars in a switch
statement into initializations from a static array.  Obviously, the values
must be constant and known at compile time and a default branch must be
provided.  For example, the following code:

	int a,b;

	switch (argc)
	{
	 case 1:
	 case 2:
		a_1 = 8;
		b_1 = 6;
		break;
	 case 3:
		a_2 = 9;
		b_2 = 5;
		break;
	 case 12:
		a_3 = 10;
		b_3 = 4;
		break;
	 default:
		a_4 = 16;
		b_4 = 1;
		break;
	}
	a_5 = PHI <a_1, a_2, a_3, a_4>
	b_5 = PHI <b_1, b_2, b_3, b_4>


is changed into:

	static const int = CSWTCH01[] = {6, 6, 5, 1, 1, 1, 1, 1, 1, 1, 1, 4};
	static const int = CSWTCH02[] = {8, 8, 9, 16, 16, 16, 16, 16, 16, 16,
				 16, 16, 10};

	if (((unsigned) argc) - 1 < 11)
	  {
	    a_6 = CSWTCH02[argc - 1];
	    b_6 = CSWTCH01[argc - 1];
	  }
	else
	  {
	    a_7 = 16;
	    b_7 = 1;
	  }
	a_5 = PHI <a_6, a_7>
	b_b = PHI <b_6, b_7>

There are further constraints.  Specifically, the range of values across all
case labels must not be bigger than param_switch_conversion_branch_ratio
(default eight) times the number of the actual switch branches.

This transformation was contributed by Martin Jambor, see this e-mail:
   http://gcc.gnu.org/ml/gcc-patches/2008-07/msg00011.html  */

/* The main structure of the pass.  */
class switch_conversion
{
public:
  /* Constructor.  */
  switch_conversion ();

  /* Destructor.  */
  ~switch_conversion ();

  /* The following function is invoked on every switch statement (the current
     one is given in SWTCH) and runs the individual phases of switch
     conversion on it one after another until one fails or the conversion
     is completed.  On success, NULL is in m_reason, otherwise points
     to a string with the reason why the conversion failed.  */
  void expand (gswitch *swtch);

  /* Collection information about SWTCH statement.  */
  void collect (gswitch *swtch);

  /* Checks whether the range given by individual case statements of the switch
     switch statement isn't too big and whether the number of branches actually
     satisfies the size of the new array.  */
  bool check_range ();

  /* Checks whether all but the final BB basic blocks are empty.  */
  bool check_all_empty_except_final ();

  /* This function checks whether all required values in phi nodes in final_bb
     are constants.  Required values are those that correspond to a basic block
     which is a part of the examined switch statement.  It returns true if the
     phi nodes are OK, otherwise false.  */
  bool check_final_bb ();

  /* The following function allocates default_values, target_{in,out}_names and
     constructors arrays.  The last one is also populated with pointers to
     vectors that will become constructors of new arrays.  */
  void create_temp_arrays ();

  /* Populate the array of default values in the order of phi nodes.
     DEFAULT_CASE is the CASE_LABEL_EXPR for the default switch branch
     if the range is non-contiguous or the default case has standard
     structure, otherwise it is the first non-default case instead.  */
  void gather_default_values (tree default_case);

  /* The following function populates the vectors in the constructors array with
     future contents of the static arrays.  The vectors are populated in the
     order of phi nodes.  */
  void build_constructors ();

  /* If all values in the constructor vector are products of a linear function
     a * x + b, then return true.  When true, COEFF_A and COEFF_B and
     coefficients of the linear function.  Note that equal values are special
     case of a linear function with a and b equal to zero.  */
  bool contains_linear_function_p (vec<constructor_elt, va_gc> *vec,
				   wide_int *coeff_a, wide_int *coeff_b);

  /* Return type which should be used for array elements, either TYPE's
     main variant or, for integral types, some smaller integral type
     that can still hold all the constants.  */
  tree array_value_type (tree type, int num);

  /* Create an appropriate array type and declaration and assemble a static
     array variable.  Also create a load statement that initializes
     the variable in question with a value from the static array.  SWTCH is
     the switch statement being converted, NUM is the index to
     arrays of constructors, default values and target SSA names
     for this particular array.  ARR_INDEX_TYPE is the type of the index
     of the new array, PHI is the phi node of the final BB that corresponds
     to the value that will be loaded from the created array.  TIDX
     is an ssa name of a temporary variable holding the index for loads from the
     new array.  */
  void build_one_array (int num, tree arr_index_type,
			gphi *phi, tree tidx);

  /* Builds and initializes static arrays initialized with values gathered from
     the switch statement.  Also creates statements that load values from
     them.  */
  void build_arrays ();

  /* Generates and appropriately inserts loads of default values at the position
     given by GSI.  Returns the last inserted statement.  */
  gassign *gen_def_assigns (gimple_stmt_iterator *gsi);

  /* Deletes the unused bbs and edges that now contain the switch statement and
     its empty branch bbs.  BBD is the now dead BB containing
     the original switch statement, FINAL is the last BB of the converted
     switch statement (in terms of succession).  */
  void prune_bbs (basic_block bbd, basic_block final, basic_block default_bb);

  /* Add values to phi nodes in final_bb for the two new edges.  E1F is the edge
     from the basic block loading values from an array and E2F from the basic
     block loading default values.  BBF is the last switch basic block (see the
     bbf description in the comment below).  */
  void fix_phi_nodes (edge e1f, edge e2f, basic_block bbf);

  /* Creates a check whether the switch expression value actually falls into the
     range given by all the cases.  If it does not, the temporaries are loaded
     with default values instead.  */
  void gen_inbound_check ();

  /* Switch statement for which switch conversion takes place.  */
  gswitch *m_switch;

  /* The expression used to decide the switch branch.  */
  tree m_index_expr;

  /* The following integer constants store the minimum and maximum value
     covered by the case labels.  */
  tree m_range_min;
  tree m_range_max;

  /* The difference between the above two numbers.  Stored here because it
     is used in all the conversion heuristics, as well as for some of the
     transformation, and it is expensive to re-compute it all the time.  */
  tree m_range_size;

  /* Basic block that contains the actual GIMPLE_SWITCH.  */
  basic_block m_switch_bb;

  /* Basic block that is the target of the default case.  */
  basic_block m_default_bb;

  /* The single successor block of all branches out of the GIMPLE_SWITCH,
     if such a block exists.  Otherwise NULL.  */
  basic_block m_final_bb;

  /* The probability of the default edge in the replaced switch.  */
  profile_probability m_default_prob;

  /* Number of phi nodes in the final bb (that we'll be replacing).  */
  int m_phi_count;

  /* Constructors of new static arrays.  */
  vec<constructor_elt, va_gc> **m_constructors;

  /* Array of default values, in the same order as phi nodes.  */
  tree *m_default_values;

  /* Array of ssa names that are initialized with a value from a new static
     array.  */
  tree *m_target_inbound_names;

  /* Array of ssa names that are initialized with the default value if the
     switch expression is out of range.  */
  tree *m_target_outbound_names;

  /* VOP SSA_NAME.  */
  tree m_target_vop;

  /* The first load statement that loads a temporary from a new static array.
   */
  gimple *m_arr_ref_first;

  /* The last load statement that loads a temporary from a new static array.  */
  gimple *m_arr_ref_last;

  /* String reason why the case wasn't a good candidate that is written to the
     dump file, if there is one.  */
  const char *m_reason;

  /* True if default case is not used for any value between range_min and
     range_max inclusive.  */
  bool m_contiguous_range;

  /* True if default case does not have the required shape for other case
     labels.  */
  bool m_default_case_nonstandard;

  /* Number of uniq labels for non-default edges.  */
  unsigned int m_uniq;

  /* Count is number of non-default edges.  */
  unsigned int m_count;

  /* True if CFG has been changed.  */
  bool m_cfg_altered;
};

void
switch_decision_tree::reset_out_edges_aux (gswitch *swtch)
{
  basic_block bb = gimple_bb (swtch);
  edge e;
  edge_iterator ei;
  FOR_EACH_EDGE (e, ei, bb->succs)
    e->aux = (void *) 0;
}

/* Release CLUSTERS vector and destruct all dynamically allocated items.  */

static inline void
release_clusters (vec<cluster *> &clusters)
{
  for (unsigned i = 0; i < clusters.length (); i++)
    delete clusters[i];
  clusters.release ();
}

} // tree_switch_conversion namespace

#endif // TREE_SWITCH_CONVERSION_H
