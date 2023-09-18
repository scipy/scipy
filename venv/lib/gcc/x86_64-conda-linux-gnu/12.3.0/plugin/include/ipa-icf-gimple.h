/* Interprocedural semantic function equality pass
   Copyright (C) 2014-2022 Free Software Foundation, Inc.

   Contributed by Jan Hubicka <hubicka@ucw.cz> and Martin Liska <mliska@suse.cz>

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

/* Gimple identical code folding (class func_checker) is an infrastructure
   capable of comparing two given functions. The class compares every
   gimple statement and uses many dictionaries to map source and target
   SSA_NAMEs, declarations and other components.

   To use the infrastructure, create an instance of func_checker and call
   a comparison function based on type of gimple statement.  */

/* Prints string STRING to a FILE with a given number of SPACE_COUNT.  */
#define FPUTS_SPACES(file, space_count, string) \
  fprintf (file, "%*s" string, space_count, " ");

/* fprintf function wrapper that transforms given FORMAT to follow given
   number for SPACE_COUNT and call fprintf for a FILE.  */
#define FPRINTF_SPACES(file, space_count, format, ...) \
  fprintf (file, "%*s" format, space_count, " ", ##__VA_ARGS__);

/* Logs a MESSAGE to dump_file if exists and returns false. FUNC is name
   of function and LINE is location in the source file.  */

static inline bool
return_false_with_message_1 (const char *message, const char *filename,
			     const char *func, unsigned int line)
{
  if (dump_file && (dump_flags & TDF_DETAILS))
    fprintf (dump_file, "  false returned: '%s' in %s at %s:%u\n", message, func,
	     filename, line);
  return false;
}

/* Logs a MESSAGE to dump_file if exists and returns false.  */
#define return_false_with_msg(message) \
  return_false_with_message_1 (message, __FILE__, __func__, __LINE__)

/* Return false and log that false value is returned.  */
#define return_false() return_false_with_msg ("")

/* Logs return value if RESULT is false. FUNC is name of function and LINE
   is location in the source file.  */

static inline bool
return_with_result (bool result, const char *filename,
		    const char *func, unsigned int line)
{
  if (!result && dump_file && (dump_flags & TDF_DETAILS))
    fprintf (dump_file, "  false returned: '' in %s at %s:%u\n", func,
	     filename, line);

  return result;
}

/* Logs return value if RESULT is false.  */
#define return_with_debug(result) return_with_result \
  (result, __FILE__, __func__, __LINE__)

/* Verbose logging function logging statements S1 and S2 of a CODE.
   FUNC is name of function and LINE is location in the source file.  */

static inline bool
return_different_stmts_1 (gimple *s1, gimple *s2, const char *code,
			  const char *func, unsigned int line)
{
  if (dump_file && (dump_flags & TDF_DETAILS))
    {
      fprintf (dump_file, "  different statement for code: %s (%s:%u):\n",
	       code, func, line);

      print_gimple_stmt (dump_file, s1, 3, TDF_DETAILS);
      print_gimple_stmt (dump_file, s2, 3, TDF_DETAILS);
    }

  return false;
}

/* Verbose logging function logging statements S1 and S2 of a CODE.  */
#define return_different_stmts(s1, s2, code) \
  return_different_stmts_1 (s1, s2, code, __func__, __LINE__)

namespace ipa_icf_gimple {

/* Basic block struct for semantic equality pass.  */
class sem_bb
{
public:
  sem_bb (basic_block bb_, unsigned nondbg_stmt_count_, unsigned edge_count_):
    bb (bb_), nondbg_stmt_count (nondbg_stmt_count_), edge_count (edge_count_) {}

  /* Basic block the structure belongs to.  */
  basic_block bb;

  /* Number of non-debug statements in the basic block.  */
  unsigned nondbg_stmt_count;

  /* Number of edges connected to the block.  */
  unsigned edge_count;
};

/* A class aggregating all connections and semantic equivalents
   for a given pair of semantic function candidates.  */
class func_checker : ao_compare
{
public:
  /* Default constructor.  */
  func_checker ():
    m_source_func_decl (NULL_TREE), m_target_func_decl (NULL_TREE),
    m_ignored_source_nodes (NULL), m_ignored_target_nodes (NULL),
    m_ignore_labels (false), m_tbaa (true)
  {
    m_source_ssa_names.create (0);
    m_target_ssa_names.create (0);
  }

  /* Initialize internal structures for a given SOURCE_FUNC_DECL and
     TARGET_FUNC_DECL. Strict polymorphic comparison is processed if
     an option COMPARE_POLYMORPHIC is true. For special cases, one can
     set IGNORE_LABELS to skip label comparison.
     Similarly, IGNORE_SOURCE_DECLS and IGNORE_TARGET_DECLS are sets
     of declarations that can be skipped.  */
  func_checker (tree source_func_decl, tree target_func_decl,
		bool ignore_labels = false,
		bool tbaa = true,
		hash_set<symtab_node *> *ignored_source_nodes = NULL,
		hash_set<symtab_node *> *ignored_target_nodes = NULL);

  /* Memory release routine.  */
  virtual ~func_checker ();

  /* Function visits all gimple labels and creates corresponding
     mapping between basic blocks and labels.  */
  void parse_labels (sem_bb *bb);

  /* Basic block equivalence comparison function that returns true if
     basic blocks BB1 and BB2 correspond.  */
  bool compare_bb (sem_bb *bb1, sem_bb *bb2);

  /* Verifies that trees T1 and T2 are equivalent from perspective of ICF.  */
  bool compare_ssa_name (const_tree t1, const_tree t2);

  /* Verification function for edges E1 and E2.  */
  bool compare_edge (edge e1, edge e2);

  /* Verifies for given GIMPLEs S1 and S2 that
     call statements are semantically equivalent.  */
  bool compare_gimple_call (gcall *s1, gcall *s2);

  /* Verifies for given GIMPLEs S1 and S2 that
     assignment statements are semantically equivalent.  */
  bool compare_gimple_assign (gimple *s1, gimple *s2);

  /* Verifies for given GIMPLEs S1 and S2 that
     condition statements are semantically equivalent.  */
  bool compare_gimple_cond (gimple *s1, gimple *s2);

  /* Verifies for given GIMPLE_LABEL stmts S1 and S2 that
     label statements are semantically equivalent.  */
  bool compare_gimple_label (const glabel *s1, const glabel *s2);

  /* Verifies for given GIMPLE_SWITCH stmts S1 and S2 that
     switch statements are semantically equivalent.  */
  bool compare_gimple_switch (const gswitch *s1, const gswitch *s2);

  /* Verifies for given GIMPLE_RETURN stmts S1 and S2 that
     return statements are semantically equivalent.  */
  bool compare_gimple_return (const greturn *s1, const greturn *s2);

  /* Verifies for given GIMPLEs S1 and S2 that
     goto statements are semantically equivalent.  */
  bool compare_gimple_goto (gimple *s1, gimple *s2);

  /* Verifies for given GIMPLE_RESX stmts S1 and S2 that
     resx statements are semantically equivalent.  */
  bool compare_gimple_resx (const gresx *s1, const gresx *s2);

  /* Verifies for given GIMPLE_ASM stmts S1 and S2 that ASM statements
     are equivalent.
     For the beginning, the pass only supports equality for
     '__asm__ __volatile__ ("", "", "", "memory")'.  */
  bool compare_gimple_asm (const gasm *s1, const gasm *s2);

  /* Verification function for declaration trees T1 and T2.  */
  bool compare_decl (const_tree t1, const_tree t2);

  /* Compute hash map MAP that determines loads and stores of STMT.  */
  enum operand_access_type {OP_MEMORY, OP_NORMAL};
  typedef hash_set<tree> operand_access_type_map;

  /* Function responsible for comparison of various operands T1 and T2.
     If these components, from functions FUNC1 and FUNC2, are equal, true
     is returned.  */
  bool compare_operand (tree t1, tree t2, operand_access_type type);

  /* Compares GIMPLE ASM inputs (or outputs) where we iterate tree chain
     and compare both TREE_PURPOSEs and TREE_VALUEs.  */
  bool compare_asm_inputs_outputs (tree t1, tree t2,
				   operand_access_type_map *map);

  /* Verifies that trees T1 and T2, representing function declarations
     are equivalent from perspective of ICF.  */
  bool compare_function_decl (tree t1, tree t2);

  /* Verifies that trees T1 and T2 do correspond.  */
  bool compare_variable_decl (const_tree t1, const_tree t2);

  /* Compare loop information for basic blocks BB1 and BB2.  */
  bool compare_loops (basic_block bb1, basic_block bb2);

  /* Return true if types are compatible for polymorphic call analysis.
     COMPARE_PTR indicates if polymorphic type comparison should be
     done for pointers, too.  */
  static bool compatible_polymorphic_types_p (tree t1, tree t2,
					      bool compare_ptr);

  /* Return true if types are compatible from perspective of ICF.
     FIRST_ARGUMENT indicates if the comparison is called for
     first parameter of a function.  */
  static bool compatible_types_p (tree t1, tree t2);

  /* Compute hash map determining access types of operands.  */
  static void classify_operands (const gimple *stmt,
				 operand_access_type_map *map);

  /* Return access type of a given operand.  */
  static operand_access_type get_operand_access_type
		 (operand_access_type_map *map, tree);
private:
  /* Vector mapping source SSA names to target ones.  */
  vec <int> m_source_ssa_names;

  /* Vector mapping target SSA names to source ones.  */
  vec <int> m_target_ssa_names;

  /* Source TREE function declaration.  */
  tree m_source_func_decl;

  /* Target TREE function declaration.  */
  tree m_target_func_decl;

  /* Source symbol nodes that should be skipped by
     declaration comparison.  */
  hash_set<symtab_node *> *m_ignored_source_nodes;

  /* Target symbol nodes that should be skipped by
     declaration comparison.  */
  hash_set<symtab_node *> *m_ignored_target_nodes;

  /* Source to target edge map.  */
  hash_map <edge, edge> m_edge_map;

  /* Source to target declaration map.  */
  hash_map <const_tree, const_tree> m_decl_map;

  /* Label to basic block index mapping.  */
  hash_map <const_tree, int> m_label_bb_map;

  /* Flag if ignore labels in comparison.  */
  bool m_ignore_labels;

  /* Flag if we should compare type based alias analysis info.  */
  bool m_tbaa;

public:
  /* Return true if two operands are equal.  The flags fields can be used
     to specify OEP flags described above.  */
  virtual bool operand_equal_p (const_tree, const_tree, unsigned int flags);

  /* Generate a hash value for an expression.  This can be used iteratively
     by passing a previous result as the HSTATE argument.  */
  virtual void hash_operand (const_tree, inchash::hash &, unsigned flags);
  void hash_operand (const_tree, inchash::hash &, unsigned flags,
		     operand_access_type access);
};

} // ipa_icf_gimple namespace
