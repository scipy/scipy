/* IPA predicates.
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

/* Representation of inline parameters that do depend on context function is
   inlined into (i.e. known constant values of function parameters.

   Conditions that are interesting for function body are collected into CONDS
   vector.  They are of simple as kind of a mathematical transformation on
   function parameter, T(function_param), in which the parameter occurs only
   once, and other operands are IPA invariant.  The conditions are then
   referred by predicates.  */


/* A simplified representation of tree node, for unary, binary and ternary
   operation.  Computations on parameter are decomposed to a series of this
   kind of structure.  */
struct GTY(()) expr_eval_op
{
  /* Result type of expression.  */
  tree type;
  /* Constant operands in expression, there are at most two.  */
  tree val[2];
  /* Index of parameter operand in expression.  */
  unsigned index : 2;
  /* Operation code of expression.  */
  ENUM_BITFIELD(tree_code) code : 16;
};

typedef vec<expr_eval_op, va_gc> *expr_eval_ops;

struct GTY(()) condition
{
  /* If agg_contents is set, this is the offset from which the used data was
     loaded.  */
  HOST_WIDE_INT offset;
  /* Type of the access reading the data (or the PARM_DECL SSA_NAME).  */
  tree type;
  tree val;
  int operand_num;
  ENUM_BITFIELD(tree_code) code : 16;
  /* Set if the used data were loaded from an aggregate parameter or from
     data received by reference.  */
  unsigned agg_contents : 1;
  /* If agg_contents is set, this differentiates between loads from data
     passed by reference and by value.  */
  unsigned by_ref : 1;
  /* A set of sequential operations on the parameter, which can be seen as
     a mathematical function on the parameter.  */
  expr_eval_ops param_ops;
};

/* Information kept about parameter of call site.  */
struct inline_param_summary
{
  /* REG_BR_PROB_BASE based probability that parameter will change in between
     two invocation of the calls.
     I.e. loop invariant parameters
     REG_BR_PROB_BASE/estimated_iterations and regular
     parameters REG_BR_PROB_BASE.

     Value 0 is reserved for compile time invariants. */
  short change_prob;
  unsigned points_to_local_or_readonly_memory : 1;
  bool equal_to (const inline_param_summary &other) const
  {
    return change_prob == other.change_prob
	   && points_to_local_or_readonly_memory
	      == other.points_to_local_or_readonly_memory;
  }
  bool useless_p (void) const
  {
    return change_prob == REG_BR_PROB_BASE
	   && !points_to_local_or_readonly_memory;
  }
};

typedef vec<condition, va_gc> *conditions;

/* Predicates are used to represent function parameters (such as runtime)
   which depend on a context function is called in.

   Predicates are logical formulas in conjunctive-disjunctive form consisting
   of clauses which are bitmaps specifying a set of condition that must
   be true for a clause to be satisfied. Physically they are represented as
   array of clauses terminated by 0.

   In order to make predicate (possibly) true, all of its clauses must
   be (possibly) true. To make clause (possibly) true, one of conditions
   it mentions must be (possibly) true.

   There are fixed bounds on number of clauses and conditions and all the
   manipulation functions are conservative in positive direction. I.e. we
   may lose precision by thinking that predicate may be true even when it
   is not.  */

typedef uint32_t clause_t;
class ipa_predicate
{
public:
  enum predicate_conditions
    {
      false_condition = 0,
      not_inlined_condition = 1,
      first_dynamic_condition = 2
    };

  /* Maximal number of conditions predicate can refer to.  This is limited
     by using clause_t to be 32bit.  */
  static const int num_conditions = 32;

  /* Special condition code we use to represent test that operand is compile
     time constant.  */
  static const tree_code is_not_constant = ERROR_MARK;

  /* Special condition code we use to represent test that operand is not changed
     across invocation of the function.  When operand IS_NOT_CONSTANT it is
     always CHANGED, however i.e. loop invariants can be NOT_CHANGED given
     percentage of executions even when they are not compile time constants.  */
  static const tree_code changed = IDENTIFIER_NODE;



  /* Initialize predicate either to true of false depending on P.  */
  inline ipa_predicate (bool p = true)
    {
      if (p)
        /* True predicate.  */
        m_clause[0] = 0;
      else
        /* False predicate. */
        set_to_cond (false_condition);
    }

  /* Sanity check that we do not mix pointers to predicates with predicates.  */
  inline ipa_predicate (ipa_predicate *)
    {
      gcc_unreachable ();
    }

  /* Return predicate testing condition I.  */
  static inline ipa_predicate predicate_testing_cond (int i)
    {
      ipa_predicate p;
      p.set_to_cond (i + first_dynamic_condition);
      return p;
    }

  /* Return predicate testing that function was not inlined.  */
  static ipa_predicate not_inlined (void)
    {
      ipa_predicate p;
      p.set_to_cond (not_inlined_condition);
      return p;
    }

  /* Compute logical and of ipa_predicates.  */
  ipa_predicate & operator &= (const ipa_predicate &);
  inline ipa_predicate operator &(const ipa_predicate &p) const
    {
      ipa_predicate ret = *this;
      ret &= p;
      return ret;
    }

  /* Compute logical or of ipa_predicates.  This is not operator because
     extra parameter CONDITIONS is needed  */
  ipa_predicate or_with (conditions, const ipa_predicate &) const;

  /* Return true if ipa_predicates are known to be equal.  */
  inline bool operator==(const ipa_predicate &p2) const
    {
      int i;
      for (i = 0; m_clause[i]; i++)
	{
	  gcc_checking_assert (i < max_clauses);
	  gcc_checking_assert (m_clause[i] > m_clause[i + 1]);
	  gcc_checking_assert (!p2.m_clause[i]
			       || p2.m_clause[i] > p2.m_clause[i + 1]);
	  if (m_clause[i] != p2.m_clause[i])
	    return false;
	}
      return !p2.m_clause[i];
    }

  /* Return true if predicates are known to be true or false depending
     on COND.  */
  inline bool operator==(const bool cond) const
    {
      if (cond)
        return !m_clause[0];
      if (m_clause[0] == (1 << false_condition))
	{
	  gcc_checking_assert (!m_clause[1]
			       && m_clause[0] == 1
				  << false_condition);
	  return true;
	}
      return false;
    }

  inline bool operator!=(const ipa_predicate &p2) const
    {
      return !(*this == p2);
    }

  inline bool operator!=(const bool cond) const
    {
      return !(*this == cond);
    }

  /* Evaluate if predicate is known to be false given the clause of possible
     truths.  */
  bool evaluate (clause_t) const;

  /* Estimate probability that predicate will be true in a given context.  */
  int probability (conditions, clause_t, vec<inline_param_summary>) const;

  /* Dump predicate to F. Output newline if nl.  */
  void dump (FILE *f, conditions, bool nl=true) const;
  void DEBUG_FUNCTION debug (conditions) const;

  /* Return ipa_predicate equal to THIS after duplication.  */
  ipa_predicate remap_after_duplication (clause_t);

  /* Return ipa_predicate equal to THIS after inlining.  */
  ipa_predicate remap_after_inlining (class ipa_fn_summary *,
				      ipa_node_params *params_summary,
				      ipa_fn_summary *,
				      const vec<int> &,
				      const vec<HOST_WIDE_INT> &,
				      clause_t, const ipa_predicate &);

  void stream_in (lto_input_block *);
  void stream_out (output_block *);

private:
  static const int max_clauses = 8;
  clause_t m_clause[max_clauses + 1];

  /* Initialize predicate to one testing single condition number COND.  */
  inline void set_to_cond (int cond)
    {
      m_clause[0] = 1 << cond;
      m_clause[1] = 0;
    }

  void add_clause (conditions conditions, clause_t);
};

void dump_condition (FILE *f, conditions conditions, int cond);
ipa_predicate add_condition (ipa_fn_summary *summary,
			     ipa_node_params *params_summary,
			     int operand_num,
			     tree type, struct agg_position_info *aggpos,
			     enum tree_code code, tree val,
			     expr_eval_ops param_ops = NULL);
