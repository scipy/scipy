/* Support for simple predicate analysis.

   Copyright (C) 2021-2022 Free Software Foundation, Inc.
   Contributed by Martin Sebor <msebor@redhat.com>

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

#ifndef GIMPLE_PREDICATE_ANALYSIS_H_INCLUDED
#define GIMPLE_PREDICATE_ANALYSIS_H_INCLUDED

#define MAX_NUM_CHAINS 8
#define MAX_CHAIN_LEN 5
#define MAX_POSTDOM_CHECK 8
#define MAX_SWITCH_CASES 40

/* Represents a simple Boolean predicate.  */
struct pred_info
{
  tree pred_lhs;
  tree pred_rhs;
  enum tree_code cond_code;
  bool invert;
};

/* The type to represent a sequence of predicates grouped
   with .AND. operation.  */
typedef vec<pred_info, va_heap, vl_ptr> pred_chain;

/* The type to represent a sequence of pred_chains grouped
   with .OR. operation.  */
typedef vec<pred_chain, va_heap, vl_ptr> pred_chain_union;

/* Represents a complex Boolean predicate expression.  */
class predicate
{
 public:
  /* Base function object type used to determine whether an expression
     is of interest.  */
  struct func_t
  {
    typedef unsigned phi_arg_set_t;

    /* Return true if the argument is an expression of interest.  */
    virtual bool operator()(tree) = 0;
    /* Return a bitset of PHI arguments of interest.  By default returns
       bitset with a bit set for each argument.  Should be called in
       the overriden function first and, if nonzero, the result then
       refined as appropriate.  */
    virtual phi_arg_set_t phi_arg_set (gphi *);

    /* Maximum number of PHI arguments supported by phi_arg_set().  */
    static constexpr unsigned max_phi_args =
      sizeof (phi_arg_set_t) * CHAR_BIT;
  };

  /* Construct with the specified EVAL object.  */
  predicate (func_t &eval)
    : m_preds (vNULL), m_eval (eval), m_use_expr () { }

  /* Copy.  */
  predicate (const predicate &rhs)
    : m_preds (vNULL), m_eval (rhs.m_eval), m_use_expr ()
    {
      *this = rhs;
    }

  predicate (basic_block, basic_block, func_t &);

  ~predicate ();

  /* Assign.  */
  predicate& operator= (const predicate &);

  bool is_empty () const
  {
    return m_preds.is_empty ();
  }

  const pred_chain_union chain () const
  {
    return m_preds;
  }

  /* Return true if the use by a statement in the basic block of
     a PHI operand is ruled out (i.e., guarded) by *THIS.  */
  bool is_use_guarded (gimple *, basic_block, gphi *, unsigned);

  void init_from_control_deps (const vec<edge> *, unsigned);

  void dump (gimple *, const char *) const;

  void normalize (gimple * = NULL, bool = false);
  void simplify (gimple * = NULL, bool = false);

  bool is_use_guarded (gimple *, basic_block, gphi *, unsigned,
		       hash_set<gphi *> *);

  /* Return the predicate expression guarding the definition of
     the interesting variable, optionally inverted.  */
  tree def_expr (bool = false) const;
  /* Return the predicate expression guarding the use of the interesting
     variable.  */
  tree use_expr () const;

  tree expr (bool = false) const;

private:
  bool includes (const pred_chain &) const;
  bool superset_of (const predicate &) const;
  bool overlap (gphi *, unsigned, hash_set<gphi *> *);
  bool use_cannot_happen (gphi *, unsigned);

  bool init_from_phi_def (gphi *);

  void push_pred (const pred_info &);

  /* Normalization functions.  */
  void normalize (pred_chain *, pred_info, tree_code, pred_chain *,
		  hash_set<tree> *);

  void normalize (const pred_info &);
  void normalize (const pred_chain &);

  /* Simplification functions.  */
  bool simplify_2 ();
  bool simplify_3 ();
  bool simplify_4 ();

private:
  /* Representation of the predicate expression(s).  */
  pred_chain_union m_preds;
  /* Callback to evaluate an operand.  Return true if it's interesting.  */
  func_t &m_eval;
  /* The predicate expression guarding the use of the interesting
     variable.  */
  tree m_use_expr;
};

/* Bit mask handling macros.  */
#define MASK_SET_BIT(mask, pos) mask |= (1 << pos)
#define MASK_TEST_BIT(mask, pos) (mask & (1 << pos))
#define MASK_EMPTY(mask) (mask == 0)

#endif // GIMPLE_PREDICATE_ANALYSIS_H_INCLUDED
