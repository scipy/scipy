/* Information about function binary interfaces.
   Copyright (C) 2019-2022 Free Software Foundation, Inc.

This file is part of GCC

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

#ifndef GCC_FUNCTION_ABI_H
#define GCC_FUNCTION_ABI_H

/* Most targets use the same ABI for all functions in a translation
   unit, but some targets support interoperability between several ABIs.
   Each such ABI has a unique 0-based identifier, with 0 always being
   the default choice of ABI.

   NUM_ABI_IDS is the maximum number of such ABIs that GCC can handle at once.
   A bitfield with this number of bits can represent any combinaion of the
   supported ABIs.  */
const size_t NUM_ABI_IDS = 8;

/* Information about one of the target's predefined ABIs.  */
class predefined_function_abi
{
public:
  /* A target-specific identifier for this ABI.  The value must be in
     the range [0, NUM_ABI_IDS - 1].  */
  unsigned int id () const { return m_id; }

  /* True if this ABI has been initialized.  */
  bool initialized_p () const { return m_initialized; }

  /* Return true if a function call is allowed to alter every bit of
     register REGNO, so that the register contains an arbitrary value
     on return.  If so, the register cannot hold any part of a value
     that is live across a call.  */
  bool
  clobbers_full_reg_p (unsigned int regno) const
  {
    return TEST_HARD_REG_BIT (m_full_reg_clobbers, regno);
  }

  /* Return true if a function call is allowed to alter some or all bits
     of register REGNO.

     This is true whenever clobbers_full_reg_p (REGNO) is true.  It is
     also true if, for example, the ABI says that a call must preserve the
     low 32 or 64 bits of REGNO, but can clobber the upper bits of REGNO.
     In the latter case, it is possible for REGNO to hold values that
     are live across a call, provided that the value occupies only the
     call-preserved part of the register.  */
  bool
  clobbers_at_least_part_of_reg_p (unsigned int regno) const
  {
    return TEST_HARD_REG_BIT (m_full_and_partial_reg_clobbers, regno);
  }

  /* Return true if a function call is allowed to clobber at least part
     of (reg:MODE REGNO).  If so, it is not possible for the register
     as a whole to be live across a call.  */
  bool
  clobbers_reg_p (machine_mode mode, unsigned int regno) const
  {
    return overlaps_hard_reg_set_p (m_mode_clobbers[mode], mode, regno);
  }

  /* Return the set of registers that a function call is allowed to
     alter completely, so that the registers contain arbitrary values
     on return.  This doesn't include registers that a call can only
     partly clobber (as per TARGET_HARD_REGNO_CALL_PART_CLOBBERED).

     These registers cannot hold any part of a value that is live across
     a call.  */
  HARD_REG_SET full_reg_clobbers () const { return m_full_reg_clobbers; }

  /* Return the set of registers that a function call is allowed to alter
     to some degree.  For example, if an ABI says that a call must preserve
     the low 32 or 64 bits of a register R, but can clobber the upper bits
     of R, R would be in this set but not in full_reg_clobbers ().

     This set is a superset of full_reg_clobbers ().  It is possible for a
     register in full_and_partial_reg_clobbers () & ~full_reg_clobbers ()
     to contain values that are live across a call, provided that the live
     value only occupies the call-preserved part of the register.  */
  HARD_REG_SET
  full_and_partial_reg_clobbers () const
  {
    return m_full_and_partial_reg_clobbers;
  }

  /* Return the set of registers that cannot be used to hold a value of
     mode MODE across a function call.  That is:

       (reg:REGNO MODE)

     might be clobbered by a call whenever:

       overlaps_hard_reg_set (mode_clobbers (MODE), MODE, REGNO)

     In allocation terms, the registers in the returned set conflict
     with any value of mode MODE that is live across a call.  */
  HARD_REG_SET
  mode_clobbers (machine_mode mode) const
  {
    return m_mode_clobbers[mode];
  }

  void initialize (unsigned int, const_hard_reg_set);
  void add_full_reg_clobber (unsigned int);

private:
  unsigned int m_id : NUM_ABI_IDS;
  unsigned int m_initialized : 1;
  HARD_REG_SET m_full_reg_clobbers;
  HARD_REG_SET m_full_and_partial_reg_clobbers;
  HARD_REG_SET m_mode_clobbers[NUM_MACHINE_MODES];
};

/* Describes either a predefined ABI or the ABI of a particular function.
   In the latter case, the ABI might make use of extra function-specific
   information, such as for -fipa-ra.  */
class function_abi
{
public:
  /* Initialize the structure for a general function with the given ABI.  */
  function_abi (const predefined_function_abi &base_abi)
    : m_base_abi (&base_abi),
      m_mask (base_abi.full_and_partial_reg_clobbers ()) {}

  /* Initialize the structure for a function that has the given ABI and
     that is known not to clobber registers outside MASK.  */
  function_abi (const predefined_function_abi &base_abi,
		const_hard_reg_set mask)
    : m_base_abi (&base_abi), m_mask (mask) {}

  /* The predefined ABI from which this ABI is derived.  */
  const predefined_function_abi &base_abi () const { return *m_base_abi; }

  /* The target-specific identifier of the predefined ABI.  */
  unsigned int id () const { return m_base_abi->id (); }

  /* See the corresponding predefined_function_abi functions for
     details about the following functions.  */

  HARD_REG_SET
  full_reg_clobbers () const
  {
    return m_mask & m_base_abi->full_reg_clobbers ();
  }

  HARD_REG_SET
  full_and_partial_reg_clobbers () const
  {
    return m_mask & m_base_abi->full_and_partial_reg_clobbers ();
  }

  HARD_REG_SET
  mode_clobbers (machine_mode mode) const
  {
    return m_mask & m_base_abi->mode_clobbers (mode);
  }

  bool
  clobbers_full_reg_p (unsigned int regno) const
  {
    return (TEST_HARD_REG_BIT (m_mask, regno)
	    & m_base_abi->clobbers_full_reg_p (regno));
  }

  bool
  clobbers_at_least_part_of_reg_p (unsigned int regno) const
  {
    return (TEST_HARD_REG_BIT (m_mask, regno)
	    & m_base_abi->clobbers_at_least_part_of_reg_p (regno));
  }

  bool
  clobbers_reg_p (machine_mode mode, unsigned int regno) const
  {
    return overlaps_hard_reg_set_p (mode_clobbers (mode), mode, regno);
  }

  bool
  operator== (const function_abi &other) const
  {
    return m_base_abi == other.m_base_abi && m_mask == other.m_mask;
  }

  bool
  operator!= (const function_abi &other) const
  {
    return !operator== (other);
  }

protected:
  const predefined_function_abi *m_base_abi;
  HARD_REG_SET m_mask;
};

/* This class collects information about the ABIs of functions that are
   called in a particular region of code.  It is mostly intended to be
   used as a local variable during an IR walk.  */
class function_abi_aggregator
{
public:
  function_abi_aggregator () : m_abi_clobbers () {}

  /* Record that the code region calls a function with the given ABI.  */
  void
  note_callee_abi (const function_abi &abi)
  {
    m_abi_clobbers[abi.id ()] |= abi.full_and_partial_reg_clobbers ();
  }

  HARD_REG_SET caller_save_regs (const function_abi &) const;

private:
  HARD_REG_SET m_abi_clobbers[NUM_ABI_IDS];
};

struct target_function_abi_info
{
  /* An array of all the target ABIs that are available in this
     translation unit.  Not all entries are used for all targets,
     but the structures are relatively small, and using a fixed-size
     array avoids extra indirection.

     There are various ways of getting an ABI descriptor:

     * fndecl_abi (FNDECL) is the ABI of function FNDECL.

     * fntype_abi (FNTYPE) is the ABI of a function with type FNTYPE.

     * crtl->abi is the ABI of the function that we are currently
       compiling to rtl.

     * insn_callee_abi (INSN) is the ABI used by the target of call insn INSN.

     * eh_edge_abi is the "ABI" used when taking an EH edge from an
       exception-throwing statement to an exception handler.  Catching
       exceptions from calls can be treated as an abnormal return from
       those calls, and this ABI therefore describes the ABI of functions
       on such an abnormal return.  Statements that throw non-call
       exceptions can be treated as being implicitly wrapped in a call
       that has such an abnormal return.

       At present, no target needs to support more than one EH ABI.

     * function_abis[N] is the ABI with identifier N.  This can be useful
       when referring back to ABIs that have been collected by number in
       a bitmask, such as after walking function calls in a particular
       region of code.

     * default_function_abi refers specifically to the target's default
       choice of ABI, regardless of which (if any) functions actually
       use it.  This ABI and data derived from it do *not* provide
       globally conservatively-correct information, so it is only
       useful in very specific circumstances.  */
  predefined_function_abi x_function_abis[NUM_ABI_IDS];
};

extern target_function_abi_info default_target_function_abi_info;
#if SWITCHABLE_TARGET
extern target_function_abi_info *this_target_function_abi_info;
#else
#define this_target_function_abi_info (&default_target_function_abi_info)
#endif

/* See the comment above x_function_abis for when these macros should be used.
   At present, eh_edge_abi is always the default ABI, but that could change
   in future if a target needs it to.  */
#define function_abis \
  (this_target_function_abi_info->x_function_abis)
#define default_function_abi \
  (this_target_function_abi_info->x_function_abis[0])
#define eh_edge_abi default_function_abi

extern HARD_REG_SET call_clobbers_in_region (unsigned int, const_hard_reg_set,
					     machine_mode mode);

/* Return true if (reg:MODE REGNO) might be clobbered by one of the
   calls in a region described by ABIS and MASK, where:

   * Bit ID of ABIS is set if the region contains a call with
     function_abi identifier ID.

   * MASK contains all the registers that are fully or partially
     clobbered by calls in the region.

   This is not quite as accurate as testing each individual call,
   but it's a close and conservatively-correct approximation.
   It's much better for some targets than:

     overlaps_hard_reg_set_p (MASK, MODE, REGNO).  */

inline bool
call_clobbered_in_region_p (unsigned int abis, const_hard_reg_set mask,
			    machine_mode mode, unsigned int regno)
{
  HARD_REG_SET clobbers = call_clobbers_in_region (abis, mask, mode);
  return overlaps_hard_reg_set_p (clobbers, mode, regno);
}

extern const predefined_function_abi &fntype_abi (const_tree);
extern function_abi fndecl_abi (const_tree);
extern function_abi insn_callee_abi (const rtx_insn *);
extern function_abi expr_callee_abi (const_tree);

#endif
