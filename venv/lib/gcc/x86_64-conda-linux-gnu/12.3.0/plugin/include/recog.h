/* Declarations for interface to insn recognizer and insn-output.cc.
   Copyright (C) 1987-2022 Free Software Foundation, Inc.

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

#ifndef GCC_RECOG_H
#define GCC_RECOG_H

/* Random number that should be large enough for all purposes.  Also define
   a type that has at least MAX_RECOG_ALTERNATIVES + 1 bits, with the extra
   bit giving an invalid value that can be used to mean "uninitialized".  */
#define MAX_RECOG_ALTERNATIVES 35
typedef uint64_t alternative_mask;  /* Keep in sync with genattrtab.cc.  */

/* A mask of all alternatives.  */
#define ALL_ALTERNATIVES ((alternative_mask) -1)

/* A mask containing just alternative X.  */
#define ALTERNATIVE_BIT(X) ((alternative_mask) 1 << (X))

/* Types of operands.  */
enum op_type {
  OP_IN,
  OP_OUT,
  OP_INOUT
};

struct operand_alternative
{
  /* Pointer to the beginning of the constraint string for this alternative,
     for easier access by alternative number.  */
  const char *constraint;

  /* The register class valid for this alternative (possibly NO_REGS).  */
  ENUM_BITFIELD (reg_class) cl : 16;

  /* "Badness" of this alternative, computed from number of '?' and '!'
     characters in the constraint string.  */
  unsigned int reject : 16;

  /* -1 if no matching constraint was found, or an operand number.  */
  int matches : 8;
  /* The same information, but reversed: -1 if this operand is not
     matched by any other, or the operand number of the operand that
     matches this one.  */
  int matched : 8;

  /* Nonzero if '&' was found in the constraint string.  */
  unsigned int earlyclobber : 1;
  /* Nonzero if TARGET_MEM_CONSTRAINT was found in the constraint
     string.  */
  unsigned int memory_ok : 1;
  /* Nonzero if 'p' was found in the constraint string.  */
  unsigned int is_address : 1;
  /* Nonzero if 'X' was found in the constraint string, or if the constraint
     string for this alternative was empty.  */
  unsigned int anything_ok : 1;

  unsigned int unused : 12;
};

/* Return the class for operand I of alternative ALT, taking matching
   constraints into account.  */

static inline enum reg_class
alternative_class (const operand_alternative *alt, int i)
{
  return alt[i].matches >= 0 ? alt[alt[i].matches].cl : alt[i].cl;
}

/* A class for substituting one rtx for another within an instruction,
   or for recursively simplifying the instruction as-is.  Derived classes
   can record or filter certain decisions.  */

class insn_propagation : public simplify_context
{
public:
  /* Assignments for RESULT_FLAGS.

     UNSIMPLIFIED is true if a substitution has been made inside an rtx
     X and if neither X nor its parent expressions could be simplified.

     FIRST_SPARE_RESULT is the first flag available for derived classes.  */
  static const uint16_t UNSIMPLIFIED = 1U << 0;
  static const uint16_t FIRST_SPARE_RESULT = 1U << 1;

  insn_propagation (rtx_insn *);
  insn_propagation (rtx_insn *, rtx, rtx, bool = true);
  bool apply_to_pattern (rtx *);
  bool apply_to_rvalue (rtx *);

  /* Return true if we should accept a substitution into the address of
     memory expression MEM.  Undoing changes OLD_NUM_CHANGES and up restores
     MEM's original address.  */
  virtual bool check_mem (int /*old_num_changes*/,
			  rtx /*mem*/) { return true; }

  /* Note that we've simplified OLD_RTX into NEW_RTX.  When substituting,
     this only happens if a substitution occured within OLD_RTX.
     Undoing OLD_NUM_CHANGES and up will restore the old form of OLD_RTX.
     OLD_RESULT_FLAGS is the value that RESULT_FLAGS had before processing
     OLD_RTX.  */
  virtual void note_simplification (int /*old_num_changes*/,
				    uint16_t /*old_result_flags*/,
				    rtx /*old_rtx*/, rtx /*new_rtx*/) {}

private:
  bool apply_to_mem_1 (rtx);
  bool apply_to_lvalue_1 (rtx);
  bool apply_to_rvalue_1 (rtx *);
  bool apply_to_pattern_1 (rtx *);

public:
  /* The instruction that we are simplifying or propagating into.  */
  rtx_insn *insn;

  /* If FROM is nonnull, we're replacing FROM with TO, otherwise we're
     just doing a recursive simplification.  */
  rtx from;
  rtx to;

  /* The number of times that we have replaced FROM with TO.  */
  unsigned int num_replacements;

  /* A bitmask of flags that describe the result of the simplificiation;
     see above for details.  */
  uint16_t result_flags : 16;

  /* True if we should unshare TO when making the next substitution,
     false if we can use TO itself.  */
  uint16_t should_unshare : 1;

  /* True if we should call check_mem after substituting into a memory.  */
  uint16_t should_check_mems : 1;

  /* True if we should call note_simplification after each simplification.  */
  uint16_t should_note_simplifications : 1;

  /* For future expansion.  */
  uint16_t spare : 13;

  /* Gives the reason that a substitution failed, for debug purposes.  */
  const char *failure_reason;
};

/* Try to replace FROM with TO in INSN.  SHARED_P is true if TO is shared
   with other instructions, false if INSN can use TO directly.  */

inline insn_propagation::insn_propagation (rtx_insn *insn, rtx from, rtx to,
					   bool shared_p)
  : insn (insn),
    from (from),
    to (to),
    num_replacements (0),
    result_flags (0),
    should_unshare (shared_p),
    should_check_mems (false),
    should_note_simplifications (false),
    spare (0),
    failure_reason (nullptr)
{
}

/* Try to simplify INSN without performing a substitution.  */

inline insn_propagation::insn_propagation (rtx_insn *insn)
  : insn_propagation (insn, NULL_RTX, NULL_RTX)
{
}

extern void init_recog (void);
extern void init_recog_no_volatile (void);
extern int check_asm_operands (rtx);
extern int asm_operand_ok (rtx, const char *, const char **);
extern bool validate_change (rtx, rtx *, rtx, bool);
extern bool validate_unshare_change (rtx, rtx *, rtx, bool);
extern bool validate_change_xveclen (rtx, rtx *, int, bool);
extern bool canonicalize_change_group (rtx_insn *insn, rtx x);
extern int insn_invalid_p (rtx_insn *, bool);
extern int verify_changes (int);
extern void confirm_change_group (void);
extern int apply_change_group (void);
extern int num_validated_changes (void);
extern void cancel_changes (int);
extern void temporarily_undo_changes (int);
extern void redo_changes (int);
extern int constrain_operands (int, alternative_mask);
extern int constrain_operands_cached (rtx_insn *, int);
extern bool memory_address_addr_space_p (machine_mode, rtx, addr_space_t);
#define memory_address_p(mode,addr) \
	memory_address_addr_space_p ((mode), (addr), ADDR_SPACE_GENERIC)
extern bool strict_memory_address_addr_space_p (machine_mode, rtx,
						addr_space_t);
#define strict_memory_address_p(mode,addr) \
	strict_memory_address_addr_space_p ((mode), (addr), ADDR_SPACE_GENERIC)
extern int validate_replace_rtx_subexp (rtx, rtx, rtx_insn *, rtx *);
extern int validate_replace_rtx (rtx, rtx, rtx_insn *);
extern int validate_replace_rtx_part (rtx, rtx, rtx *, rtx_insn *);
extern int validate_replace_rtx_part_nosimplify (rtx, rtx, rtx *, rtx_insn *);
extern void validate_replace_rtx_group (rtx, rtx, rtx_insn *);
extern void validate_replace_src_group (rtx, rtx, rtx_insn *);
extern bool validate_simplify_insn (rtx_insn *insn);
extern int num_changes_pending (void);
extern bool reg_fits_class_p (const_rtx, reg_class_t, int, machine_mode);
extern bool valid_insn_p (rtx_insn *);

extern bool offsettable_memref_p (rtx);
extern bool offsettable_nonstrict_memref_p (rtx);
extern bool offsettable_address_addr_space_p (int, machine_mode, rtx,
					     addr_space_t);
#define offsettable_address_p(strict,mode,addr) \
	offsettable_address_addr_space_p ((strict), (mode), (addr), \
					  ADDR_SPACE_GENERIC)
extern bool mode_dependent_address_p (rtx, addr_space_t);

extern int recog (rtx, rtx_insn *, int *);
#ifndef GENERATOR_FILE
static inline int recog_memoized (rtx_insn *insn);
#endif
extern void add_clobbers (rtx, int);
extern int added_clobbers_hard_reg_p (int);
extern void insn_extract (rtx_insn *);
extern void extract_insn (rtx_insn *);
extern void extract_constrain_insn (rtx_insn *insn);
extern void extract_constrain_insn_cached (rtx_insn *);
extern void extract_insn_cached (rtx_insn *);
extern void preprocess_constraints (int, int, const char **,
				    operand_alternative *, rtx **);
extern const operand_alternative *preprocess_insn_constraints (unsigned int);
extern void preprocess_constraints (rtx_insn *);
extern rtx_insn *peep2_next_insn (int);
extern int peep2_regno_dead_p (int, int);
extern int peep2_reg_dead_p (int, rtx);
#ifdef HARD_CONST
extern rtx peep2_find_free_register (int, int, const char *,
				     machine_mode, HARD_REG_SET *);
#endif
extern rtx_insn *peephole2_insns (rtx, rtx_insn *, int *);

extern int store_data_bypass_p (rtx_insn *, rtx_insn *);
extern int if_test_bypass_p (rtx_insn *, rtx_insn *);

extern void copy_frame_info_to_split_insn (rtx_insn *, rtx_insn *);

#ifndef GENERATOR_FILE
/* Try recognizing the instruction INSN,
   and return the code number that results.
   Remember the code so that repeated calls do not
   need to spend the time for actual rerecognition.

   This function is the normal interface to instruction recognition.
   The automatically-generated function `recog' is normally called
   through this one.  */

static inline int
recog_memoized (rtx_insn *insn)
{
  if (INSN_CODE (insn) < 0)
    INSN_CODE (insn) = recog (PATTERN (insn), insn, 0);
  return INSN_CODE (insn);
}
#endif

/* Skip chars until the next ',' or the end of the string.  This is
   useful to skip alternatives in a constraint string.  */
static inline const char *
skip_alternative (const char *p)
{
  const char *r = p;
  while (*r != '\0' && *r != ',')
    r++;
  if (*r == ',')
    r++;
  return r;
}

/* Nonzero means volatile operands are recognized.  */
extern int volatile_ok;

/* RAII class for temporarily setting volatile_ok.  */

class temporary_volatile_ok
{
public:
  temporary_volatile_ok (int value) : save_volatile_ok (volatile_ok)
  {
    volatile_ok = value;
  }

  ~temporary_volatile_ok () { volatile_ok = save_volatile_ok; }

private:
  temporary_volatile_ok (const temporary_volatile_ok &);
  int save_volatile_ok;
};

/* Set by constrain_operands to the number of the alternative that
   matched.  */
extern int which_alternative;

/* The following vectors hold the results from insn_extract.  */

struct recog_data_d
{
  /* It is very tempting to make the 5 operand related arrays into a
     structure and index on that.  However, to be source compatible
     with all of the existing md file insn constraints and output
     templates, we need `operand' as a flat array.  Without that
     member, making an array for the rest seems pointless.  */

  /* Gives value of operand N.  */
  rtx operand[MAX_RECOG_OPERANDS];

  /* Gives location where operand N was found.  */
  rtx *operand_loc[MAX_RECOG_OPERANDS];

  /* Gives the constraint string for operand N.  */
  const char *constraints[MAX_RECOG_OPERANDS];

  /* Nonzero if operand N is a match_operator or a match_parallel.  */
  char is_operator[MAX_RECOG_OPERANDS];

  /* Gives the mode of operand N.  */
  machine_mode operand_mode[MAX_RECOG_OPERANDS];

  /* Gives the type (in, out, inout) for operand N.  */
  enum op_type operand_type[MAX_RECOG_OPERANDS];

  /* Gives location where the Nth duplicate-appearance of an operand
     was found.  This is something that matched MATCH_DUP.  */
  rtx *dup_loc[MAX_DUP_OPERANDS];

  /* Gives the operand number that was duplicated in the Nth
     duplicate-appearance of an operand.  */
  char dup_num[MAX_DUP_OPERANDS];

  /* ??? Note that these are `char' instead of `unsigned char' to (try to)
     avoid certain lossage from K&R C, wherein `unsigned char' default
     promotes to `unsigned int' instead of `int' as in ISO C.  As of 1999,
     the most common places to bootstrap from K&R C are SunOS and HPUX,
     both of which have signed characters by default.  The only other
     supported natives that have both K&R C and unsigned characters are
     ROMP and Irix 3, and neither have been seen for a while, but do
     continue to consider unsignedness when performing arithmetic inside
     a comparison.  */

  /* The number of operands of the insn.  */
  char n_operands;

  /* The number of MATCH_DUPs in the insn.  */
  char n_dups;

  /* The number of alternatives in the constraints for the insn.  */
  char n_alternatives;

  /* True if insn is ASM_OPERANDS.  */
  bool is_asm;

  /* In case we are caching, hold insn data was generated for.  */
  rtx_insn *insn;
};

extern struct recog_data_d recog_data;

extern const operand_alternative *recog_op_alt;

/* Return a pointer to an array in which index OP describes the constraints
   on operand OP of the current instruction alternative (which_alternative).
   Only valid after calling preprocess_constraints and constrain_operands.  */

inline static const operand_alternative *
which_op_alt ()
{
  gcc_checking_assert (IN_RANGE (which_alternative, 0,
				 recog_data.n_alternatives - 1));
  return &recog_op_alt[which_alternative * recog_data.n_operands];
}

/* A table defined in insn-output.cc that give information about
   each insn-code value.  */

typedef bool (*insn_operand_predicate_fn) (rtx, machine_mode);
typedef const char * (*insn_output_fn) (rtx *, rtx_insn *);

struct insn_gen_fn
{
  typedef void (*stored_funcptr) (void);

  template<typename ...Ts>
  rtx_insn *operator() (Ts... args) const
  {
    typedef rtx_insn *(*funcptr) (decltype ((void) args, NULL_RTX)...);
    return ((funcptr) func) (args...);
  }

  // This is for compatibility of code that invokes functions like
  //   (*funcptr) (arg)
  insn_gen_fn operator * (void) const { return *this; }

  // The wrapped function pointer must be public and there must not be any
  // constructors.  Otherwise the insn_data_d struct initializers generated
  // by genoutput.cc will result in static initializer functions, which defeats
  // the purpose of the generated insn_data_d array.
  stored_funcptr func;
};

struct insn_operand_data
{
  const insn_operand_predicate_fn predicate;

  const char *const constraint;

  ENUM_BITFIELD(machine_mode) const mode : 16;

  const char strict_low;

  const char is_operator;

  const char eliminable;

  const char allows_mem;
};

/* Legal values for insn_data.output_format.  Indicate what type of data
   is stored in insn_data.output.  */
#define INSN_OUTPUT_FORMAT_NONE		0	/* abort */
#define INSN_OUTPUT_FORMAT_SINGLE	1	/* const char * */
#define INSN_OUTPUT_FORMAT_MULTI	2	/* const char * const * */
#define INSN_OUTPUT_FORMAT_FUNCTION	3	/* const char * (*)(...) */

struct insn_data_d
{
  const char *const name;
#if HAVE_DESIGNATED_UNION_INITIALIZERS
  union {
    const char *single;
    const char *const *multi;
    insn_output_fn function;
  } output;
#else
  struct {
    const char *single;
    const char *const *multi;
    insn_output_fn function;
  } output;
#endif
  const insn_gen_fn genfun;
  const struct insn_operand_data *const operand;

  const char n_generator_args;
  const char n_operands;
  const char n_dups;
  const char n_alternatives;
  const char output_format;
};

extern const struct insn_data_d insn_data[];
extern int peep2_current_count;

#ifndef GENERATOR_FILE
#include "insn-codes.h"

/* An enum of boolean attributes that may only depend on the current
   subtarget, not on things like operands or compiler phase.  */
enum bool_attr {
  BA_ENABLED,
  BA_PREFERRED_FOR_SPEED,
  BA_PREFERRED_FOR_SIZE,
  BA_LAST = BA_PREFERRED_FOR_SIZE
};

/* Target-dependent globals.  */
struct target_recog {
  bool x_initialized;
  alternative_mask x_bool_attr_masks[NUM_INSN_CODES][BA_LAST + 1];
  operand_alternative *x_op_alt[NUM_INSN_CODES];
};

extern struct target_recog default_target_recog;
#if SWITCHABLE_TARGET
extern struct target_recog *this_target_recog;
#else
#define this_target_recog (&default_target_recog)
#endif

alternative_mask get_enabled_alternatives (rtx_insn *);
alternative_mask get_preferred_alternatives (rtx_insn *);
alternative_mask get_preferred_alternatives (rtx_insn *, basic_block);
bool check_bool_attrs (rtx_insn *);

void recog_init ();

/* This RAII class can help to undo tentative insn changes on failure.
   When an object of the class goes out of scope, it undoes all group
   changes that have been made via the validate_change machinery and
   not yet confirmed via confirm_change_group.

   For example:

      insn_change_watermark watermark;
      validate_change (..., true); // A
      ...
      if (test)
	// Undoes change A.
	return false;
      ...
      validate_change (..., true); // B
      ...
      if (test)
	// Undoes changes A and B.
	return false;
      ...
      confirm_change_group ();

   Code that wants to avoid this behavior can use keep ():

      insn_change_watermark watermark;
      validate_change (..., true); // A
      ...
      if (test)
	// Undoes change A.
	return false;
      ...
      watermark.keep ();
      validate_change (..., true); // B
      ...
      if (test)
	// Undoes change B, but not A.
	return false;
      ...
      confirm_change_group ();  */
class insn_change_watermark
{
public:
  insn_change_watermark () : m_old_num_changes (num_validated_changes ()) {}
  ~insn_change_watermark ();
  void keep () { m_old_num_changes = num_validated_changes (); }

private:
  int m_old_num_changes;
};

inline insn_change_watermark::~insn_change_watermark ()
{
  if (m_old_num_changes < num_validated_changes ())
    cancel_changes (m_old_num_changes);
}

#endif

#endif /* GCC_RECOG_H */
