/* Profile counter container type.
   Copyright (C) 2017-2022 Free Software Foundation, Inc.
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

#ifndef GCC_PROFILE_COUNT_H
#define GCC_PROFILE_COUNT_H

struct function;
struct profile_count;
class sreal;

/* Quality of the profile count.  Because gengtype does not support enums
   inside of classes, this is in global namespace.  */
enum profile_quality {
  /* Uninitialized value.  */
  UNINITIALIZED_PROFILE,

  /* Profile is based on static branch prediction heuristics and may
     or may not match reality.  It is local to function and cannot be compared
     inter-procedurally.  Never used by probabilities (they are always local).
   */
  GUESSED_LOCAL,

  /* Profile was read by feedback and was 0, we used local heuristics to guess
     better.  This is the case of functions not run in profile feedback.
     Never used by probabilities.  */
  GUESSED_GLOBAL0,

  /* Same as GUESSED_GLOBAL0 but global count is adjusted 0.  */
  GUESSED_GLOBAL0_ADJUSTED,

  /* Profile is based on static branch prediction heuristics.  It may or may
     not reflect the reality but it can be compared interprocedurally
     (for example, we inlined function w/o profile feedback into function
      with feedback and propagated from that).
     Never used by probabilities.  */
  GUESSED,

  /* Profile was determined by autofdo.  */
  AFDO,

  /* Profile was originally based on feedback but it was adjusted
     by code duplicating optimization.  It may not precisely reflect the
     particular code path.  */
  ADJUSTED,

  /* Profile was read from profile feedback or determined by accurate static
     method.  */
  PRECISE
};

extern const char *profile_quality_as_string (enum profile_quality);
extern bool parse_profile_quality (const char *value,
				   profile_quality *quality);

/* The base value for branch probability notes and edge probabilities.  */
#define REG_BR_PROB_BASE  10000

#define RDIV(X,Y) (((X) + (Y) / 2) / (Y))

bool slow_safe_scale_64bit (uint64_t a, uint64_t b, uint64_t c, uint64_t *res);

/* Compute RES=(a*b + c/2)/c capping and return false if overflow happened.  */

inline bool
safe_scale_64bit (uint64_t a, uint64_t b, uint64_t c, uint64_t *res)
{
#if (GCC_VERSION >= 5000)
  uint64_t tmp;
  if (!__builtin_mul_overflow (a, b, &tmp)
      && !__builtin_add_overflow (tmp, c/2, &tmp))
    {
      *res = tmp / c;
      return true;
    }
  if (c == 1)
    {
      *res = (uint64_t) -1;
      return false;
    }
#else
  if (a < ((uint64_t)1 << 31)
      && b < ((uint64_t)1 << 31)
      && c < ((uint64_t)1 << 31))
    {
      *res = (a * b + (c / 2)) / c;
      return true;
    }
#endif
  return slow_safe_scale_64bit (a, b, c, res);
}

/* Data type to hold probabilities.  It implements fixed point arithmetics
   with capping so probability is always in range [0,1] and scaling requiring
   values greater than 1 needs to be represented otherwise.

   In addition to actual value the quality of profile is tracked and propagated
   through all operations.  Special value UNINITIALIZED_PROFILE is used for probabilities
   that has not been determined yet (for example because of
   -fno-guess-branch-probability)

   Typically probabilities are derived from profile feedback (via
   probability_in_gcov_type), autoFDO or guessed statically and then propagated
   thorough the compilation.

   Named probabilities are available:
     - never           (0 probability)
     - guessed_never
     - very_unlikely   (1/2000 probability)
     - unlikely        (1/5 probability)
     - even            (1/2 probability)
     - likely          (4/5 probability)
     - very_likely     (1999/2000 probability)
     - guessed_always
     - always

   Named probabilities except for never/always are assumed to be statically
   guessed and thus not necessarily accurate.  The difference between never
   and guessed_never is that the first one should be used only in case that
   well behaving program will very likely not execute the "never" path.
   For example if the path is going to abort () call or it exception handling.

   Always and guessed_always probabilities are symmetric.

   For legacy code we support conversion to/from REG_BR_PROB_BASE based fixpoint
   integer arithmetics. Once the code is converted to branch probabilities,
   these conversions will probably go away because they are lossy.
*/

class GTY((user)) profile_probability
{
  static const int n_bits = 29;
  /* We can technically use ((uint32_t) 1 << (n_bits - 1)) - 2 but that
     will lead to harder multiplication sequences.  */
  static const uint32_t max_probability = (uint32_t) 1 << (n_bits - 2);
  static const uint32_t uninitialized_probability
		 = ((uint32_t) 1 << (n_bits - 1)) - 1;

  uint32_t m_val : 29;
  enum profile_quality m_quality : 3;

  friend struct profile_count;
public:
  profile_probability (): m_val (uninitialized_probability),
    m_quality (GUESSED)
  {}

  profile_probability (uint32_t val, profile_quality quality):
    m_val (val), m_quality (quality)
  {}

  /* Named probabilities.  */
  static profile_probability never ()
    {
      profile_probability ret;
      ret.m_val = 0;
      ret.m_quality = PRECISE;
      return ret;
    }

  static profile_probability guessed_never ()
    {
      profile_probability ret;
      ret.m_val = 0;
      ret.m_quality = GUESSED;
      return ret;
    }

  static profile_probability very_unlikely ()
    {
      /* Be consistent with PROB_VERY_UNLIKELY in predict.h.  */
      profile_probability r = guessed_always ().apply_scale (1, 2000);
      r.m_val--;
      return r;
    }

  static profile_probability unlikely ()
    {
      /* Be consistent with PROB_VERY_LIKELY in predict.h.  */
      profile_probability r = guessed_always ().apply_scale (1, 5);
      r.m_val--;
      return r;
    }

  static profile_probability even ()
    {
      return guessed_always ().apply_scale (1, 2);
    }

  static profile_probability very_likely ()
    {
      return always () - very_unlikely ();
    }

  static profile_probability likely ()
    {
      return always () - unlikely ();
    }

  static profile_probability guessed_always ()
    {
      profile_probability ret;
      ret.m_val = max_probability;
      ret.m_quality = GUESSED;
      return ret;
    }

  static profile_probability always ()
    {
      profile_probability ret;
      ret.m_val = max_probability;
      ret.m_quality = PRECISE;
      return ret;
    }

  /* Probabilities which has not been initialized. Either because
     initialization did not happen yet or because profile is unknown.  */
  static profile_probability uninitialized ()
    {
      profile_probability c;
      c.m_val = uninitialized_probability;
      c.m_quality = GUESSED;
      return c;
    }

  /* Return true if value has been initialized.  */
  bool initialized_p () const
    {
      return m_val != uninitialized_probability;
    }

  /* Return true if value can be trusted.  */
  bool reliable_p () const
    {
      return m_quality >= ADJUSTED;
    }

  /* Conversion from and to REG_BR_PROB_BASE integer fixpoint arithmetics.
     this is mostly to support legacy code and should go away.  */
  static profile_probability from_reg_br_prob_base (int v)
    {
      profile_probability ret;
      gcc_checking_assert (v >= 0 && v <= REG_BR_PROB_BASE);
      ret.m_val = RDIV (v * (uint64_t) max_probability, REG_BR_PROB_BASE);
      ret.m_quality = GUESSED;
      return ret;
    }

  /* Return THIS with quality set to ADJUSTED.  */
  profile_probability adjusted () const
    {
      profile_probability ret = *this;
      if (!initialized_p ())
	return *this;
      ret.m_quality = ADJUSTED;
      return ret;
    }

  int to_reg_br_prob_base () const
    {
      gcc_checking_assert (initialized_p ());
      return RDIV (m_val * (uint64_t) REG_BR_PROB_BASE, max_probability);
    }

  /* Conversion to and from RTL representation of profile probabilities.  */
  static profile_probability from_reg_br_prob_note (int v)
    {
      profile_probability ret;
      ret.m_val = ((unsigned int)v) / 8;
      ret.m_quality = (enum profile_quality)(v & 7);
      return ret;
    }

  int to_reg_br_prob_note () const
    {
      gcc_checking_assert (initialized_p ());
      int ret = m_val * 8 + m_quality;
      gcc_checking_assert (from_reg_br_prob_note (ret) == *this);
      return ret;
    }

  /* Return VAL1/VAL2.  */
  static profile_probability probability_in_gcov_type
				 (gcov_type val1, gcov_type val2)
    {
      profile_probability ret;
      gcc_checking_assert (val1 >= 0 && val2 > 0);
      if (val1 > val2)
	ret.m_val = max_probability;
      else
	{
	  uint64_t tmp;
	  safe_scale_64bit (val1, max_probability, val2, &tmp);
	  gcc_checking_assert (tmp <= max_probability);
	  ret.m_val = tmp;
	}
      ret.m_quality = PRECISE;
      return ret;
    }

  /* Basic operations.  */
  bool operator== (const profile_probability &other) const
    {
      return m_val == other.m_val && m_quality == other.m_quality;
    }

  profile_probability operator+ (const profile_probability &other) const
    {
      if (other == never ())
	return *this;
      if (*this == never ())
	return other;
      if (!initialized_p () || !other.initialized_p ())
	return uninitialized ();

      profile_probability ret;
      ret.m_val = MIN ((uint32_t)(m_val + other.m_val), max_probability);
      ret.m_quality = MIN (m_quality, other.m_quality);
      return ret;
    }

  profile_probability &operator+= (const profile_probability &other)
    {
      if (other == never ())
	return *this;
      if (*this == never ())
	{
	  *this = other;
	  return *this;
	}
      if (!initialized_p () || !other.initialized_p ())
	return *this = uninitialized ();
      else
	{
	  m_val = MIN ((uint32_t)(m_val + other.m_val), max_probability);
	  m_quality = MIN (m_quality, other.m_quality);
	}
      return *this;
    }

  profile_probability operator- (const profile_probability &other) const
    {
      if (*this == never ()
	  || other == never ())
	return *this;
      if (!initialized_p () || !other.initialized_p ())
	return uninitialized ();
      profile_probability ret;
      ret.m_val = m_val >= other.m_val ? m_val - other.m_val : 0;
      ret.m_quality = MIN (m_quality, other.m_quality);
      return ret;
    }

  profile_probability &operator-= (const profile_probability &other)
    {
      if (*this == never ()
	  || other == never ())
	return *this;
      if (!initialized_p () || !other.initialized_p ())
	return *this = uninitialized ();
      else
	{
	  m_val = m_val >= other.m_val ? m_val - other.m_val : 0;
	  m_quality = MIN (m_quality, other.m_quality);
	}
      return *this;
    }

  profile_probability operator* (const profile_probability &other) const
    {
      if (*this == never ()
	  || other == never ())
	return never ();
      if (!initialized_p () || !other.initialized_p ())
	return uninitialized ();
      profile_probability ret;
      ret.m_val = RDIV ((uint64_t)m_val * other.m_val, max_probability);
      ret.m_quality = MIN (MIN (m_quality, other.m_quality), ADJUSTED);
      return ret;
    }

  profile_probability &operator*= (const profile_probability &other)
    {
      if (*this == never ()
	  || other == never ())
	return *this = never ();
      if (!initialized_p () || !other.initialized_p ())
	return *this = uninitialized ();
      else
	{
	  m_val = RDIV ((uint64_t)m_val * other.m_val, max_probability);
	  m_quality = MIN (MIN (m_quality, other.m_quality), ADJUSTED);
	}
      return *this;
    }

  profile_probability operator/ (const profile_probability &other) const
    {
      if (*this == never ())
	return never ();
      if (!initialized_p () || !other.initialized_p ())
	return uninitialized ();
      profile_probability ret;
      /* If we get probability above 1, mark it as unreliable and return 1. */
      if (m_val >= other.m_val)
	{
	  ret.m_val = max_probability;
          ret.m_quality = MIN (MIN (m_quality, other.m_quality),
			       GUESSED);
	  return ret;
	}
      else if (!m_val)
	ret.m_val = 0;
      else
	{
	  gcc_checking_assert (other.m_val);
	  ret.m_val = MIN (RDIV ((uint64_t)m_val * max_probability,
				 other.m_val),
			   max_probability);
	}
      ret.m_quality = MIN (MIN (m_quality, other.m_quality), ADJUSTED);
      return ret;
    }

  profile_probability &operator/= (const profile_probability &other)
    {
      if (*this == never ())
	return *this = never ();
      if (!initialized_p () || !other.initialized_p ())
	return *this = uninitialized ();
      else
	{
          /* If we get probability above 1, mark it as unreliable
	     and return 1. */
	  if (m_val > other.m_val)
	    {
	      m_val = max_probability;
              m_quality = MIN (MIN (m_quality, other.m_quality),
			       GUESSED);
	      return *this;
	    }
	  else if (!m_val)
	    ;
	  else
	    {
	      gcc_checking_assert (other.m_val);
	      m_val = MIN (RDIV ((uint64_t)m_val * max_probability,
				 other.m_val),
			   max_probability);
	    }
	  m_quality = MIN (MIN (m_quality, other.m_quality), ADJUSTED);
	}
      return *this;
    }

  /* Split *THIS (ORIG) probability into 2 probabilities, such that
     the returned one (FIRST) is *THIS * CPROB and *THIS is
     adjusted (SECOND) so that FIRST + FIRST.invert () * SECOND
     == ORIG.  This is useful e.g. when splitting a conditional
     branch like:
     if (cond)
       goto lab; // ORIG probability
     into
     if (cond1)
       goto lab; // FIRST = ORIG * CPROB probability
     if (cond2)
       goto lab; // SECOND probability
     such that the overall probability of jumping to lab remains
     the same.  CPROB gives the relative probability between the
     branches.  */
  profile_probability split (const profile_probability &cprob)
    {
      profile_probability ret = *this * cprob;
      /* The following is equivalent to:
         *this = cprob.invert () * *this / ret.invert ();
	 Avoid scaling when overall outcome is supposed to be always.
	 Without knowing that one is inverse of other, the result would be
	 conservative.  */
      if (!(*this == always ()))
        *this = (*this - ret) / ret.invert ();
      return ret;
    }

  gcov_type apply (gcov_type val) const
    {
      if (*this == uninitialized ())
	return val / 2;
      return RDIV (val * m_val, max_probability);
    }

  /* Return 1-*THIS.  */
  profile_probability invert () const
    {
      return always() - *this;
    }

  /* Return THIS with quality dropped to GUESSED.  */
  profile_probability guessed () const
    {
      profile_probability ret = *this;
      ret.m_quality = GUESSED;
      return ret;
    }

  /* Return THIS with quality dropped to AFDO.  */
  profile_probability afdo () const
    {
      profile_probability ret = *this;
      ret.m_quality = AFDO;
      return ret;
    }

  /* Return *THIS * NUM / DEN.  */
  profile_probability apply_scale (int64_t num, int64_t den) const
    {
      if (*this == never ())
	return *this;
      if (!initialized_p ())
	return uninitialized ();
      profile_probability ret;
      uint64_t tmp;
      safe_scale_64bit (m_val, num, den, &tmp);
      ret.m_val = MIN (tmp, max_probability);
      ret.m_quality = MIN (m_quality, ADJUSTED);
      return ret;
    }

  /* Return true when the probability of edge is reliable.

     The profile guessing code is good at predicting branch outcome (i.e.
     taken/not taken), that is predicted right slightly over 75% of time.
     It is however notoriously poor on predicting the probability itself.
     In general the profile appear a lot flatter (with probabilities closer
     to 50%) than the reality so it is bad idea to use it to drive optimization
     such as those disabling dynamic branch prediction for well predictable
     branches.

     There are two exceptions - edges leading to noreturn edges and edges
     predicted by number of iterations heuristics are predicted well.  This macro
     should be able to distinguish those, but at the moment it simply check for
     noreturn heuristic that is only one giving probability over 99% or bellow
     1%.  In future we might want to propagate reliability information across the
     CFG if we find this information useful on multiple places.   */
  bool probably_reliable_p () const
    {
      if (m_quality >= ADJUSTED)
	return true;
      if (!initialized_p ())
	return false;
      return m_val < max_probability / 100
	     || m_val > max_probability - max_probability / 100;
    }

  /* Return false if profile_probability is bogus.  */
  bool verify () const
    {
      gcc_checking_assert (m_quality != UNINITIALIZED_PROFILE);
      if (m_val == uninitialized_probability)
	return m_quality == GUESSED;
      else if (m_quality < GUESSED)
	return false;
      return m_val <= max_probability;
    }

  /* Comparisons are three-state and conservative.  False is returned if
     the inequality cannot be decided.  */
  bool operator< (const profile_probability &other) const
    {
      return initialized_p () && other.initialized_p () && m_val < other.m_val;
    }

  bool operator> (const profile_probability &other) const
    {
      return initialized_p () && other.initialized_p () && m_val > other.m_val;
    }

  bool operator<= (const profile_probability &other) const
    {
      return initialized_p () && other.initialized_p () && m_val <= other.m_val;
    }

  bool operator>= (const profile_probability &other) const
    {
      return initialized_p () && other.initialized_p () && m_val >= other.m_val;
    }

  /* Get the value of the count.  */
  uint32_t value () const { return m_val; }

  /* Get the quality of the count.  */
  enum profile_quality quality () const { return m_quality; }

  /* Output THIS to F.  */
  void dump (FILE *f) const;

  /* Output THIS to BUFFER.  */
  void dump (char *buffer) const;

  /* Print THIS to stderr.  */
  void debug () const;

  /* Return true if THIS is known to differ significantly from OTHER.  */
  bool differs_from_p (profile_probability other) const;

  /* Return if difference is greater than 50%.  */
  bool differs_lot_from_p (profile_probability other) const;

  /* COUNT1 times event happens with *THIS probability, COUNT2 times OTHER
     happens with COUNT2 probability. Return probability that either *THIS or
     OTHER happens.  */
  profile_probability combine_with_count (profile_count count1,
					  profile_probability other,
					  profile_count count2) const;

  /* Return probability as sreal.  */
  sreal to_sreal () const;
  /* LTO streaming support.  */
  static profile_probability stream_in (class lto_input_block *);
  void stream_out (struct output_block *);
  void stream_out (struct lto_output_stream *);
};

/* Main data type to hold profile counters in GCC. Profile counts originate
   either from profile feedback, static profile estimation or both.  We do not
   perform whole program profile propagation and thus profile estimation
   counters are often local to function, while counters from profile feedback
   (or special cases of profile estimation) can be used inter-procedurally.

   There are 3 basic types
     1) local counters which are result of intra-procedural static profile
        estimation.
     2) ipa counters which are result of profile feedback or special case
        of static profile estimation (such as in function main).
     3) counters which counts as 0 inter-procedurally (because given function
        was never run in train feedback) but they hold local static profile
        estimate.

   Counters of type 1 and 3 cannot be mixed with counters of different type
   within operation (because whole function should use one type of counter)
   with exception that global zero mix in most operations where outcome is
   well defined.

   To take local counter and use it inter-procedurally use ipa member function
   which strips information irrelevant at the inter-procedural level.

   Counters are 61bit integers representing number of executions during the
   train run or normalized frequency within the function.

   As the profile is maintained during the compilation, many adjustments are
   made.  Not all transformations can be made precisely, most importantly
   when code is being duplicated.  It also may happen that part of CFG has
   profile counts known while other do not - for example when LTO optimizing
   partly profiled program or when profile was lost due to COMDAT merging.

   For this reason profile_count tracks more information than
   just unsigned integer and it is also ready for profile mismatches.
   The API of this data type represent operations that are natural
   on profile counts - sum, difference and operation with scales and
   probabilities.  All operations are safe by never getting negative counts
   and they do end up in uninitialized scale if any of the parameters is
   uninitialized.

   All comparisons that are three state and handling of probabilities.  Thus
   a < b is not equal to !(a >= b).

   The following pre-defined counts are available:

   profile_count::zero ()  for code that is known to execute zero times at
      runtime (this can be detected statically i.e. for paths leading to
      abort ();
   profile_count::one () for code that is known to execute once (such as
      main () function
   profile_count::uninitialized ()  for unknown execution count.

 */

struct GTY(()) profile_count
{
public:
  /* Use 62bit to hold basic block counters.  Should be at least
     64bit.  Although a counter cannot be negative, we use a signed
     type to hold various extra stages.  */

  static const int n_bits = 61;
  static const uint64_t max_count = ((uint64_t) 1 << n_bits) - 2;
private:
  static const uint64_t uninitialized_count = ((uint64_t) 1 << n_bits) - 1;

#if defined (__arm__) && (__GNUC__ >= 6 && __GNUC__ <= 8)
  /* Work-around for PR88469.  A bug in the gcc-6/7/8 PCS layout code
     incorrectly detects the alignment of a structure where the only
     64-bit aligned object is a bit-field.  We force the alignment of
     the entire field to mitigate this.  */
#define UINT64_BIT_FIELD_ALIGN __attribute__ ((aligned(8)))
#else
#define UINT64_BIT_FIELD_ALIGN
#endif
  uint64_t UINT64_BIT_FIELD_ALIGN m_val : n_bits;
#undef UINT64_BIT_FIELD_ALIGN
  enum profile_quality m_quality : 3;
public:

  /* Return true if both values can meaningfully appear in single function
     body.  We have either all counters in function local or global, otherwise
     operations between them are not really defined well.  */
  bool compatible_p (const profile_count other) const
    {
      if (!initialized_p () || !other.initialized_p ())
	return true;
      if (*this == zero ()
	  || other == zero ())
	return true;
      /* Do not allow nonzero global profile together with local guesses
	 that are globally0.  */
      if (ipa ().nonzero_p ()
	  && !(other.ipa () == other))
	return false;
      if (other.ipa ().nonzero_p ()
	  && !(ipa () == *this))
	return false;
	
      return ipa_p () == other.ipa_p ();
    }

  /* Used for counters which are expected to be never executed.  */
  static profile_count zero ()
    {
      return from_gcov_type (0);
    }

  static profile_count adjusted_zero ()
    {
      profile_count c;
      c.m_val = 0;
      c.m_quality = ADJUSTED;
      return c;
    }

  static profile_count guessed_zero ()
    {
      profile_count c;
      c.m_val = 0;
      c.m_quality = GUESSED;
      return c;
    }

  static profile_count one ()
    {
      return from_gcov_type (1);
    }

  /* Value of counters which has not been initialized. Either because
     initialization did not happen yet or because profile is unknown.  */
  static profile_count uninitialized ()
    {
      profile_count c;
      c.m_val = uninitialized_count;
      c.m_quality = GUESSED_LOCAL;
      return c;
    }

  /* Conversion to gcov_type is lossy.  */
  gcov_type to_gcov_type () const
    {
      gcc_checking_assert (initialized_p ());
      return m_val;
    }

  /* Return true if value has been initialized.  */
  bool initialized_p () const
    {
      return m_val != uninitialized_count;
    }

  /* Return true if value can be trusted.  */
  bool reliable_p () const
    {
      return m_quality >= ADJUSTED;
    }

  /* Return true if value can be operated inter-procedurally.  */
  bool ipa_p () const
    {
      return !initialized_p () || m_quality >= GUESSED_GLOBAL0;
    }

  /* Return true if quality of profile is precise.  */
  bool precise_p () const
    {
      return m_quality == PRECISE;
    }

  /* Get the value of the count.  */
  uint64_t value () const { return m_val; }

  /* Get the quality of the count.  */
  enum profile_quality quality () const { return m_quality; }

  /* When merging basic blocks, the two different profile counts are unified.
     Return true if this can be done without losing info about profile.
     The only case we care about here is when first BB contains something
     that makes it terminate in a way not visible in CFG.  */
  bool ok_for_merging (profile_count other) const
    {
      if (m_quality < ADJUSTED
	  || other.m_quality < ADJUSTED)
	return true;
      return !(other < *this);
    }

  /* When merging two BBs with different counts, pick common count that looks
     most representative.  */
  profile_count merge (profile_count other) const
    {
      if (*this == other || !other.initialized_p ()
	  || m_quality > other.m_quality)
	return *this;
      if (other.m_quality > m_quality
	  || other > *this)
	return other;
      return *this;
    }

  /* Basic operations.  */
  bool operator== (const profile_count &other) const
    {
      return m_val == other.m_val && m_quality == other.m_quality;
    }

  profile_count operator+ (const profile_count &other) const
    {
      if (other == zero ())
	return *this;
      if (*this == zero ())
	return other;
      if (!initialized_p () || !other.initialized_p ())
	return uninitialized ();

      profile_count ret;
      gcc_checking_assert (compatible_p (other));
      ret.m_val = m_val + other.m_val;
      ret.m_quality = MIN (m_quality, other.m_quality);
      return ret;
    }

  profile_count &operator+= (const profile_count &other)
    {
      if (other == zero ())
	return *this;
      if (*this == zero ())
	{
	  *this = other;
	  return *this;
	}
      if (!initialized_p () || !other.initialized_p ())
	return *this = uninitialized ();
      else
	{
          gcc_checking_assert (compatible_p (other));
	  m_val += other.m_val;
	  m_quality = MIN (m_quality, other.m_quality);
	}
      return *this;
    }

  profile_count operator- (const profile_count &other) const
    {
      if (*this == zero () || other == zero ())
	return *this;
      if (!initialized_p () || !other.initialized_p ())
	return uninitialized ();
      gcc_checking_assert (compatible_p (other));
      profile_count ret;
      ret.m_val = m_val >= other.m_val ? m_val - other.m_val : 0;
      ret.m_quality = MIN (m_quality, other.m_quality);
      return ret;
    }

  profile_count &operator-= (const profile_count &other)
    {
      if (*this == zero () || other == zero ())
	return *this;
      if (!initialized_p () || !other.initialized_p ())
	return *this = uninitialized ();
      else
	{
          gcc_checking_assert (compatible_p (other));
	  m_val = m_val >= other.m_val ? m_val - other.m_val: 0;
	  m_quality = MIN (m_quality, other.m_quality);
	}
      return *this;
    }

  /* Return false if profile_count is bogus.  */
  bool verify () const
    {
      gcc_checking_assert (m_quality != UNINITIALIZED_PROFILE);
      return m_val != uninitialized_count || m_quality == GUESSED_LOCAL;
    }

  /* Comparisons are three-state and conservative.  False is returned if
     the inequality cannot be decided.  */
  bool operator< (const profile_count &other) const
    {
      if (!initialized_p () || !other.initialized_p ())
	return false;
      if (*this == zero ())
	return !(other == zero ());
      if (other == zero ())
	return false;
      gcc_checking_assert (compatible_p (other));
      return m_val < other.m_val;
    }

  bool operator> (const profile_count &other) const
    {
      if (!initialized_p () || !other.initialized_p ())
	return false;
      if (*this  == zero ())
	return false;
      if (other == zero ())
	return !(*this == zero ());
      gcc_checking_assert (compatible_p (other));
      return initialized_p () && other.initialized_p () && m_val > other.m_val;
    }

  bool operator< (const gcov_type other) const
    {
      gcc_checking_assert (ipa_p ());
      gcc_checking_assert (other >= 0);
      return ipa ().initialized_p () && ipa ().m_val < (uint64_t) other;
    }

  bool operator> (const gcov_type other) const
    {
      gcc_checking_assert (ipa_p ());
      gcc_checking_assert (other >= 0);
      return ipa ().initialized_p () && ipa ().m_val > (uint64_t) other;
    }

  bool operator<= (const profile_count &other) const
    {
      if (!initialized_p () || !other.initialized_p ())
	return false;
      if (*this == zero ())
	return true;
      if (other == zero ())
	return (*this == zero ());
      gcc_checking_assert (compatible_p (other));
      return m_val <= other.m_val;
    }

  bool operator>= (const profile_count &other) const
    {
      if (!initialized_p () || !other.initialized_p ())
	return false;
      if (other == zero ())
	return true;
      if (*this == zero ())
	return (other == zero ());
      gcc_checking_assert (compatible_p (other));
      return m_val >= other.m_val;
    }

  bool operator<= (const gcov_type other) const
    {
      gcc_checking_assert (ipa_p ());
      gcc_checking_assert (other >= 0);
      return ipa ().initialized_p () && ipa ().m_val <= (uint64_t) other;
    }

  bool operator>= (const gcov_type other) const
    {
      gcc_checking_assert (ipa_p ());
      gcc_checking_assert (other >= 0);
      return ipa ().initialized_p () && ipa ().m_val >= (uint64_t) other;
    }

  /* Return true when value is not zero and can be used for scaling. 
     This is different from *this > 0 because that requires counter to
     be IPA.  */
  bool nonzero_p () const
    {
      return initialized_p () && m_val != 0;
    }

  /* Make counter forcibly nonzero.  */
  profile_count force_nonzero () const
    {
      if (!initialized_p ())
	return *this;
      profile_count ret = *this;
      if (ret.m_val == 0)
	{
	  ret.m_val = 1;
          ret.m_quality = MIN (m_quality, ADJUSTED);
	}
      return ret;
    }

  profile_count max (profile_count other) const
    {
      profile_count val = *this;

      /* Always prefer nonzero IPA counts over local counts.  */
      if (ipa ().nonzero_p () || other.ipa ().nonzero_p ())
	{
	  val = ipa ();
	  other = other.ipa ();
	}
      if (!initialized_p ())
	return other;
      if (!other.initialized_p ())
	return *this;
      if (*this == zero ())
	return other;
      if (other == zero ())
	return *this;
      gcc_checking_assert (compatible_p (other));
      if (val.m_val < other.m_val || (m_val == other.m_val
				      && val.m_quality < other.m_quality))
	return other;
      return *this;
    }

  /* PROB is a probability in scale 0...REG_BR_PROB_BASE.  Scale counter
     accordingly.  */
  profile_count apply_probability (int prob) const
    {
      gcc_checking_assert (prob >= 0 && prob <= REG_BR_PROB_BASE);
      if (m_val == 0)
	return *this;
      if (!initialized_p ())
	return uninitialized ();
      profile_count ret;
      ret.m_val = RDIV (m_val * prob, REG_BR_PROB_BASE);
      ret.m_quality = MIN (m_quality, ADJUSTED);
      return ret;
    }

  /* Scale counter according to PROB.  */
  profile_count apply_probability (profile_probability prob) const
    {
      if (*this == zero ())
	return *this;
      if (prob == profile_probability::never ())
	return zero ();
      if (!initialized_p ())
	return uninitialized ();
      profile_count ret;
      uint64_t tmp;
      safe_scale_64bit (m_val, prob.m_val, profile_probability::max_probability,
			&tmp);
      ret.m_val = tmp;
      ret.m_quality = MIN (m_quality, prob.m_quality);
      return ret;
    }

  /* Return *THIS * NUM / DEN.  */
  profile_count apply_scale (int64_t num, int64_t den) const
    {
      if (m_val == 0)
	return *this;
      if (!initialized_p ())
	return uninitialized ();
      profile_count ret;
      uint64_t tmp;

      gcc_checking_assert (num >= 0 && den > 0);
      safe_scale_64bit (m_val, num, den, &tmp);
      ret.m_val = MIN (tmp, max_count);
      ret.m_quality = MIN (m_quality, ADJUSTED);
      return ret;
    }

  profile_count apply_scale (profile_count num, profile_count den) const
    {
      if (*this == zero ())
	return *this;
      if (num == zero ())
	return num;
      if (!initialized_p () || !num.initialized_p () || !den.initialized_p ())
	return uninitialized ();
      if (num == den)
	return *this;
      gcc_checking_assert (den.m_val);

      profile_count ret;
      uint64_t val;
      safe_scale_64bit (m_val, num.m_val, den.m_val, &val);
      ret.m_val = MIN (val, max_count);
      ret.m_quality = MIN (MIN (MIN (m_quality, ADJUSTED),
			        num.m_quality), den.m_quality);
      /* Be sure that ret is not local if num is global.
	 Also ensure that ret is not global0 when num is global.  */
      if (num.ipa_p ())
	ret.m_quality = MAX (ret.m_quality,
			     num == num.ipa () ? GUESSED : num.m_quality);
      return ret;
    }

  /* Return THIS with quality dropped to GUESSED_LOCAL.  */
  profile_count guessed_local () const
    {
      profile_count ret = *this;
      if (!initialized_p ())
	return *this;
      ret.m_quality = GUESSED_LOCAL;
      return ret;
    }

  /* We know that profile is globally 0 but keep local profile if present.  */
  profile_count global0 () const
    {
      profile_count ret = *this;
      if (!initialized_p ())
	return *this;
      ret.m_quality = GUESSED_GLOBAL0;
      return ret;
    }

  /* We know that profile is globally adjusted 0 but keep local profile
     if present.  */
  profile_count global0adjusted () const
    {
      profile_count ret = *this;
      if (!initialized_p ())
	return *this;
      ret.m_quality = GUESSED_GLOBAL0_ADJUSTED;
      return ret;
    }

  /* Return THIS with quality dropped to GUESSED.  */
  profile_count guessed () const
    {
      profile_count ret = *this;
      ret.m_quality = MIN (ret.m_quality, GUESSED);
      return ret;
    }

  /* Return variant of profile count which is always safe to compare
     across functions.  */
  profile_count ipa () const
    {
      if (m_quality > GUESSED_GLOBAL0_ADJUSTED)
	return *this;
      if (m_quality == GUESSED_GLOBAL0)
	return zero ();
      if (m_quality == GUESSED_GLOBAL0_ADJUSTED)
	return adjusted_zero ();
      return uninitialized ();
    }

  /* Return THIS with quality dropped to AFDO.  */
  profile_count afdo () const
    {
      profile_count ret = *this;
      ret.m_quality = AFDO;
      return ret;
    }

  /* Return probability of event with counter THIS within event with counter
     OVERALL.  */
  profile_probability probability_in (const profile_count overall) const
    {
      if (*this == zero ()
	  && !(overall == zero ()))
	return profile_probability::never ();
      if (!initialized_p () || !overall.initialized_p ()
	  || !overall.m_val)
	return profile_probability::uninitialized ();
      if (*this == overall && m_quality == PRECISE)
	return profile_probability::always ();
      profile_probability ret;
      gcc_checking_assert (compatible_p (overall));

      if (overall.m_val < m_val)
	{
	  ret.m_val = profile_probability::max_probability;
	  ret.m_quality = GUESSED;
	  return ret;
	}
      else
	ret.m_val = RDIV (m_val * profile_probability::max_probability,
			  overall.m_val);
      ret.m_quality = MIN (MAX (MIN (m_quality, overall.m_quality),
				GUESSED), ADJUSTED);
      return ret;
    }

  int to_frequency (struct function *fun) const;
  int to_cgraph_frequency (profile_count entry_bb_count) const;
  sreal to_sreal_scale (profile_count in, bool *known = NULL) const;

  /* Output THIS to F.  */
  void dump (FILE *f) const;

  /* Output THIS to BUFFER.  */
  void dump (char *buffer) const;

  /* Print THIS to stderr.  */
  void debug () const;

  /* Return true if THIS is known to differ significantly from OTHER.  */
  bool differs_from_p (profile_count other) const;

  /* We want to scale profile across function boundary from NUM to DEN.
     Take care of the side case when NUM and DEN are zeros of incompatible
     kinds.  */
  static void adjust_for_ipa_scaling (profile_count *num, profile_count *den);

  /* THIS is a count of bb which is known to be executed IPA times.
     Combine this information into bb counter.  This means returning IPA
     if it is nonzero, not changing anything if IPA is uninitialized
     and if IPA is zero, turning THIS into corresponding local profile with
     global0.  */
  profile_count combine_with_ipa_count (profile_count ipa);

  /* Same as combine_with_ipa_count but inside function with count IPA2.  */
  profile_count combine_with_ipa_count_within
		 (profile_count ipa, profile_count ipa2);

  /* The profiling runtime uses gcov_type, which is usually 64bit integer.
     Conversions back and forth are used to read the coverage and get it
     into internal representation.  */
  static profile_count from_gcov_type (gcov_type v,
				       profile_quality quality = PRECISE);

  /* LTO streaming support.  */
  static profile_count stream_in (class lto_input_block *);
  void stream_out (struct output_block *);
  void stream_out (struct lto_output_stream *);
};
#endif
