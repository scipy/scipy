/* Generated automatically by the program `genattr'
   from the machine description file `md'.  */

#ifndef GCC_INSN_ATTR_H
#define GCC_INSN_ATTR_H

#include "insn-attr-common.h"

#define HAVE_ATTR_cpu 1
extern enum attr_cpu get_attr_cpu (void);

#define HAVE_ATTR_type 1
extern enum attr_type get_attr_type (rtx_insn *);

#define HAVE_ATTR_mode 1
extern enum attr_mode get_attr_mode (rtx_insn *);

#define HAVE_ATTR_unit 1
extern enum attr_unit get_attr_unit (rtx_insn *);

#define HAVE_ATTR_length_immediate 1
extern int get_attr_length_immediate (rtx_insn *);
#define HAVE_ATTR_length_address 1
extern int get_attr_length_address (rtx_insn *);
#define HAVE_ATTR_prefix_data16 1
extern int get_attr_prefix_data16 (rtx_insn *);
#define HAVE_ATTR_prefix_rep 1
extern int get_attr_prefix_rep (rtx_insn *);
#define HAVE_ATTR_prefix_0f 1
extern int get_attr_prefix_0f (rtx_insn *);
#define HAVE_ATTR_prefix_rex 1
extern int get_attr_prefix_rex (rtx_insn *);
#define HAVE_ATTR_prefix_extra 1
extern int get_attr_prefix_extra (rtx_insn *);
#define HAVE_ATTR_prefix 1
extern enum attr_prefix get_attr_prefix (rtx_insn *);

#define HAVE_ATTR_prefix_vex_w 1
extern int get_attr_prefix_vex_w (rtx_insn *);
#define HAVE_ATTR_length_vex 1
extern int get_attr_length_vex (rtx_insn *);
#define HAVE_ATTR_length_evex 1
extern int get_attr_length_evex (rtx_insn *);
#define HAVE_ATTR_modrm 1
extern int get_attr_modrm (rtx_insn *);
#define HAVE_ATTR_length 1
extern int get_attr_length (rtx_insn *);
extern void shorten_branches (rtx_insn *);
extern int insn_default_length (rtx_insn *);
extern int insn_min_length (rtx_insn *);
extern int insn_variable_length_p (rtx_insn *);
extern int insn_current_length (rtx_insn *);

#include "insn-addr.h"

#define HAVE_ATTR_memory 1
extern enum attr_memory get_attr_memory (rtx_insn *);

#define HAVE_ATTR_imm_disp 1
extern enum attr_imm_disp get_attr_imm_disp (rtx_insn *);

#define HAVE_ATTR_fp_int_src 1
extern enum attr_fp_int_src get_attr_fp_int_src (rtx_insn *);

#define HAVE_ATTR_i387_cw 1
extern enum attr_i387_cw get_attr_i387_cw (rtx_insn *);

#define HAVE_ATTR_avx_partial_xmm_update 1
extern enum attr_avx_partial_xmm_update get_attr_avx_partial_xmm_update (rtx_insn *);

#define HAVE_ATTR_use_carry 1
extern enum attr_use_carry get_attr_use_carry (rtx_insn *);

#define HAVE_ATTR_movu 1
extern enum attr_movu get_attr_movu (rtx_insn *);

#define HAVE_ATTR_isa 1
extern enum attr_isa get_attr_isa (rtx_insn *);

#define HAVE_ATTR_mmx_isa 1
extern enum attr_mmx_isa get_attr_mmx_isa (rtx_insn *);

#define HAVE_ATTR_enabled 1
extern int get_attr_enabled (rtx_insn *);
#define HAVE_ATTR_preferred_for_size 1
extern int get_attr_preferred_for_size (rtx_insn *);
#define HAVE_ATTR_preferred_for_speed 1
extern int get_attr_preferred_for_speed (rtx_insn *);
#define HAVE_ATTR_pent_prefix 1
extern enum attr_pent_prefix get_attr_pent_prefix (rtx_insn *);

#define HAVE_ATTR_pent_pair 1
extern enum attr_pent_pair get_attr_pent_pair (rtx_insn *);

#define HAVE_ATTR_athlon_decode 1
extern enum attr_athlon_decode get_attr_athlon_decode (rtx_insn *);

#define HAVE_ATTR_amdfam10_decode 1
extern enum attr_amdfam10_decode get_attr_amdfam10_decode (rtx_insn *);

#define HAVE_ATTR_bdver1_decode 1
extern enum attr_bdver1_decode get_attr_bdver1_decode (rtx_insn *);

#define HAVE_ATTR_btver2_decode 1
extern enum attr_btver2_decode get_attr_btver2_decode (rtx_insn *);

#define HAVE_ATTR_btver2_sse_attr 1
extern enum attr_btver2_sse_attr get_attr_btver2_sse_attr (rtx_insn *);

#define HAVE_ATTR_znver1_decode 1
extern enum attr_znver1_decode get_attr_znver1_decode (rtx_insn *);

#define HAVE_ATTR_znver4_decode 1
extern enum attr_znver4_decode get_attr_znver4_decode (rtx_insn *);

#define HAVE_ATTR_atom_unit 1
extern enum attr_atom_unit get_attr_atom_unit (rtx_insn *);

#define HAVE_ATTR_atom_sse_attr 1
extern enum attr_atom_sse_attr get_attr_atom_sse_attr (rtx_insn *);

#define HAVE_ATTR_i7_domain 1
extern enum attr_i7_domain get_attr_i7_domain (rtx_insn *);

#define HAVE_ATTR_hsw_domain 1
extern enum attr_hsw_domain get_attr_hsw_domain (rtx_insn *);

#define HAVE_ATTR_unord_subst 1
extern enum attr_unord_subst get_attr_unord_subst (rtx_insn *);

#define HAVE_ATTR_mask 1
extern enum attr_mask get_attr_mask (rtx_insn *);

#define HAVE_ATTR_maskc 1
extern enum attr_maskc get_attr_maskc (rtx_insn *);

#define HAVE_ATTR_merge_mask 1
extern enum attr_merge_mask get_attr_merge_mask (rtx_insn *);

#define HAVE_ATTR_mask_scalar_merge 1
extern enum attr_mask_scalar_merge get_attr_mask_scalar_merge (rtx_insn *);

#define HAVE_ATTR_sd 1
extern enum attr_sd get_attr_sd (rtx_insn *);

#define HAVE_ATTR_maskz_scalar 1
extern enum attr_maskz_scalar get_attr_maskz_scalar (rtx_insn *);

#define HAVE_ATTR_sdc 1
extern enum attr_sdc get_attr_sdc (rtx_insn *);

#define HAVE_ATTR_round 1
extern enum attr_round get_attr_round (rtx_insn *);

#define HAVE_ATTR_round_saeonly 1
extern enum attr_round_saeonly get_attr_round_saeonly (rtx_insn *);

#define HAVE_ATTR_round_expand 1
extern enum attr_round_expand get_attr_round_expand (rtx_insn *);

#define HAVE_ATTR_round_saeonly_expand 1
extern enum attr_round_saeonly_expand get_attr_round_saeonly_expand (rtx_insn *);

#define HAVE_ATTR_mask_expand4 1
extern enum attr_mask_expand4 get_attr_mask_expand4 (rtx_insn *);

#define HAVE_ATTR_mask_scalar 1
extern enum attr_mask_scalar get_attr_mask_scalar (rtx_insn *);

#define HAVE_ATTR_mask_scalarcz 1
extern enum attr_mask_scalarcz get_attr_mask_scalarcz (rtx_insn *);

#define HAVE_ATTR_mask_scalarc 1
extern enum attr_mask_scalarc get_attr_mask_scalarc (rtx_insn *);

#define HAVE_ATTR_round_scalar 1
extern enum attr_round_scalar get_attr_round_scalar (rtx_insn *);

#define HAVE_ATTR_round_scalarcz 1
extern enum attr_round_scalarcz get_attr_round_scalarcz (rtx_insn *);

#define HAVE_ATTR_round_saeonly_scalar 1
extern enum attr_round_saeonly_scalar get_attr_round_saeonly_scalar (rtx_insn *);

#define HAVE_ATTR_maskz_half 1
extern enum attr_maskz_half get_attr_maskz_half (rtx_insn *);

extern int num_delay_slots (rtx_insn *);
extern int eligible_for_delay (rtx_insn *, int, rtx_insn *, int);

extern int const_num_delay_slots (rtx_insn *);

#define ANNUL_IFTRUE_SLOTS 0
extern int eligible_for_annul_true (rtx_insn *, int, rtx_insn *, int);
#define ANNUL_IFFALSE_SLOTS 0
extern int eligible_for_annul_false (rtx_insn *, int, rtx_insn *, int);

/* DFA based pipeline interface.  */
#ifndef AUTOMATON_ALTS
#define AUTOMATON_ALTS 0
#endif


#ifndef AUTOMATON_STATE_ALTS
#define AUTOMATON_STATE_ALTS 0
#endif

#ifndef CPU_UNITS_QUERY
#define CPU_UNITS_QUERY 0
#endif

/* Initialize fn pointers for internal_dfa_insn_code
   and insn_default_latency.  */
extern void init_sched_attrs (void);

/* Internal insn code number used by automata.  */
extern int (*internal_dfa_insn_code) (rtx_insn *);

/* Insn latency time defined in define_insn_reservation. */
extern int (*insn_default_latency) (rtx_insn *);

/* Return nonzero if there is a bypass for given insn
   which is a data producer.  */
extern int bypass_p (rtx_insn *);

/* Insn latency time on data consumed by the 2nd insn.
   Use the function if bypass_p returns nonzero for
   the 1st insn. */
extern int insn_latency (rtx_insn *, rtx_insn *);

/* Maximal insn latency time possible of all bypasses for this insn.
   Use the function if bypass_p returns nonzero for
   the 1st insn. */
extern int maximal_insn_latency (rtx_insn *);


#if AUTOMATON_ALTS
/* The following function returns number of alternative
   reservations of given insn.  It may be used for better
   insns scheduling heuristics. */
extern int insn_alts (rtx);

#endif

/* Maximal possible number of insns waiting results being
   produced by insns whose execution is not finished. */
extern const int max_insn_queue_index;

/* Pointer to data describing current state of DFA.  */
typedef void *state_t;

/* Size of the data in bytes.  */
extern int state_size (void);

/* Initiate given DFA state, i.e. Set up the state
   as all functional units were not reserved.  */
extern void state_reset (state_t);
/* The following function returns negative value if given
   insn can be issued in processor state described by given
   DFA state.  In this case, the DFA state is changed to
   reflect the current and future reservations by given
   insn.  Otherwise the function returns minimal time
   delay to issue the insn.  This delay may be zero
   for superscalar or VLIW processors.  If the second
   parameter is NULL the function changes given DFA state
   as new processor cycle started.  */
extern int state_transition (state_t, rtx);

#if AUTOMATON_STATE_ALTS
/* The following function returns number of possible
   alternative reservations of given insn in given
   DFA state.  It may be used for better insns scheduling
   heuristics.  By default the function is defined if
   macro AUTOMATON_STATE_ALTS is defined because its
   implementation may require much memory.  */
extern int state_alts (state_t, rtx);
#endif

extern int min_issue_delay (state_t, rtx_insn *);
/* The following function returns nonzero if no one insn
   can be issued in current DFA state. */
extern int state_dead_lock_p (state_t);
/* The function returns minimal delay of issue of the 2nd
   insn after issuing the 1st insn in given DFA state.
   The 1st insn should be issued in given state (i.e.
    state_transition should return negative value for
    the insn and the state).  Data dependencies between
    the insns are ignored by the function.  */
extern int min_insn_conflict_delay (state_t, rtx_insn *, rtx_insn *);
/* The following function outputs reservations for given
   insn as they are described in the corresponding
   define_insn_reservation.  */
extern void print_reservation (FILE *, rtx_insn *);

#if CPU_UNITS_QUERY
/* The following function returns code of functional unit
   with given name (see define_cpu_unit). */
extern int get_cpu_unit_code (const char *);
/* The following function returns nonzero if functional
   unit with given code is currently reserved in given
   DFA state.  */
extern int cpu_unit_reservation_p (state_t, int);
#endif

/* The following function returns true if insn
   has a dfa reservation.  */
extern bool insn_has_dfa_reservation_p (rtx_insn *);

/* Clean insn code cache.  It should be called if there
   is a chance that condition value in a
   define_insn_reservation will be changed after
   last call of dfa_start.  */
extern void dfa_clean_insn_cache (void);

extern void dfa_clear_single_insn_cache (rtx_insn *);

/* Initiate and finish work with DFA.  They should be
   called as the first and the last interface
   functions.  */
extern void dfa_start (void);
extern void dfa_finish (void);
#ifndef HAVE_ATTR_length
#define HAVE_ATTR_length 0
#endif
#ifndef HAVE_ATTR_enabled
#define HAVE_ATTR_enabled 0
#endif
#ifndef HAVE_ATTR_preferred_for_size
#define HAVE_ATTR_preferred_for_size 0
#endif
#ifndef HAVE_ATTR_preferred_for_speed
#define HAVE_ATTR_preferred_for_speed 0
#endif
#if !HAVE_ATTR_length
extern int hook_int_rtx_insn_unreachable (rtx_insn *);
#define insn_default_length hook_int_rtx_insn_unreachable
#define insn_min_length hook_int_rtx_insn_unreachable
#define insn_variable_length_p hook_int_rtx_insn_unreachable
#define insn_current_length hook_int_rtx_insn_unreachable
#include "insn-addr.h"
#endif
extern int hook_int_rtx_1 (rtx);
#if !HAVE_ATTR_enabled
#define get_attr_enabled hook_int_rtx_1
#endif
#if !HAVE_ATTR_preferred_for_size
#define get_attr_preferred_for_size hook_int_rtx_1
#endif
#if !HAVE_ATTR_preferred_for_speed
#define get_attr_preferred_for_speed hook_int_rtx_1
#endif


#define ATTR_FLAG_forward	0x1
#define ATTR_FLAG_backward	0x2

#endif /* GCC_INSN_ATTR_H */
