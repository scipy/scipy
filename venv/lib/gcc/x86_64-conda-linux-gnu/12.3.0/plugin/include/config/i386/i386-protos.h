/* Definitions of target machine for GCC for IA-32.
   Copyright (C) 1988-2022 Free Software Foundation, Inc.

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

/* In i386-common.cc.  */
extern bool ix86_handle_option (struct gcc_options *opts,
				struct gcc_options *opts_set ATTRIBUTE_UNUSED,
				const struct cl_decoded_option *decoded,
				location_t loc);

/* Functions in i386.cc */
extern bool ix86_target_stack_probe (void);
extern bool ix86_can_use_return_insn_p (void);
extern bool ix86_function_ms_hook_prologue (const_tree fn);
extern void ix86_setup_frame_addresses (void);
extern bool ix86_rip_relative_addr_p (struct ix86_address *parts);

extern HOST_WIDE_INT ix86_initial_elimination_offset (int, int);
extern void ix86_expand_prologue (void);
extern void ix86_maybe_emit_epilogue_vzeroupper (void);
extern void ix86_expand_epilogue (int);
extern void ix86_expand_split_stack_prologue (void);

extern void ix86_output_addr_vec_elt (FILE *, int);
extern void ix86_output_addr_diff_elt (FILE *, int, int);

extern const char *ix86_output_ssemov (rtx_insn *, rtx *);

extern enum calling_abi ix86_cfun_abi (void);
extern enum calling_abi ix86_function_type_abi (const_tree);

extern bool ix86_use_pseudo_pic_reg (void);

extern void ix86_reset_previous_fndecl (void);

extern bool ix86_using_red_zone (void);

extern rtx ix86_gen_scratch_sse_rtx (machine_mode);

extern unsigned int ix86_regmode_natural_size (machine_mode);
extern bool ix86_check_builtin_isa_match (unsigned int fcode);
#ifdef RTX_CODE
extern int standard_80387_constant_p (rtx);
extern const char *standard_80387_constant_opcode (rtx);
extern rtx standard_80387_constant_rtx (int);
extern int standard_sse_constant_p (rtx, machine_mode);
extern const char *standard_sse_constant_opcode (rtx_insn *, rtx *);
extern bool ix86_standard_x87sse_constant_load_p (const rtx_insn *, rtx);
extern bool ix86_pre_reload_split (void);
extern bool symbolic_reference_mentioned_p (rtx);
extern bool extended_reg_mentioned_p (rtx);
extern bool x86_extended_QIreg_mentioned_p (rtx_insn *);
extern bool x86_extended_reg_mentioned_p (rtx);
extern bool x86_maybe_negate_const_int (rtx *, machine_mode);
extern machine_mode ix86_cc_mode (enum rtx_code, rtx, rtx);

extern int avx_vpermilp_parallel (rtx par, machine_mode mode);
extern int avx_vperm2f128_parallel (rtx par, machine_mode mode);

extern bool ix86_expand_strlen (rtx, rtx, rtx, rtx);
extern bool ix86_expand_set_or_cpymem (rtx, rtx, rtx, rtx, rtx, rtx,
				       rtx, rtx, rtx, rtx, bool);
extern bool ix86_expand_cmpstrn_or_cmpmem (rtx, rtx, rtx, rtx, rtx, bool);

extern bool constant_address_p (rtx);
extern bool legitimate_pic_operand_p (rtx);
extern bool legitimate_pic_address_disp_p (rtx);
extern bool ix86_force_load_from_GOT_p (rtx, bool = false);
extern void print_reg (rtx, int, FILE*);
extern void ix86_print_operand (FILE *, rtx, int);

extern void split_double_mode (machine_mode, rtx[], int, rtx[], rtx[]);

extern const char *output_set_got (rtx, rtx);
extern const char *output_387_binary_op (rtx_insn *, rtx*);
extern const char *output_387_reg_move (rtx_insn *, rtx*);
extern const char *output_fix_trunc (rtx_insn *, rtx*, bool);
extern const char *output_fp_compare (rtx_insn *, rtx*, bool, bool);
extern const char *output_adjust_stack_and_probe (rtx);
extern const char *output_probe_stack_range (rtx, rtx);

extern void ix86_output_patchable_area (unsigned int, bool);

extern void ix86_expand_clear (rtx);
extern void ix86_expand_move (machine_mode, rtx[]);
extern void ix86_expand_vector_move (machine_mode, rtx[]);
extern void ix86_expand_vector_move_misalign (machine_mode, rtx[]);
extern rtx ix86_fixup_binary_operands (enum rtx_code,
				       machine_mode, rtx[]);
extern void ix86_fixup_binary_operands_no_copy (enum rtx_code,
						machine_mode, rtx[]);
extern void ix86_expand_binary_operator (enum rtx_code,
					 machine_mode, rtx[]);
extern void ix86_expand_vector_logical_operator (enum rtx_code,
						 machine_mode, rtx[]);
extern bool ix86_binary_operator_ok (enum rtx_code, machine_mode, rtx[]);
extern bool ix86_avoid_lea_for_add (rtx_insn *, rtx[]);
extern bool ix86_use_lea_for_mov (rtx_insn *, rtx[]);
extern bool ix86_avoid_lea_for_addr (rtx_insn *, rtx[]);
extern void ix86_split_lea_for_addr (rtx_insn *, rtx[], machine_mode);
extern bool ix86_lea_for_add_ok (rtx_insn *, rtx[]);
extern int ix86_last_zero_store_uid;
extern bool ix86_vec_interleave_v2df_operator_ok (rtx operands[3], bool high);
extern bool ix86_dep_by_shift_count (const_rtx set_insn, const_rtx use_insn);
extern bool ix86_agi_dependent (rtx_insn *set_insn, rtx_insn *use_insn);
extern void ix86_expand_unary_operator (enum rtx_code, machine_mode,
					rtx[]);
extern rtx ix86_build_const_vector (machine_mode, bool, rtx);
extern rtx ix86_build_signbit_mask (machine_mode, bool, bool);
extern void ix86_split_convert_uns_si_sse (rtx[]);
extern void ix86_expand_convert_uns_didf_sse (rtx, rtx);
extern void ix86_expand_convert_uns_sixf_sse (rtx, rtx);
extern void ix86_expand_convert_uns_sidf_sse (rtx, rtx);
extern void ix86_expand_convert_uns_sisf_sse (rtx, rtx);
extern void ix86_expand_convert_sign_didf_sse (rtx, rtx);
extern void ix86_expand_vector_convert_uns_vsivsf (rtx, rtx);
extern rtx ix86_expand_adjust_ufix_to_sfix_si (rtx, rtx *);
extern enum ix86_fpcmp_strategy ix86_fp_comparison_strategy (enum rtx_code);
extern void ix86_expand_fp_absneg_operator (enum rtx_code, machine_mode,
					    rtx[]);
extern void ix86_split_fp_absneg_operator (enum rtx_code, machine_mode,
					   rtx[]);
extern void ix86_expand_copysign (rtx []);
extern void ix86_expand_xorsign (rtx []);
extern bool ix86_unary_operator_ok (enum rtx_code, machine_mode, rtx[]);
extern bool ix86_match_ccmode (rtx, machine_mode);
extern void ix86_expand_branch (enum rtx_code, rtx, rtx, rtx);
extern void ix86_expand_setcc (rtx, enum rtx_code, rtx, rtx);
extern bool ix86_expand_int_movcc (rtx[]);
extern bool ix86_expand_fp_movcc (rtx[]);
extern bool ix86_expand_fp_vcond (rtx[]);
extern bool ix86_expand_int_vcond (rtx[]);
extern void ix86_expand_vec_perm (rtx[]);
extern bool ix86_expand_mask_vec_cmp (rtx, enum rtx_code, rtx, rtx);
extern bool ix86_expand_int_vec_cmp (rtx[]);
extern bool ix86_expand_fp_vec_cmp (rtx[]);
extern void ix86_expand_sse_movcc (rtx, rtx, rtx, rtx);
extern void ix86_expand_sse_unpack (rtx, rtx, bool, bool);
extern void ix86_expand_fp_spaceship (rtx, rtx, rtx);
extern bool ix86_expand_int_addcc (rtx[]);
extern rtx_insn *ix86_expand_call (rtx, rtx, rtx, rtx, rtx, bool);
extern bool ix86_call_use_plt_p (rtx);
extern void ix86_split_call_vzeroupper (rtx, rtx);
extern void x86_initialize_trampoline (rtx, rtx, rtx);
extern rtx ix86_zero_extend_to_Pmode (rtx);
extern void ix86_split_long_move (rtx[]);
extern void ix86_split_ashl (rtx *, rtx, machine_mode);
extern void ix86_split_ashr (rtx *, rtx, machine_mode);
extern void ix86_split_lshr (rtx *, rtx, machine_mode);
extern void ix86_expand_v1ti_shift (enum rtx_code, rtx[]);
extern void ix86_expand_v1ti_rotate (enum rtx_code, rtx[]);
extern void ix86_expand_v1ti_ashiftrt (rtx[]);
extern rtx ix86_find_base_term (rtx);
extern bool ix86_check_movabs (rtx, int);
extern bool ix86_check_no_addr_space (rtx);
extern void ix86_split_idivmod (machine_mode, rtx[], bool);
extern bool ix86_hardreg_mov_ok (rtx, rtx);

extern rtx assign_386_stack_local (machine_mode, enum ix86_stack_slot);
extern int ix86_attr_length_immediate_default (rtx_insn *, bool);
extern int ix86_attr_length_address_default (rtx_insn *);
extern int ix86_attr_length_vex_default (rtx_insn *, bool, bool);

extern rtx ix86_libcall_value (machine_mode);
extern bool ix86_function_arg_regno_p (int);
extern void ix86_asm_output_function_label (FILE *, const char *, tree);
extern void ix86_call_abi_override (const_tree);
extern int ix86_reg_parm_stack_space (const_tree);

extern bool ix86_libc_has_function (enum function_class fn_class);

extern void x86_order_regs_for_local_alloc (void);
extern void x86_function_profiler (FILE *, int);
extern void x86_emit_floatuns (rtx [2]);
extern void ix86_emit_fp_unordered_jump (rtx);

extern void ix86_emit_i387_sinh (rtx, rtx);
extern void ix86_emit_i387_cosh (rtx, rtx);
extern void ix86_emit_i387_tanh (rtx, rtx);
extern void ix86_emit_i387_asinh (rtx, rtx);
extern void ix86_emit_i387_acosh (rtx, rtx);
extern void ix86_emit_i387_atanh (rtx, rtx);
extern void ix86_emit_i387_log1p (rtx, rtx);
extern void ix86_emit_i387_round (rtx, rtx);
extern void ix86_emit_swdivsf (rtx, rtx, rtx, machine_mode);
extern void ix86_emit_swsqrtsf (rtx, rtx, machine_mode, bool);

extern enum rtx_code ix86_reverse_condition (enum rtx_code, machine_mode);

extern void ix86_expand_lround (rtx, rtx);
extern void ix86_expand_lfloorceil (rtx, rtx, bool);
extern void ix86_expand_rint (rtx, rtx);
extern void ix86_expand_floorceil (rtx, rtx, bool);
extern void ix86_expand_floorceildf_32 (rtx, rtx, bool);
extern void ix86_expand_trunc (rtx, rtx);
extern void ix86_expand_truncdf_32 (rtx, rtx);
extern void ix86_expand_round (rtx, rtx);
extern void ix86_expand_rounddf_32 (rtx, rtx);
extern void ix86_expand_round_sse4 (rtx, rtx);

extern void ix86_expand_vecop_qihi (enum rtx_code, rtx, rtx, rtx);
extern rtx ix86_split_stack_guard (void);

extern void ix86_move_vector_high_sse_to_mmx (rtx);
extern void ix86_split_mmx_pack (rtx[], enum rtx_code);
extern void ix86_split_mmx_punpck (rtx[], bool);
extern void ix86_expand_avx_vzeroupper (void);
extern void ix86_expand_atomic_fetch_op_loop (rtx, rtx, rtx, enum rtx_code,
					      bool, bool);
extern void ix86_expand_cmpxchg_loop (rtx *, rtx, rtx, rtx, rtx, rtx,
				      bool, rtx_code_label *);

#ifdef TREE_CODE
extern void init_cumulative_args (CUMULATIVE_ARGS *, tree, rtx, tree, int);
#endif	/* TREE_CODE  */

#endif	/* RTX_CODE  */

#ifdef TREE_CODE
extern int ix86_data_alignment (tree, unsigned int, bool);
extern unsigned int ix86_local_alignment (tree, machine_mode,
					  unsigned int, bool = false);
extern unsigned int ix86_minimum_alignment (tree, machine_mode,
					    unsigned int);
extern tree ix86_handle_shared_attribute (tree *, tree, tree, int, bool *);
extern tree ix86_handle_selectany_attribute (tree *, tree, tree, int, bool *);
extern int x86_field_alignment (tree, int);
extern tree ix86_valid_target_attribute_tree (tree, tree,
					      struct gcc_options *,
					      struct gcc_options *, bool);
extern unsigned int ix86_get_callcvt (const_tree);

#endif

extern rtx ix86_tls_module_base (void);
extern bool ix86_gpr_tls_address_pattern_p (rtx);
extern bool ix86_tls_address_pattern_p (rtx);
extern rtx ix86_rewrite_tls_address (rtx);

extern void ix86_expand_vector_init (bool, rtx, rtx);
extern void ix86_expand_vector_set (bool, rtx, rtx, int);
extern void ix86_expand_vector_set_var (rtx, rtx, rtx);
extern void ix86_expand_vector_extract (bool, rtx, rtx, int);
extern void ix86_expand_reduc (rtx (*)(rtx, rtx, rtx), rtx, rtx);

extern void ix86_expand_vec_extract_even_odd (rtx, rtx, rtx, unsigned);
extern bool ix86_expand_pextr (rtx *);
extern bool ix86_expand_pinsr (rtx *);
extern void ix86_expand_mul_widen_evenodd (rtx, rtx, rtx, bool, bool);
extern void ix86_expand_mul_widen_hilo (rtx, rtx, rtx, bool, bool);
extern void ix86_expand_sse2_mulv4si3 (rtx, rtx, rtx);
extern void ix86_expand_sse2_mulvxdi3 (rtx, rtx, rtx);
extern void ix86_expand_sse2_abs (rtx, rtx);
extern bool ix86_expand_vector_init_duplicate (bool, machine_mode, rtx,
					       rtx);
extern bool ix86_extract_perm_from_pool_constant (int*, rtx);

/* In i386-c.cc  */
extern void ix86_target_macros (void);
extern void ix86_register_pragmas (void);

/* In winnt.cc  */
extern void i386_pe_unique_section (tree, int);
extern void i386_pe_declare_function_type (FILE *, const char *, int);
extern void i386_pe_record_external_function (tree, const char *);
extern void i386_pe_maybe_record_exported_symbol (tree, const char *, int);
extern void i386_pe_encode_section_info (tree, rtx, int);
extern bool i386_pe_binds_local_p (const_tree);
extern const char *i386_pe_strip_name_encoding_full (const char *);
extern bool i386_pe_valid_dllimport_attribute_p (const_tree);
extern unsigned int i386_pe_section_type_flags (tree, const char *, int);
extern void i386_pe_asm_named_section (const char *, unsigned int, tree);
extern void i386_pe_asm_output_aligned_decl_common (FILE *, tree,
						    const char *,
						    HOST_WIDE_INT,
						    HOST_WIDE_INT);
extern void i386_pe_file_end (void);
extern void i386_pe_asm_lto_start (void);
extern void i386_pe_asm_lto_end (void);
extern void i386_pe_start_function (FILE *, const char *, tree);
extern void i386_pe_end_function (FILE *, const char *, tree);
extern void i386_pe_end_cold_function (FILE *, const char *, tree);
extern void i386_pe_assemble_visibility (tree, int);
extern tree i386_pe_mangle_decl_assembler_name (tree, tree);
extern tree i386_pe_mangle_assembler_name (const char *);
extern void i386_pe_record_stub (const char *);

extern void i386_pe_seh_init (FILE *);
extern void i386_pe_seh_end_prologue (FILE *);
extern void i386_pe_seh_cold_init (FILE *, const char *);
extern void i386_pe_seh_unwind_emit (FILE *, rtx_insn *);
extern void i386_pe_seh_emit_except_personality (rtx);
extern void i386_pe_seh_init_sections (void);

/* In winnt-cxx.cc and winnt-stubs.cc  */
extern void i386_pe_adjust_class_at_definition (tree);
extern bool i386_pe_type_dllimport_p (tree);
extern bool i386_pe_type_dllexport_p (tree);

extern int i386_pe_reloc_rw_mask (void);

extern char internal_label_prefix[16];
extern int internal_label_prefix_len;

extern bool ix86_epilogue_uses (int);

struct ix86_address
{
  rtx base, index, disp;
  HOST_WIDE_INT scale;
  addr_space_t seg;
};

extern bool ix86_decompose_address (rtx, struct ix86_address *);
extern int memory_address_length (rtx, bool);
extern void x86_output_aligned_bss (FILE *, tree, const char *,
				    unsigned HOST_WIDE_INT, unsigned);
extern void x86_elf_aligned_decl_common (FILE *, tree, const char *,
					 unsigned HOST_WIDE_INT, unsigned);

#ifdef RTX_CODE
extern void ix86_fp_comparison_codes (enum rtx_code code, enum rtx_code *,
				      enum rtx_code *, enum rtx_code *);
extern enum rtx_code ix86_fp_compare_code_to_integer (enum rtx_code);
#endif
extern int asm_preferred_eh_data_format (int, int);

#ifdef HAVE_ATTR_cpu
extern enum attr_cpu ix86_schedule;
#endif

extern const char * ix86_output_call_insn (rtx_insn *insn, rtx call_op);
extern const char * ix86_output_indirect_jmp (rtx call_op);
extern const char * ix86_output_function_return (bool long_p);
extern const char * ix86_output_indirect_function_return (rtx ret_op);
extern void ix86_split_simple_return_pop_internal (rtx);
extern bool ix86_operands_ok_for_move_multiple (rtx *operands, bool load,
						machine_mode mode);
extern int ix86_min_insn_size (rtx_insn *);

extern int ix86_issue_rate (void);
extern int ix86_adjust_cost (rtx_insn *insn, int dep_type, rtx_insn *dep_insn,
			     int cost, unsigned int);
extern int ia32_multipass_dfa_lookahead (void);
extern bool ix86_macro_fusion_p (void);
extern bool ix86_macro_fusion_pair_p (rtx_insn *condgen, rtx_insn *condjmp);

extern bool ix86_bd_has_dispatch (rtx_insn *insn, int action);
extern void ix86_bd_do_dispatch (rtx_insn *insn, int mode);

extern void ix86_core2i7_init_hooks (void);

extern int ix86_atom_sched_reorder (FILE *, int, rtx_insn **, int *, int);

extern poly_int64 ix86_push_rounding (poly_int64);

#ifdef RTX_CODE
/* Target data for multipass lookahead scheduling.
   Currently used for Core 2/i7 tuning.  */
struct ix86_first_cycle_multipass_data_
{
  /* The length (in bytes) of ifetch block in this solution.  */
  int ifetch_block_len;
  /* Number of instructions in ifetch block in this solution.  */
  int ifetch_block_n_insns;
  /* Bitmap to remember changes to ready_try for backtracking.  */
  sbitmap ready_try_change;
  /* Size of the bitmap.  */
  int ready_try_change_size;
};
# define TARGET_SCHED_FIRST_CYCLE_MULTIPASS_DATA_T	\
  struct ix86_first_cycle_multipass_data_
#endif /* RTX_CODE */

const addr_space_t ADDR_SPACE_SEG_FS = 1;
const addr_space_t ADDR_SPACE_SEG_GS = 2;

namespace gcc { class context; }
class rtl_opt_pass;

extern rtl_opt_pass *make_pass_insert_vzeroupper (gcc::context *);
extern rtl_opt_pass *make_pass_stv (gcc::context *);
extern rtl_opt_pass *make_pass_insert_endbr_and_patchable_area
  (gcc::context *);
extern rtl_opt_pass *make_pass_remove_partial_avx_dependency
  (gcc::context *);

extern bool ix86_has_no_direct_extern_access;

/* In i386-expand.cc.  */
bool ix86_check_builtin_isa_match (unsigned int, HOST_WIDE_INT*,
				   HOST_WIDE_INT*);
