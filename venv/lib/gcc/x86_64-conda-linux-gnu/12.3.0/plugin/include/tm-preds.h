/* Generated automatically by the program 'build/genpreds'
   from the machine description file '../../gcc/config/i386/i386.md'.  */

#ifndef GCC_TM_PREDS_H
#define GCC_TM_PREDS_H

#ifdef HAVE_MACHINE_MODES
extern bool general_operand (rtx, machine_mode);
extern bool address_operand (rtx, machine_mode);
extern bool register_operand (rtx, machine_mode);
extern bool pmode_register_operand (rtx, machine_mode);
extern bool scratch_operand (rtx, machine_mode);
extern bool immediate_operand (rtx, machine_mode);
extern bool const_int_operand (rtx, machine_mode);
extern bool const_scalar_int_operand (rtx, machine_mode);
extern bool const_double_operand (rtx, machine_mode);
extern bool nonimmediate_operand (rtx, machine_mode);
extern bool nonmemory_operand (rtx, machine_mode);
extern bool push_operand (rtx, machine_mode);
extern bool pop_operand (rtx, machine_mode);
extern bool memory_operand (rtx, machine_mode);
extern bool indirect_operand (rtx, machine_mode);
extern bool ordered_comparison_operator (rtx, machine_mode);
extern bool comparison_operator (rtx, machine_mode);
extern bool any_fp_register_operand (rtx, machine_mode);
extern bool fp_register_operand (rtx, machine_mode);
extern bool general_reg_operand (rtx, machine_mode);
extern bool nonimmediate_gr_operand (rtx, machine_mode);
extern bool general_gr_operand (rtx, machine_mode);
extern bool mmx_reg_operand (rtx, machine_mode);
extern bool register_mmxmem_operand (rtx, machine_mode);
extern bool sse_reg_operand (rtx, machine_mode);
extern bool any_QIreg_operand (rtx, machine_mode);
extern bool QIreg_operand (rtx, machine_mode);
extern bool ext_QIreg_operand (rtx, machine_mode);
extern bool ax_reg_operand (rtx, machine_mode);
extern bool flags_reg_operand (rtx, machine_mode);
extern bool mask_reg_operand (rtx, machine_mode);
extern bool int_nonimmediate_operand (rtx, machine_mode);
extern bool register_ssemem_operand (rtx, machine_mode);
extern bool nonimm_ssenomem_operand (rtx, machine_mode);
extern bool x87nonimm_ssenomem_operand (rtx, machine_mode);
extern bool register_sse4nonimm_operand (rtx, machine_mode);
extern bool symbol_operand (rtx, machine_mode);
extern bool ix86_endbr_immediate_operand (rtx, machine_mode);
extern bool x86_64_immediate_operand (rtx, machine_mode);
extern bool x86_64_zext_immediate_operand (rtx, machine_mode);
extern bool x86_64_hilo_int_operand (rtx, machine_mode);
extern bool x86_64_dwzext_immediate_operand (rtx, machine_mode);
extern bool x86_64_immediate_size_operand (rtx, machine_mode);
extern bool x86_64_general_operand (rtx, machine_mode);
extern bool x86_64_hilo_general_operand (rtx, machine_mode);
extern bool x86_64_sext_operand (rtx, machine_mode);
extern bool sext_operand (rtx, machine_mode);
extern bool x86_64_zext_operand (rtx, machine_mode);
extern bool x86_64_szext_general_operand (rtx, machine_mode);
extern bool x86_64_nonmemory_operand (rtx, machine_mode);
extern bool x86_64_szext_nonmemory_operand (rtx, machine_mode);
extern bool pic_32bit_operand (rtx, machine_mode);
extern bool x86_64_movabs_operand (rtx, machine_mode);
extern bool symbolic_operand (rtx, machine_mode);
extern bool local_symbolic_operand (rtx, machine_mode);
extern bool gotoff_operand (rtx, machine_mode);
extern bool tls_symbolic_operand (rtx, machine_mode);
extern bool tls_modbase_operand (rtx, machine_mode);
extern bool tls_address_pattern (rtx, machine_mode);
extern bool constant_call_address_operand (rtx, machine_mode);
extern bool call_register_no_elim_operand (rtx, machine_mode);
extern bool register_no_elim_operand (rtx, machine_mode);
extern bool index_register_operand (rtx, machine_mode);
extern bool general_no_elim_operand (rtx, machine_mode);
extern bool nonmemory_no_elim_operand (rtx, machine_mode);
extern bool indirect_branch_operand (rtx, machine_mode);
extern bool sibcall_memory_operand (rtx, machine_mode);
extern bool GOT_memory_operand (rtx, machine_mode);
extern bool call_insn_operand (rtx, machine_mode);
extern bool sibcall_insn_operand (rtx, machine_mode);
extern bool GOT32_symbol_operand (rtx, machine_mode);
extern bool const0_operand (rtx, machine_mode);
extern bool const1_operand (rtx, machine_mode);
extern bool constm1_operand (rtx, machine_mode);
extern bool const8_operand (rtx, machine_mode);
extern bool const128_operand (rtx, machine_mode);
extern bool const_32bit_mask (rtx, machine_mode);
extern bool const248_operand (rtx, machine_mode);
extern bool const123_operand (rtx, machine_mode);
extern bool const2367_operand (rtx, machine_mode);
extern bool const1248_operand (rtx, machine_mode);
extern bool const359_operand (rtx, machine_mode);
extern bool const_4_or_8_to_11_operand (rtx, machine_mode);
extern bool const48_operand (rtx, machine_mode);
extern bool const_0_to_1_operand (rtx, machine_mode);
extern bool const_0_to_3_operand (rtx, machine_mode);
extern bool const_0_to_4_operand (rtx, machine_mode);
extern bool const_0_to_5_operand (rtx, machine_mode);
extern bool const_0_to_7_operand (rtx, machine_mode);
extern bool const_0_to_15_operand (rtx, machine_mode);
extern bool const_0_to_31_operand (rtx, machine_mode);
extern bool const_0_to_63_operand (rtx, machine_mode);
extern bool const_0_to_255_operand (rtx, machine_mode);
extern bool const_0_to_255_mul_8_operand (rtx, machine_mode);
extern bool const_1_to_31_operand (rtx, machine_mode);
extern bool const_1_to_63_operand (rtx, machine_mode);
extern bool const_2_to_3_operand (rtx, machine_mode);
extern bool const_4_to_5_operand (rtx, machine_mode);
extern bool const_4_to_7_operand (rtx, machine_mode);
extern bool const_6_to_7_operand (rtx, machine_mode);
extern bool const_8_to_9_operand (rtx, machine_mode);
extern bool const_8_to_11_operand (rtx, machine_mode);
extern bool const_8_to_15_operand (rtx, machine_mode);
extern bool const_10_to_11_operand (rtx, machine_mode);
extern bool const_12_to_13_operand (rtx, machine_mode);
extern bool const_12_to_15_operand (rtx, machine_mode);
extern bool const_14_to_15_operand (rtx, machine_mode);
extern bool const_16_to_19_operand (rtx, machine_mode);
extern bool const_16_to_31_operand (rtx, machine_mode);
extern bool const_20_to_23_operand (rtx, machine_mode);
extern bool const_24_to_27_operand (rtx, machine_mode);
extern bool const_28_to_31_operand (rtx, machine_mode);
extern bool incdec_operand (rtx, machine_mode);
extern bool vec_setm_sse41_operand (rtx, machine_mode);
extern bool vec_setm_avx2_operand (rtx, machine_mode);
extern bool vec_setm_mmx_operand (rtx, machine_mode);
extern bool reg_or_pm1_operand (rtx, machine_mode);
extern bool regmem_or_bitnot_regmem_operand (rtx, machine_mode);
extern bool shiftdi_operand (rtx, machine_mode);
extern bool ashldi_input_operand (rtx, machine_mode);
extern bool zero_extended_scalar_load_operand (rtx, machine_mode);
extern bool float_vector_all_ones_operand (rtx, machine_mode);
extern bool vector_all_ones_operand (rtx, machine_mode);
extern bool vector_memory_operand (rtx, machine_mode);
extern bool vector_operand (rtx, machine_mode);
extern bool bcst_mem_operand (rtx, machine_mode);
extern bool bcst_vector_operand (rtx, machine_mode);
extern bool nonimmediate_or_const_vector_operand (rtx, machine_mode);
extern bool nonimmediate_or_const_vec_dup_operand (rtx, machine_mode);
extern bool reg_or_const_vector_operand (rtx, machine_mode);
extern bool nonimmediate_or_sse_const_operand (rtx, machine_mode);
extern bool reg_or_0_operand (rtx, machine_mode);
extern bool nonimm_or_0_operand (rtx, machine_mode);
extern bool norex_memory_operand (rtx, machine_mode);
extern bool SImode_address_operand (rtx, machine_mode);
extern bool address_no_seg_operand (rtx, machine_mode);
extern bool vsib_address_operand (rtx, machine_mode);
extern bool vsib_mem_operator (rtx, machine_mode);
extern bool aligned_operand (rtx, machine_mode);
extern bool memory_displacement_operand (rtx, machine_mode);
extern bool memory_displacement_only_operand (rtx, machine_mode);
extern bool long_memory_operand (rtx, machine_mode);
extern bool fcmov_comparison_operator (rtx, machine_mode);
extern bool sse_comparison_operator (rtx, machine_mode);
extern bool ix86_comparison_int_operator (rtx, machine_mode);
extern bool ix86_comparison_uns_operator (rtx, machine_mode);
extern bool bt_comparison_operator (rtx, machine_mode);
extern bool shr_comparison_operator (rtx, machine_mode);
extern bool add_comparison_operator (rtx, machine_mode);
extern bool ix86_comparison_operator (rtx, machine_mode);
extern bool ix86_carry_flag_operator (rtx, machine_mode);
extern bool ix86_carry_flag_unset_operator (rtx, machine_mode);
extern bool ix86_trivial_fp_comparison_operator (rtx, machine_mode);
extern bool ix86_fp_comparison_operator (rtx, machine_mode);
extern bool cmp_fp_expander_operand (rtx, machine_mode);
extern bool binary_fp_operator (rtx, machine_mode);
extern bool mult_operator (rtx, machine_mode);
extern bool div_operator (rtx, machine_mode);
extern bool logic_operator (rtx, machine_mode);
extern bool plusminuslogic_operator (rtx, machine_mode);
extern bool arith_or_logical_operator (rtx, machine_mode);
extern bool commutative_operator (rtx, machine_mode);
extern bool promotable_binary_operator (rtx, machine_mode);
extern bool compare_operator (rtx, machine_mode);
extern bool misaligned_operand (rtx, machine_mode);
extern bool movq_parallel (rtx, machine_mode);
extern bool vzeroall_operation (rtx, machine_mode);
extern bool vzeroall_pattern (rtx, machine_mode);
extern bool vzeroupper_pattern (rtx, machine_mode);
extern bool addsub_vm_operator (rtx, machine_mode);
extern bool addsub_vs_operator (rtx, machine_mode);
extern bool addsub_vs_parallel (rtx, machine_mode);
extern bool permvar_truncate_operand (rtx, machine_mode);
extern bool pshufb_truncv4siv4hi_operand (rtx, machine_mode);
extern bool pshufb_truncv8hiv8qi_operand (rtx, machine_mode);
extern bool pmovzx_parallel (rtx, machine_mode);
extern bool const_vector_duplicate_operand (rtx, machine_mode);
extern bool avx_vbroadcast_operand (rtx, machine_mode);
extern bool palignr_operand (rtx, machine_mode);
extern bool avx2_pblendw_operand (rtx, machine_mode);
extern bool general_vector_operand (rtx, machine_mode);
extern bool register_or_constm1_operand (rtx, machine_mode);
extern bool save_multiple (rtx, machine_mode);
extern bool restore_multiple (rtx, machine_mode);
extern bool encodekey128_operation (rtx, machine_mode);
extern bool encodekey256_operation (rtx, machine_mode);
extern bool aeswidekl_operation (rtx, machine_mode);
#endif /* HAVE_MACHINE_MODES */

#define CONSTRAINT_NUM_DEFINED_P 1
enum constraint_num
{
  CONSTRAINT__UNKNOWN = 0,
  CONSTRAINT_r,
  CONSTRAINT_R,
  CONSTRAINT_q,
  CONSTRAINT_Q,
  CONSTRAINT_l,
  CONSTRAINT_a,
  CONSTRAINT_b,
  CONSTRAINT_c,
  CONSTRAINT_d,
  CONSTRAINT_S,
  CONSTRAINT_D,
  CONSTRAINT_A,
  CONSTRAINT_U,
  CONSTRAINT_f,
  CONSTRAINT_t,
  CONSTRAINT_u,
  CONSTRAINT_Yk,
  CONSTRAINT_k,
  CONSTRAINT_y,
  CONSTRAINT_x,
  CONSTRAINT_v,
  CONSTRAINT_Yz,
  CONSTRAINT_Yd,
  CONSTRAINT_Yp,
  CONSTRAINT_Ya,
  CONSTRAINT_Yb,
  CONSTRAINT_Yf,
  CONSTRAINT_Yr,
  CONSTRAINT_Yv,
  CONSTRAINT_Yw,
  CONSTRAINT_YW,
  CONSTRAINT_I,
  CONSTRAINT_J,
  CONSTRAINT_K,
  CONSTRAINT_L,
  CONSTRAINT_M,
  CONSTRAINT_N,
  CONSTRAINT_O,
  CONSTRAINT_m,
  CONSTRAINT_o,
  CONSTRAINT_Bk,
  CONSTRAINT_Bm,
  CONSTRAINT_Bc,
  CONSTRAINT_Bn,
  CONSTRAINT_Br,
  CONSTRAINT_p,
  CONSTRAINT_Tv,
  CONSTRAINT_Ts,
  CONSTRAINT_Bz,
  CONSTRAINT_Wb,
  CONSTRAINT_Ww,
  CONSTRAINT_G,
  CONSTRAINT_e,
  CONSTRAINT_We,
  CONSTRAINT_Wz,
  CONSTRAINT_Wd,
  CONSTRAINT_Wf,
  CONSTRAINT_Z,
  CONSTRAINT_Bf,
  CONSTRAINT_V,
  CONSTRAINT__l,
  CONSTRAINT__g,
  CONSTRAINT_BF,
  CONSTRAINT_BM,
  CONSTRAINT_i,
  CONSTRAINT_s,
  CONSTRAINT_n,
  CONSTRAINT_E,
  CONSTRAINT_F,
  CONSTRAINT_X,
  CONSTRAINT_Bg,
  CONSTRAINT_Bs,
  CONSTRAINT_Bw,
  CONSTRAINT_BC,
  CONSTRAINT_C,
  CONSTRAINT__LIMIT
};

extern enum constraint_num lookup_constraint_1 (const char *);
extern const unsigned char lookup_constraint_array[];

/* Return the constraint at the beginning of P, or CONSTRAINT__UNKNOWN if it
   isn't recognized.  */

static inline enum constraint_num
lookup_constraint (const char *p)
{
  unsigned int index = lookup_constraint_array[(unsigned char) *p];
  return (index == UCHAR_MAX
          ? lookup_constraint_1 (p)
          : (enum constraint_num) index);
}

extern bool (*constraint_satisfied_p_array[]) (rtx);

/* Return true if X satisfies constraint C.  */

static inline bool
constraint_satisfied_p (rtx x, enum constraint_num c)
{
  int i = (int) c - (int) CONSTRAINT_I;
  return i >= 0 && constraint_satisfied_p_array[i] (x);
}

static inline bool
insn_extra_register_constraint (enum constraint_num c)
{
  return c >= CONSTRAINT_r && c <= CONSTRAINT_YW;
}

static inline bool
insn_extra_memory_constraint (enum constraint_num c)
{
  return c >= CONSTRAINT_m && c <= CONSTRAINT_Bk;
}

static inline bool
insn_extra_special_memory_constraint (enum constraint_num c)
{
  return c >= CONSTRAINT_Bm && c <= CONSTRAINT_Br;
}

static inline bool
insn_extra_relaxed_memory_constraint (enum constraint_num)
{
  return false;
}

static inline bool
insn_extra_address_constraint (enum constraint_num c)
{
  return c >= CONSTRAINT_p && c <= CONSTRAINT_Ts;
}

static inline void
insn_extra_constraint_allows_reg_mem (enum constraint_num c,
				      bool *allows_reg, bool *allows_mem)
{
  if (c >= CONSTRAINT_Bz && c <= CONSTRAINT_Z)
    return;
  if (c >= CONSTRAINT_Bf && c <= CONSTRAINT_Bf)
    {
      *allows_reg = true;
      return;
    }
  if (c >= CONSTRAINT_V && c <= CONSTRAINT_BM)
    {
      *allows_mem = true;
      return;
    }
  (void) c;
  *allows_reg = true;
  *allows_mem = true;
}

static inline size_t
insn_constraint_len (char fc, const char *str ATTRIBUTE_UNUSED)
{
  switch (fc)
    {
    case 'B': return 2;
    case 'T': return 2;
    case 'W': return 2;
    case 'Y': return 2;
    default: break;
    }
  return 1;
}

#define CONSTRAINT_LEN(c_,s_) insn_constraint_len (c_,s_)

extern enum reg_class reg_class_for_constraint_1 (enum constraint_num);

static inline enum reg_class
reg_class_for_constraint (enum constraint_num c)
{
  if (insn_extra_register_constraint (c))
    return reg_class_for_constraint_1 (c);
  return NO_REGS;
}

extern bool insn_const_int_ok_for_constraint (HOST_WIDE_INT, enum constraint_num);
#define CONST_OK_FOR_CONSTRAINT_P(v_,c_,s_) \
    insn_const_int_ok_for_constraint (v_, lookup_constraint (s_))

enum constraint_type
{
  CT_REGISTER,
  CT_CONST_INT,
  CT_MEMORY,
  CT_SPECIAL_MEMORY,
  CT_RELAXED_MEMORY,
  CT_ADDRESS,
  CT_FIXED_FORM
};

static inline enum constraint_type
get_constraint_type (enum constraint_num c)
{
  if (c >= CONSTRAINT_Bm)
    {
      if (c >= CONSTRAINT_Bz)
        return CT_FIXED_FORM;
      if (c >= CONSTRAINT_p)
        return CT_ADDRESS;
      return CT_SPECIAL_MEMORY;
    }
  if (c >= CONSTRAINT_m)
    return CT_MEMORY;
  if (c >= CONSTRAINT_I)
    return CT_CONST_INT;
  return CT_REGISTER;
}
#endif /* tm-preds.h */
