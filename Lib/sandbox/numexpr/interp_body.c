{
#define VEC_LOOP(expr) for(j = 0; j < VECTOR_SIZE; j++) {       \
        expr;                                       \
    }
#define VEC_ARG1(expr)                          \
    BOUNDS_CHECK(store_in);                     \
    BOUNDS_CHECK(arg1);                         \
    {                                           \
        char *dest = params.mem[store_in];      \
        char *x1 = params.mem[arg1];            \
        VEC_LOOP(expr);                         \
    } break
#define VEC_ARG2(expr)                          \
    BOUNDS_CHECK(store_in);                     \
    BOUNDS_CHECK(arg1);                         \
    BOUNDS_CHECK(arg2);                         \
    {                                           \
        char *dest = params.mem[store_in];      \
        char *x1 = params.mem[arg1];            \
        char *x2 = params.mem[arg2];            \
        VEC_LOOP(expr);                         \
    } break

#define VEC_ARG3(expr)                          \
    BOUNDS_CHECK(store_in);                     \
    BOUNDS_CHECK(arg1);                         \
    BOUNDS_CHECK(arg2);                         \
    BOUNDS_CHECK(arg3);                         \
    {                                           \
        char *dest = params.mem[store_in];      \
        char *x1 = params.mem[arg1];            \
        char *x2 = params.mem[arg2];            \
        char *x3 = params.mem[arg3];            \
        VEC_LOOP(expr);                         \
    } break


    unsigned int pc, j, r;
    /* set up pointers to next block of inputs and outputs */
    params.mem[0] = params.output + index * params.memsteps[0];
    for (r = 0; r < params.n_inputs; r++) {
        params.mem[1+r] = params.inputs[r] + index * params.memsteps[1+r];
    }
    for (pc = 0; pc < params.prog_len; pc += 4) {
        unsigned char op = params.program[pc];
        unsigned int store_in = params.program[pc+1];
        unsigned int arg1 = params.program[pc+2];
        unsigned int arg2 = params.program[pc+3];
        #define      arg3   params.program[pc+5]
        #define b_dest ((char *)dest)[j]
        #define i_dest ((long *)dest)[j]
        #define f_dest ((double *)dest)[j]
        #define cr_dest ((double *)dest)[2*j]
        #define ci_dest ((double *)dest)[2*j+1]
        #define b1    ((char   *)x1)[j]
        #define i1    ((long   *)x1)[j]
        #define f1    ((double *)x1)[j]
        #define c1r   ((double *)x1)[2*j]
        #define c1i   ((double *)x1)[2*j+1]
        #define b2    ((char   *)x2)[j]
        #define i2    ((long   *)x2)[j]
        #define f2    ((double *)x2)[j]
        #define c2r   ((double *)x2)[2*j]
        #define c2i   ((double *)x2)[2*j+1]
        #define b3    ((char   *)x3)[j]
        #define i3    ((long   *)x3)[j]
        #define f3    ((double *)x3)[j]
        #define c3r   ((double *)x3)[2*j]
        #define c3i   ((double *)x3)[2*j+1]
        double fa, fb;
        cdouble ca, cb;

        switch (op) {

        case OP_NOOP: break;

        /* The COPY are the only ops that depend on stride */
        case OP_COPY_BB: {
            char *dst = params.mem[store_in];
            char *src = params.mem[arg1];
            intp str1 = ((arg1 <= params.n_inputs) ? params.memsteps[arg1]
                                            : params.memsizes[arg1]);
            VEC_ARG1(memcpy(dst, src, sizeof(char));
                     dst += sizeof(char); src += str1);
            }
        case OP_COPY_II: {
            char *dst = params.mem[store_in];
            char *src = params.mem[arg1];
            intp str1 = ((arg1 <= params.n_inputs) ? params.memsteps[arg1]
                                            : params.memsizes[arg1]);
            VEC_ARG1(memcpy(dst, src, sizeof(long));
                     dst += sizeof(long); src += str1);
            }
        case OP_COPY_FF: {
            char *dst = params.mem[store_in];
            char *src = params.mem[arg1];
            intp str1 = ((arg1 <= params.n_inputs) ? params.memsteps[arg1]
                                            : params.memsizes[arg1]);
            VEC_ARG1(memcpy(dst, src, sizeof(double));
                     dst += sizeof(double); src += str1);
            }
        case OP_COPY_CC: {
            char *dst = params.mem[store_in];
            char *src = params.mem[arg1];
            intp str1 = ((arg1 <= params.n_inputs) ? params.memsteps[arg1]
                                            : params.memsizes[arg1]);
            VEC_ARG1(memcpy(dst, src, sizeof(double)*2);
                     dst += sizeof(double)*2; src += str1);
            }

        case OP_INVERT_BB: VEC_ARG1(b_dest = !b1);

        case OP_AND_BBB: VEC_ARG2(b_dest = b1 && b2);
        case OP_OR_BBB: VEC_ARG2(b_dest = b1 || b2);

        case OP_GT_BII: VEC_ARG2(b_dest = (i1 > i2) ? 1 : 0);
        case OP_GE_BII: VEC_ARG2(b_dest = (i1 >= i2) ? 1 : 0);
        case OP_EQ_BII: VEC_ARG2(b_dest = (i1 == i2) ? 1 : 0);
        case OP_NE_BII: VEC_ARG2(b_dest = (i1 != i2) ? 1 : 0);

        case OP_GT_BFF: VEC_ARG2(b_dest = (f1 > f2) ? 1 : 0);
        case OP_GE_BFF: VEC_ARG2(b_dest = (f1 >= f2) ? 1 : 0);
        case OP_EQ_BFF: VEC_ARG2(b_dest = (f1 == f2) ? 1 : 0);
        case OP_NE_BFF: VEC_ARG2(b_dest = (f1 != f2) ? 1 : 0);

        case OP_CAST_IB: VEC_ARG1(i_dest = (long)b1);
        case OP_ONES_LIKE_II: VEC_ARG1(i_dest = 1);
        case OP_NEG_II: VEC_ARG1(i_dest = -i1);

        case OP_ADD_III: VEC_ARG2(i_dest = i1 + i2);
        case OP_SUB_III: VEC_ARG2(i_dest = i1 - i2);
        case OP_MUL_III: VEC_ARG2(i_dest = i1 * i2);
        case OP_DIV_III: VEC_ARG2(i_dest = i1 / i2);
        case OP_POW_III: VEC_ARG2(i_dest = (i2 < 0) ? (1 / i1) : (long)pow(i1, i2));
        case OP_MOD_III: VEC_ARG2(i_dest = i1 % i2);

        case OP_WHERE_IFII: VEC_ARG3(i_dest = f1 ? i2 : i3);

        case OP_CAST_FB: VEC_ARG1(f_dest = (long)b1);
        case OP_CAST_FI: VEC_ARG1(f_dest = (double)(i1));
        case OP_ONES_LIKE_FF: VEC_ARG1(f_dest = 1.0);
        case OP_NEG_FF: VEC_ARG1(f_dest = -f1);

        case OP_ADD_FFF: VEC_ARG2(f_dest = f1 + f2);
        case OP_SUB_FFF: VEC_ARG2(f_dest = f1 - f2);
        case OP_MUL_FFF: VEC_ARG2(f_dest = f1 * f2);
        case OP_DIV_FFF: VEC_ARG2(f_dest = f1 / f2);
        case OP_POW_FFF: VEC_ARG2(f_dest = pow(f1, f2));
        case OP_MOD_FFF: VEC_ARG2(f_dest = f1 - floor(f1/f2) * f2);

        case OP_SIN_FF: VEC_ARG1(f_dest = sin(f1));
        case OP_COS_FF: VEC_ARG1(f_dest = cos(f1));
        case OP_TAN_FF: VEC_ARG1(f_dest = tan(f1));
        case OP_SQRT_FF: VEC_ARG1(f_dest = sqrt(f1));
        case OP_ARCTAN2_FFF: VEC_ARG2(f_dest = atan2(f1, f2));

        case OP_WHERE_FFFF: VEC_ARG3(f_dest = f1 ? f2 : f3);

        case OP_FUNC_FF: VEC_ARG1(f_dest = functions_f[arg2](f1));
        case OP_FUNC_FFF: VEC_ARG2(f_dest = functions_ff[arg3](f1, f2));

        case OP_CAST_CB: VEC_ARG1(cr_dest = (double)b1;
                                  ci_dest = 0);
        case OP_CAST_CI: VEC_ARG1(cr_dest = (double)(i1);
                                  ci_dest = 0);
        case OP_CAST_CF: VEC_ARG1(cr_dest = f1;
                                  ci_dest = 0);
        case OP_ONES_LIKE_CC: VEC_ARG1(cr_dest = 1;
                                  ci_dest = 0);
        case OP_NEG_CC: VEC_ARG1(cr_dest = -c1r;
                                 ci_dest = -c1i);

        case OP_ADD_CCC: VEC_ARG2(cr_dest = c1r + c2r;
                                  ci_dest = c1i + c2i);
        case OP_SUB_CCC: VEC_ARG2(cr_dest = c1r - c2r;
                                  ci_dest = c1i - c2i);
        case OP_MUL_CCC: VEC_ARG2(fa = c1r*c2r - c1i*c2i;
                                  ci_dest = c1r*c2i + c1i*c2r;
                                  cr_dest = fa);
        case OP_DIV_CCC: VEC_ARG2(fa = c2r*c2r + c2i*c2i;
                                  fb = (c1r*c2r + c1i*c2i) / fa;
                                  ci_dest = (c1i*c2r - c1r*c2i) / fa;
                                  cr_dest = fb);

        case OP_EQ_BCC: VEC_ARG2(b_dest = (c1r == c2r && c1i == c2i) ? 1 : 0);
        case OP_NE_BCC: VEC_ARG2(b_dest = (c1r != c2r || c1i != c2i) ? 1 : 0);

        case OP_WHERE_CFCC: VEC_ARG3(cr_dest = f1 ? c2r : c3r;
                                     ci_dest = f1 ? c2i : c3i);
        case OP_FUNC_CC: VEC_ARG1(ca.real = c1r;
                                  ca.imag = c1i;
                                  functions_cc[arg2](&ca, &ca);
                                  cr_dest = ca.real;
                                  ci_dest = ca.imag);
        case OP_FUNC_CCC: VEC_ARG2(ca.real = c1r;
                                   ca.imag = c1i;
                                   cb.real = c2r;
                                   cb.imag = c2i;
                                   functions_ccc[arg3](&ca, &cb, &ca);
                                   cr_dest = ca.real;
                                   ci_dest = ca.imag);

        case OP_REAL_FC: VEC_ARG1(f_dest = c1r);
        case OP_IMAG_FC: VEC_ARG1(f_dest = c1i);
        case OP_COMPLEX_CFF: VEC_ARG2(cr_dest = f1;
                                      ci_dest = f2);

        default:
            *pc_error = pc;
            return -3;
            break;
        }
    }

#undef VEC_LOOP
#undef VEC_ARG1
#undef VEC_ARG2
#undef VEC_ARG3

#undef b_dest
#undef i_dest
#undef f_dest
#undef cr_dest
#undef ci_dest
#undef b1
#undef i1
#undef f1
#undef c1r
#undef c1i
#undef b2
#undef i2
#undef f2
#undef c2r
#undef c2i
#undef b3
#undef i3
#undef f3
#undef c3r
#undef c3i
}
