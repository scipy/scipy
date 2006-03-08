{
#if UNROLL
    unsigned int jj;
#define VEC_LOOP(expr) for(jj = 0; jj < VECTOR_SIZE; jj += 2) { \
        double d1, d2;                                          \
        j = jj;                                                 \
        d1 = expr;                                              \
        j = jj+1;                                               \
        d2 = expr;                                              \
        p_dest[jj] = d1;                                        \
        p_dest[j] = d2;                                         \
    }
#else
#define VEC_LOOP(expr) for(j = 0; j < VECTOR_SIZE; j++) {       \
        p_dest[j] = expr;                                       \
    }
#endif
#define VEC_ARG1(expr)                          \
    {                                           \
        VEC_LOOP(expr);                         \
    } break
#define VEC_ARG2_C(expr)                        \
    {                                           \
        double c = params.mem[arg2][0];         \
        VEC_LOOP(expr);                         \
    } break
#define VEC_ARG2(expr)                          \
    {                                           \
        double *p2 = params.mem[arg2];          \
        VEC_LOOP(expr);                         \
    } break

    unsigned int pc, j, r;
    /* set up pointers to next block of inputs and outputs */
    params.mem[0] = params.output + index;
    for (r = 0; r < params.n_inputs; r++) {
        params.mem[1+r] = params.inputs[r] + index;
    }
    for (pc = 0; pc < params.prog_len; pc += 4) {
        unsigned char op = params.program[pc];
        unsigned int store_in = params.program[pc+1];
        BOUNDS_CHECK(store_in);
        unsigned int arg1 = params.program[pc+2];
        BOUNDS_CHECK(arg1);
        double *p_dest = params.mem[store_in];
        double *p1 = params.mem[arg1];
        unsigned int arg2 = params.program[pc+3];
        BOUNDS_CHECK(arg2);
        switch (op) {
        case OP_NOOP:
            break;
        case OP_COPY: VEC_ARG1(p1[j]);
        case OP_NEG: VEC_ARG1(-p1[j]);
        case OP_ADD: VEC_ARG2(p1[j] + p2[j]);
        case OP_SUB: VEC_ARG2(p1[j] - p2[j]);
        case OP_MUL: VEC_ARG2(p1[j] * p2[j]);
        case OP_DIV: VEC_ARG2(p1[j] / p2[j]);
        case OP_POW: VEC_ARG2(pow(p1[j], p2[j]));
        case OP_MOD: VEC_ARG2(fmod(p1[j], p2[j]));
        case OP_DIV_C: VEC_ARG2_C(c / p1[j]);
        case OP_GT: VEC_ARG2((p1[j] > p2[j]) ? 1 : 0);
        case OP_GE: VEC_ARG2((p1[j] >= p2[j]) ? 1 : 0);
        case OP_EQ: VEC_ARG2((p1[j] == p2[j]) ? 1 : 0);
        case OP_NE: VEC_ARG2((p1[j] != p2[j]) ? 1 : 0);
        case OP_SIN: VEC_ARG1(sin(p1[j]));
        case OP_COS: VEC_ARG1(cos(p1[j]));
        case OP_TAN: VEC_ARG1(tan(p1[j]));
        case OP_ARCTAN2: VEC_ARG2(atan2(p1[j], p2[j]));
        case OP_WHERE:
        {
            unsigned int arg3 = params.program[pc+5];
            double *p2 = params.mem[arg2];
            double *p3 = params.mem[arg3];
            VEC_LOOP(p1[j] ? p2[j] : p3[j]);
            break;
        }
        case OP_FUNC_1:
        {
            Func1Ptr func = functions_1[arg2];
            VEC_LOOP(func(p1[j]));
            break;
        }
        case OP_FUNC_2:
        {
            unsigned int arg3 = params.program[pc+5];
            double *p2 = params.mem[arg2];
            Func2Ptr func = functions_2[arg3];
            VEC_LOOP(func(p1[j], p2[j]));
            break;
        }
        default:
            *pc_error = pc;
            return -3;
            break;
        }
    }

#undef VEC_LOOP
#undef VEC_ARG1
#undef VEC_ARG2_C
#undef VEC_ARG2

}
