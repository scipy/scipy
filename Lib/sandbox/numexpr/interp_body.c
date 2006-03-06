{
#define VEC_LOOP(expr) for(j = 0; j < VECTOR_SIZE; j++) {       \
        p_dest[j] = expr;                                       \
    }
#define VEC_ARG1(expr) do {                     \
        VEC_LOOP(expr);                         \
        break;                                  \
    } while (0)
#define VEC_ARG2_C(expr)                           \
    do {                                           \
        double c = constants[arg2];                \
        VEC_LOOP(expr);                            \
        break;                                     \
    } while (0)
#define VEC_ARG2(expr) do {                    \
        double *p2 = mem[arg2];                \
        VEC_LOOP(expr);                        \
        break;                                 \
    } while (0)

    unsigned int p, j, r;
    /* set up pointers to next block of inputs and outputs */
    mem[0] = output + index;
    for (r = 0; r < n_inputs; r++) {
        mem[1+r] = inputs[r] + index;;
    }
    for (p = 0; p < prog_len; p += 4) {
        char op = program[p];
        int store_in = program[p+1];
        int arg1 = program[p+2];
        double *p_dest = mem[store_in];
        double *p1 = mem[arg1];
        int arg2 = program[p+3];
        switch (op) {
        case OP_NOOP:
            break;
        case OP_COPY: VEC_ARG1(p1[j]);
        case OP_COPY_C: VEC_ARG2_C(c);
        case OP_NEG: VEC_ARG1(-p1[j]);
        case OP_ADD: VEC_ARG2(p1[j] + p2[j]);
        case OP_SUB: VEC_ARG2(p1[j] - p2[j]);
        case OP_MUL: VEC_ARG2(p1[j] * p2[j]);
        case OP_DIV: VEC_ARG2(p1[j] / p2[j]);
        case OP_POW: VEC_ARG2(pow(p1[j], p2[j]));
        case OP_MOD: VEC_ARG2(fmod(p1[j], p2[j]));
        case OP_ADD_C: VEC_ARG2_C(c + p1[j]);
        case OP_SUB_C: VEC_ARG2_C(c - p1[j]);
        case OP_MUL_C: VEC_ARG2_C(c * p1[j]);
        case OP_DIV_C: VEC_ARG2_C(c / p1[j]);
        case OP_POW_C: VEC_ARG2_C(pow(p1[j], c));
        case OP_MOD_C: VEC_ARG2_C(fmod(p1[j], c));
        case OP_GT: VEC_ARG2((p1[j] > p2[j]) ? 1 : 0);
        case OP_GE: VEC_ARG2((p1[j] >= p2[j]) ? 1 : 0);
        case OP_EQ: VEC_ARG2((p1[j] == p2[j]) ? 1 : 0);
        case OP_NE: VEC_ARG2((p1[j] != p2[j]) ? 1 : 0);
        case OP_GT_C: VEC_ARG2_C((p1[j] > c) ? 1 : 0);
        case OP_GE_C: VEC_ARG2_C((p1[j] >= c) ? 1 : 0);
        case OP_EQ_C: VEC_ARG2_C((p1[j] == c) ? 1 : 0);
        case OP_NE_C: VEC_ARG2_C((p1[j] != c) ? 1 : 0);
        case OP_LT_C: VEC_ARG2_C((p1[j] < c) ? 1 : 0);
        case OP_LE_C: VEC_ARG2_C((p1[j] <= c) ? 1 : 0);
        case OP_SIN: VEC_ARG1(sin(p1[j]));
        case OP_COS: VEC_ARG1(cos(p1[j]));
        case OP_TAN: VEC_ARG1(tan(p1[j]));
        case OP_ARCTAN2: VEC_ARG2(atan2(p1[j], p2[j]));
        case OP_WHERE:
        {
            char next_op = program[p+4];
            int arg3 = program[p+5];
            double *p2 = mem[arg2];
            double *p3 = mem[arg3];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = p1[j] ? p2[j] : p3[j];
            }
            break;
        }
        case OP_WHERE_XXC:
        {
            char next_op = program[p+4];
            int arg3 = program[p+5];
            double *p2 = mem[arg2];
            double c = constants[arg3];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = p1[j] ? p2[j] : c;
            }
            break;
        }
        case OP_WHERE_XCX:
        {
            char next_op = program[p+4];
            int arg3 = program[p+5];
            double *p2 = mem[arg2];
            double c = constants[arg3];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = p1[j] ? c : p2[j];
            }
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
            char next_op = program[p+4];
            int arg3 = program[p+5];
            double *p2 = mem[arg2];
            Func2Ptr func = functions_2[arg3];
            VEC_LOOP(func(p1[j], p2[j]));
            break;
        }
        default:
            break;
        }
    }

#undef VEC_LOOP
#undef VEC_ARG1
#undef VEC_ARG2_C
#undef VEC_ARG2

}
