{
#define VEC_LOOP(expr) for(j = 0; j < VECTOR_SIZE; j++) {       \
        expr;                                       \
    }
#define VEC_ARG1(expr)                          \
    {                                           \
        VEC_LOOP(expr);                         \
    } break
#define VEC_ARG2(expr)                          \
    {                                           \
        long   *i2 = params.mem[arg2];          \
        double *f2 = params.mem[arg2];          \
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
        long   *i_dest = params.mem[store_in];
        double *f_dest = params.mem[store_in]; 
        long   *i1 = params.mem[arg1];
        double *f1 = params.mem[arg1];
        unsigned int arg2 = params.program[pc+3];
        BOUNDS_CHECK(store_in);
        BOUNDS_CHECK(arg1);
        BOUNDS_CHECK(arg2);
        switch (op) {
            
        case OP_NOOP: break;
        
        case OP_COPY_II: VEC_ARG2(i_dest[j] = i1[j]);
        case OP_ONES_LIKE_II: VEC_ARG1(i_dest[j] = 1);
        case OP_NEG_II: VEC_ARG1(i_dest[j] = -i1[j]);
        case OP_ADD_III: VEC_ARG2(i_dest[j] = i1[j] + i2[j]);
        case OP_SUB_III: VEC_ARG2(i_dest[j] = i1[j] - i2[j]);
        case OP_MUL_III: VEC_ARG2(i_dest[j] = i1[j] * i2[j]);
        case OP_DIV_III: VEC_ARG2(i_dest[j] = i1[j] / i2[j]);
        case OP_POW_III: VEC_ARG2(i_dest[j] = (long)pow(f1[j], f2[j])); /* XXX just a placeholder for now */
        case OP_MOD_III: VEC_ARG2(i_dest[j] = i1[j] % i2[j]); 
        case OP_WHERE_IFII:
        {
            unsigned int arg3 = params.program[pc+5];
            long *i2 = params.mem[arg2];
            long *i3 = params.mem[arg3];
            BOUNDS_CHECK(arg2);
            BOUNDS_CHECK(arg3);
            VEC_LOOP(i_dest[j] = f1[j] ? i2[j] : i3[j]);
            break;
        }
        
        case OP_CAST_FI: VEC_ARG1(f_dest[j] = (double)(i1[j]));
        case OP_COPY_FF: VEC_ARG1(f_dest[j] = f1[j]);
        case OP_ONES_LIKE_FF: VEC_ARG1(f_dest[j] = 1.0);
        case OP_NEG_FF: VEC_ARG1(f_dest[j] = -f1[j]);
        case OP_ADD_FFF: VEC_ARG2(f_dest[j] = f1[j] + f2[j]);
        case OP_SUB_FFF: VEC_ARG2(f_dest[j] = f1[j] - f2[j]);
        case OP_MUL_FFF: VEC_ARG2(f_dest[j] = f1[j] * f2[j]);
        case OP_DIV_FFF: VEC_ARG2(f_dest[j] = f1[j] / f2[j]);
        case OP_POW_FFF: VEC_ARG2(f_dest[j] = pow(f1[j], f2[j]));
        case OP_MOD_FFF: VEC_ARG2(f_dest[j] = f1[j] - floor(f1[j]/f2[j]) * f2[j]);
        case OP_GT_IFF: VEC_ARG2(i_dest[j] = (f1[j] > f2[j]) ? 1 : 0);
        case OP_GE_IFF: VEC_ARG2(i_dest[j] = (f1[j] >= f2[j]) ? 1 : 0);
        case OP_EQ_IFF: VEC_ARG2(i_dest[j] = (f1[j] == f2[j]) ? 1 : 0);
        case OP_NE_IFF: VEC_ARG2(i_dest[j] = (f1[j] != f2[j]) ? 1 : 0);
        case OP_SIN_FF: VEC_ARG1(f_dest[j] = sin(f1[j]));
        case OP_COS_FF: VEC_ARG1(f_dest[j] = cos(f1[j]));
        case OP_TAN_FF: VEC_ARG1(f_dest[j] = tan(f1[j]));
        case OP_SQRT_FF: VEC_ARG1(f_dest[j] = sqrt(f1[j]));
        case OP_ARCTAN2_FFF: VEC_ARG2(f_dest[j] = atan2(f1[j], f2[j]));
        case OP_WHERE_FFFF:
        {
            unsigned int arg3 = params.program[pc+5];
            double *f2 = params.mem[arg2];
            double *f3 = params.mem[arg3];
            BOUNDS_CHECK(arg2);
            BOUNDS_CHECK(arg3);
            VEC_LOOP(f_dest[j] = f1[j] ? f2[j] : f3[j]);
            break;
        }
        case OP_FUNC_FF:
        {
            FuncFFPtr func = functions_f[arg2];
            VEC_LOOP(f_dest[j] = func(f1[j]));
            break;
        }
        case OP_FUNC_FFF:
        {
            unsigned int arg3 = params.program[pc+5];
            double *f2 = params.mem[arg2];
            FuncFFFPtr func = functions_ff[arg3];
            BOUNDS_CHECK(arg2);
            VEC_LOOP(f_dest[j] = func(f1[j], f2[j]));
            break;
        }
        case OP_ADD_CCC: VEC_ARG2(f_dest[2*j] = f1[2*j] + f2[2*j];
                                  f_dest[2*j+1] = f1[2*j+1] + f2[2*j+1]);
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
