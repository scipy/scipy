{
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
        case OP_COPY:
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = p1[j];
            }
            break;
        case OP_COPY_C:
        {
            double c = constants[arg2];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = c;
            }
            break;
        }
        case OP_NEG:
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = -p1[j];
            }
            break;
        case OP_ADD:
        {
            double *p2 = mem[arg2];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = p1[j] + p2[j];
            }
            break;
        }
        case OP_SUB:
        {
            double *p2 = mem[arg2];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = p1[j] - p2[j];
            }
            break;
        }
        case OP_MUL:
        {
            double *p2 = mem[arg2];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = p1[j] * p2[j];
            }
            break;
        }
        case OP_DIV:
        {
            double *p2 = mem[arg2];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = p1[j] / p2[j];
            }
            break;
        }
        case OP_ADD_C:
        {
            double c = constants[arg2];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = c + p1[j];
            }
            break;
        }
        case OP_SUB_C:
        {
            double c = constants[arg2];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = c - p1[j];
            }
            break;
        }
        case OP_MUL_C:
        {
            double c = constants[arg2];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = c * p1[j];
            }
            break;
        }
        case OP_DIV_C:
        {
            double c = constants[arg2];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = c / p1[j];
            }
            break;
        }
        case OP_GT:
        {
            double *p2 = mem[arg2];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = (p1[j] > p2[j]) ? 1 : 0;
            }
            break;
        }
        case OP_GE:
        {
            double *p2 = mem[arg2];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = (p1[j] >= p2[j]) ? 1 : 0;
            }
            break;
        }
        case OP_EQ:
        {
            double *p2 = mem[arg2];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = (p1[j] == p2[j]) ? 1 : 0;
            }
            break;
        }
        case OP_NE:
        {
            double *p2 = mem[arg2];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = (p1[j] != p2[j]) ? 1 : 0;
            }
            break;
        }
        case OP_GT_C:
        {
            double c = constants[arg2];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = (p1[j] > c) ? 1 : 0;
            }
            break;
        }
        case OP_GE_C:
        {
            double c = constants[arg2];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = (p1[j] >= c) ? 1 : 0;
            }
            break;
        }
        case OP_EQ_C:
        {
            double c = constants[arg2];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = (p1[j] == c) ? 1 : 0;
            }
            break;
        }
        case OP_NE_C:
        {
            double c = constants[arg2];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = (p1[j] == c) ? 1 : 0;
            }
            break;
        }
        case OP_LT_C:
        {
            double c = constants[arg2];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = (p1[j] < c) ? 1 : 0;
            }
            break;
        }
        case OP_LE_C:
        {
            double c = constants[arg2];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = (p1[j] <= c) ? 1 : 0;
            }
            break;
        }
        case OP_SIN:
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = sin(p1[j]);
            }
            break;
        case OP_COS:
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = cos(p1[j]);
            }
            break;
        case OP_TAN:
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = tan(p1[j]);
            }
            break;
        case OP_ARCTAN2:
        {
            double *p2 = mem[arg2];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = atan2(p1[j], p2[j]);
            }
            break;
        }
        case OP_ARCTAN2_1C:
        {
            double c = constants[arg2];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = atan2(p1[j], c);
            }
            break;
        }
        case OP_ARCTAN2_C1:
        {
            double c = constants[arg2];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = atan2(c, 1);
            }
            break;
        }
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
        case OP_WHERE_11C:
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
        case OP_WHERE_1C1:
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
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = func(p1[j]);
            }
            break;
        }
        case OP_FUNC_2:
        {
            char next_op = program[p+4];
            int arg3 = program[p+5];
            double *p2 = mem[arg2];
            Func2Ptr func = functions_2[arg3];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = func(p1[j], p2[j]);
            }
            break;
        }
        case OP_FUNC_1C:
        {
            char next_op = program[p+4];
            int arg3 = program[p+5];
            double c = constants[arg2];
            Func2Ptr func = functions_2[arg3];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = func(p1[j], c);
            }
            break;
        }
        case OP_FUNC_C1:
        {
            char next_op = program[p+4];
            int arg3 = program[p+5];
            double c = constants[arg2];
            Func2Ptr func = functions_2[arg3];
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = func(c, p1[j]);
            }
            break;
        }
        default:
            break;
        }
    }
}
