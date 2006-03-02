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
        case OP_COPY:
            for (j = 0; j < VECTOR_SIZE; j++) {
                p_dest[j] = p1[j];
            }
            break;
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
        default:
            break;
        }
    }
}
