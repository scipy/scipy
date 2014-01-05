"""
Ad-hoc code generation for firfilter.c 2d convolution functions.

Use python code to do loop unswitching instead of using the C compiler.

"""
from __future__ import division, print_function, absolute_import

from itertools import product

# outsize options
FULL = 2
SAME = 1
VALID = 0

# boundary options
CIRCULAR = 8
REFLECT = 4
PAD = 0 

# there are only two flip (convolve) options
FLIPPED = 16
NOT_FLIPPED = 0


def case_codegen(outsize, boundary, convolve):
    lines = []
    def addline(line):
        lines.append(line)

    # loop over image rows
    addline('for (m=0; m < Os0; m++) {')
    if outsize == FULL:
        rhs = 'm' if convolve else 'm - Nwin0 + 1'
    elif outsize == SAME:
        rhs = 'm + ((Nwin0-1) >> 1)' if convolve else 'm - ((Nwin0-1) >> 1)'
    elif outsize == VALID:
        rhs = 'm + Nwin0 - 1' if convolve else 'm'
    else:
        raise ValueError
    addline('new_m = %s;' % rhs)

    # loop over image columns
    addline('for (n=0; n < Os1; n++) {')
    addline('sum = out + m*outstr0 + n*outstr1;')
    addline('memset(sum, 0, type_size);')
    if outsize == FULL:
        rhs = 'n' if convolve else 'n - Nwin1 + 1'
    elif outsize == SAME:
        rhs = 'n + ((Nwin1-1) >> 1)' if convolve else 'n - ((Nwin1-1) >> 1)'
    elif outsize == VALID:
        rhs = 'n + Nwin1 - 1' if convolve else 'n'
    else:
        raise ValueError
    addline('new_n = %s;' % rhs)

    # sum over kernel
    addline('for (j=0; j < Nwin0; j++) {')
    rhs = 'new_m - j' if convolve else 'new_m + j'
    addline('ind0 = %s;' % rhs)
    if boundary == PAD:
        addline('bounds_pad_flag = 0;')
    addline('if (ind0 < 0) {')
    if boundary == REFLECT:
        addline('ind0 = -1 - ind0;')
    elif boundary == CIRCULAR:
        addline('ind0 = Ns0 + ind0;')
    elif boundary == PAD:
        addline('bounds_pad_flag = 1;')
    else:
        raise ValueError
    addline('} else if (ind0 >= Ns0) {')
    if boundary == REFLECT:
        addline('ind0 = Ns0 + Ns0 - 1 - ind0;')
    elif boundary == CIRCULAR:
        addline('ind0 = ind0 - Ns0;')
    elif boundary == PAD:
        addline('bounds_pad_flag = 1;')
    else:
        raise ValueError
    addline('}')
    addline('h0_memory = hvals + j * hstr0;')
    if boundary == PAD:
        addline('if (!bounds_pad_flag) {')
    addline('ind0_memory = in + ind0 * instr0;')
    if boundary == PAD:
        addline('}')

    # kernel columns
    addline('for (k=0; k < Nwin1; k++) {')
    if boundary == PAD:
        addline('if (bounds_pad_flag) {')
        addline('value = fillvalue;')
        addline('} else {')
    rhs = 'new_n - k' if convolve else 'new_n + k'
    addline('ind1 = %s;' % rhs)
    addline('if (ind1 < 0) {')
    if boundary == REFLECT:
        addline('ind1 = -1 - ind1;')
    elif boundary == CIRCULAR:
        addline('ind1 = Ns1 + ind1;')
    elif boundary == PAD:
        addline('bounds_pad_flag = 1;')
    else:
        raise ValueError
    addline('} else if (ind1 >= Ns1) {')
    if boundary == REFLECT:
        addline('ind1 = Ns1 + Ns1 - 1 - ind1;')
    elif boundary == CIRCULAR:
        addline('ind1 = ind1 - Ns1;')
    elif boundary == PAD:
        addline('bounds_pad_flag = 1;')
    else:
        raise ValueError
    addline('}')
    if boundary == PAD:
        addline('if (bounds_pad_flag) {')
        addline('value = fillvalue;')
        addline('} else {')
    addline('value = ind0_memory + ind1*instr1;')
    if boundary == PAD:
        addline('}')
        addline('bounds_pad_flag = 0;')
        addline('}')

    # multiply and accumulate unless we are looking at a padding of zero
    if boundary == PAD:
        addline('if (!(fillvalue_is_zero && value==fillvalue)) {')
    addline('mult_and_add(sum, h0_memory + k*hstr1, value);')
    if boundary == PAD:
        addline('}')

    # close the loops over the 2d kernel within the 2d image
    addline('}')
    addline('}')
    addline('}')
    addline('}')
    
    return lines


def ad_hoc_indent(lines, depth, spacing):
    indented_lines = []
    for line in lines:
        if line.startswith('}'):
            depth -= 1
        indented_line = (spacing * depth) + line
        indented_lines.append(indented_line)
        if line.endswith('{'):
            depth += 1
    return '\n'.join(indented_lines)


def main():

    # Define some description strings.
    d_outsize = {VALID : 'VALID', SAME : 'SAME', FULL : 'FULL'}
    d_boundary = {PAD : 'PAD', REFLECT : 'REFLECT', CIRCULAR : 'CIRCULAR'}
    d_convolve = {NOT_FLIPPED : 'NOT_FLIPPED', FLIPPED : 'FLIPPED'}

    outsizes = (FULL, SAME, VALID)
    boundaries = (CIRCULAR, REFLECT, PAD)
    convolves = (FLIPPED, NOT_FLIPPED)

    lines = []
    for triple in product(outsizes, boundaries, convolves):
        outsize, boundary, convolve = triple
        description = ' | '.join((
            'outsize=' + d_outsize[outsize],
            'boundary=' + d_boundary[boundary],
            'convolve=' + d_convolve[convolve]))
        lines.append('case %s: /* %s */' % (sum(triple), description))
        lines.append('{')
        lines.extend(case_codegen(outsize, boundary, convolve))
        lines.append('break;')
        lines.append('}')
    print(ad_hoc_indent(lines, 0, '  '))


if __name__ == '__main__':
    main()

