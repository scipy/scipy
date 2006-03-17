# Copyright (C) 2003-2005 Peter J. Verveer
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above
#    copyright notice, this list of conditions and the following
#    disclaimer in the documentation and/or other materials provided
#    with the distribution.
#
# 3. The name of the author may not be used to endorse or promote
#    products derived from this software without specific prior
#    written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
# OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import types
import numpy as numarray

def _extend_mode_to_code(mode):
    """Convert an extension mode to the corresponding integer code.
    """
    if mode == 'nearest':
        return 0
    elif mode == 'wrap':
        return 1
    elif mode == 'reflect':
        return 2
    elif mode == 'constant':
        return 4
    else:
        raise RuntimeError, 'boundary mode not supported'

def _normalize_sequence(input, rank, array_type = None):
    """If input is a scalar, create a sequence of length equal to the
    rank by duplicating the input. If input is a sequence,
    check if its length is equal to the length of array.
    """
    if (isinstance(input, (types.IntType, types.LongType,
                           types.FloatType))):
        normalized = [input] * rank
    else:
        normalized = list(input)
        if len(normalized) != rank:
            err = "sequence argument must have length equal to input rank"
            raise RuntimeError, err
    return normalized

import warnings
def _get_output(output, input, output_type = None, shape = None):
    if output_type != None:
        msg = "'output_type' argument is decrepated."
        msg += " Assign type to 'output' instead."
        raise RuntimeError, msg
        warnings.warn(msg, DeprecationWarning)
        if output == None:
            output = output_type
        elif ((type(output) is not type(types.TypeType)) or
              output.dtype != output_type):
            raise RuntimeError, "'output' type and 'output_type' not equal"
    if shape is None:
        shape = input.shape
    if output is None:
        output = numarray.zeros(shape, dtype = input.dtype)
        return_value = output
    elif type(output) in [type(types.TypeType), type(numarray.zeros((4,)).dtype)]:
        output = numarray.zeros(shape, dtype = output)
        return_value = output
    elif type(output) is types.StringType:
        output = numarray.typeDict[output]
        output = numarray.zeros(shape, dtype = output)
        return_value = output
    else:
        if output.shape != shape:
            raise RuntimeError, "output shape not correct"
        return_value = None
    return output, return_value

def _check_axis(axis, rank):
    if axis < 0:
        axis += rank
    if axis < 0 or axis >= rank:
        raise ValueError, 'invalid axis'
    return axis
