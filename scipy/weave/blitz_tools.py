from __future__ import absolute_import, print_function

import parser
import sys
import warnings
import copy

import numpy

from . import ast_tools
from . import slice_handler
from . import size_check
from . import converters
from . import inline_tools
from .inline_tools import attempt_function_call
function_catalog = inline_tools.function_catalog
function_cache = inline_tools.function_cache


class BlitzWarning(UserWarning):
    """Warns about compilation failures etc."""
    pass


def blitz(expr,local_dict=None, global_dict=None,check_size=1,verbose=0,**kw):
    # this could call inline, but making a copy of the
    # code here is more efficient for several reasons.
    global function_catalog

    # this grabs the local variables from the *previous* call
    # frame -- that is the locals from the function that called
    # inline.
    call_frame = sys._getframe().f_back
    if local_dict is None:
        local_dict = call_frame.f_locals
    if global_dict is None:
        global_dict = call_frame.f_globals

    # 1. Check the sizes of the arrays and make sure they are compatible.
    #    This is expensive, so unsetting the check_size flag can save a lot
    #    of time.  It also can cause core-dumps if the sizes of the inputs
    #    aren't compatible.
    if check_size and not size_check.check_expr(expr,local_dict,global_dict):
        raise ValueError("inputs failed to pass size check.")

    # 2. try local cache
    try:
        results = apply(function_cache[expr],(local_dict,global_dict))
        return results
    except:
        pass
    try:
        results = attempt_function_call(expr,local_dict,global_dict)
    # 3. build the function
    except ValueError:
        # This section is pretty much the only difference
        # between blitz and inline
        ast = parser.suite(expr)
        ast_list = ast.tolist()
        expr_code = ast_to_blitz_expr(ast_list)
        arg_names = ast_tools.harvest_variables(ast_list)
        module_dir = global_dict.get('__file__',None)
        func = inline_tools.compile_function(expr_code,arg_names,local_dict,
                                             global_dict,module_dir,
                                             compiler='gcc',auto_downcast=1,
                                             verbose=verbose,
                                             type_converters=converters.blitz,
                                             **kw)
        function_catalog.add_function(expr,func,module_dir)
        try:
            results = attempt_function_call(expr,local_dict,global_dict)
        except ValueError:
            warnings.warn('compilation failed. Executing as python code',
                          BlitzWarning)
            exec(expr, global_dict, local_dict)


def ast_to_blitz_expr(ast_seq):
    """Convert an ast_sequence to a blitz expression."""
    # Don't overwrite orignal sequence in call to transform slices.
    ast_seq = copy.deepcopy(ast_seq)
    slice_handler.transform_slices(ast_seq)

    # Build the actual program statement from ast_seq
    expr = ast_tools.ast_to_string(ast_seq)

    # Now find and replace specific symbols to convert this to
    # a blitz++ compatible statement.
    # I'm doing this with string replacement here.  It could
    # also be done on the actual ast tree (and probably should from
    # a purest standpoint...).

    # this one isn't necessary but it helps code readability
    # and compactness. It requires that
    #   Range _all = blitz::Range::all();
    # be included in the generated code.
    # These could all alternatively be done to the ast in
    # build_slice_atom()
    expr = expr.replace('slice(_beg,_end)', '_all')
    expr = expr.replace('slice', 'blitz::Range')
    expr = expr.replace('[','(')
    expr = expr.replace(']', ')')
    expr = expr.replace('_stp', '1')

    # Instead of blitz::fromStart and blitz::toEnd.  This requires
    # the following in the generated code.
    #   Range _beg = blitz::fromStart;
    #   Range _end = blitz::toEnd;
    #expr = expr.replace('_beg', 'blitz::fromStart' )
    #expr = expr.replace('_end', 'blitz::toEnd' )

    return expr + ';\n'
