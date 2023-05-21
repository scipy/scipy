import numpy as np


def deprecate_ufunc(ufunc, old_name, new_name, message):
    """
    Deprecate a ufunc in scipy.special.

    This function is based on `np.deprecate`, but it changes how the
    beginning of the docstring is modified.  The ufunc signature
    remains the first line of the deprecated docstring.
    """
    new_func = np.deprecate(ufunc, old_name=old_name, new_name=new_name,
                            message=message)
    if new_func.__doc__ is not None:
        # Rearrange the first few lines of the docstring so the ufunc-generated
        # signature remains the first line.  This avoids a Sphinx warning when
        # the docs are build.
        lines = new_func.__doc__.split('\n')
        num_msg_lines = message.count('\n') + 1
        ufunc_sig_loc = 2 + num_msg_lines
        new_func.__doc__ = '\n'.join(lines[ufunc_sig_loc:ufunc_sig_loc + 4] +
                                     lines[:ufunc_sig_loc] +
                                     lines[ufunc_sig_loc + 4:])
    return new_func
