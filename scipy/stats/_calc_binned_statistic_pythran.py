import numpy as np
import builtins


# pythran export _create_binned_data(int64[:] or intc[:] or int[:],
#                                    int[:],
#                                    int[:,:] or float[:,:],
#                                    int)
def _create_binned_data(bin_numbers, unique_bin_numbers, values, vv):
    """ Create hashmap of bin ids to values in bins
    key: bin number
    value: list of binned data
    """
    bin_map = dict()
    for i in unique_bin_numbers:
        bin_map[i] = []
    for i in builtins.range(len(bin_numbers)):
        bin_map[bin_numbers[i]].append(values[vv, i])
    return bin_map


# pythran export _calc_binned_statistic(int,
#                                       int64[:] or intc[:] or int[:],
#                                       float[:,:],
#                                       int[:,:] or float[:,:],
#                                       str)
def _calc_binned_statistic(Vdim, bin_numbers, result, values, stat_func):
    unique_bin_numbers = np.unique(bin_numbers)
    if stat_func == 'std':
        func = np.std
    elif stat_func == 'median':
        func = np.median
    elif stat_func == 'min':
        func = np.min
    elif stat_func == 'max':
        func = np.max
    else:
        raise Exception('Exception: {stat_func} is not supported')
    
    for vv in builtins.range(Vdim):
        bin_map = _create_binned_data(bin_numbers, unique_bin_numbers,
                                      values, vv)
        for i in unique_bin_numbers:
            # if the stat_func is np.std, calc std only when binned data is 2
            # or more for speed up.
            if not (stat_func == 'std' and len(bin_map[i]) < 2):
                result[vv, i] = func(np.array(bin_map[i]))
