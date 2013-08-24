from __future__ import division, print_function, absolute_import

from numpy import asarray
import numpy as np


def abut(source, *args):
    # comment: except for the repetition, this is equivalent to hstack.
    """\nLike the |Stat abut command.  It concatenates two arrays column-wise
    and returns the result.  CAUTION:  If one array is shorter, it will be
    repeated until it is as long as the other.

    Format:  abut (source, args)    where args=any # of arrays
    Returns: an array as long as the LONGEST array past, source appearing on the
    'left', arrays in <args> attached on the 'right'.\n"""

    source = asarray(source)
    if len(source.shape) == 1:
        width = 1
        source = np.resize(source,[source.shape[0],width])
    else:
        width = source.shape[1]
    for addon in args:
        if len(addon.shape) == 1:
            width = 1
            addon = np.resize(addon,[source.shape[0],width])
        else:
            width = source.shape[1]
        if len(addon) < len(source):
            addon = np.resize(addon,[source.shape[0],addon.shape[1]])
        elif len(source) < len(addon):
            source = np.resize(source,[addon.shape[0],source.shape[1]])
        source = np.concatenate((source,addon),1)
    return source


def unique(inarray):
    """Returns unique items in the FIRST dimension of the passed array. Only
    works on arrays NOT including string items (e.g., type 'O' or 'c').
    """
    inarray = asarray(inarray)
    uniques = np.array([inarray[0]])
    if len(uniques.shape) == 1:            # IF IT'S A 1D ARRAY
        for item in inarray[1:]:
            if np.add.reduce(np.equal(uniques,item).flat) == 0:
                try:
                    uniques = np.concatenate([uniques,np.array[np.newaxis,:]])
                except TypeError:
                    uniques = np.concatenate([uniques,np.array([item])])
    else:                                  # IT MUST BE A 2+D ARRAY
        if inarray.dtype.char != 'O':  # not an Object array
            for item in inarray[1:]:
                if not np.sum(np.alltrue(np.equal(uniques,item),1),axis=0):
                    try:
                        uniques = np.concatenate([uniques,item[np.newaxis,:]])
                    except TypeError:    # the item to add isn't a list
                        uniques = np.concatenate([uniques,np.array([item])])
                else:
                    pass  # this item is already in the uniques array
        else:   # must be an Object array, alltrue/equal functions don't work
            for item in inarray[1:]:
                newflag = 1
                for unq in uniques:  # NOTE: cmp --> 0=same, -1=<, 1=>
                    test = np.sum(abs(np.array(list(map(cmp,item,unq)))),axis=0)
                    if test == 0:   # if item identical to any 1 row in uniques
                        newflag = 0  # then not a novel item to add
                        break
                if newflag == 1:
                    try:
                        uniques = np.concatenate([uniques,item[np.newaxis,:]])
                    except TypeError:    # the item to add isn't a list
                        uniques = np.concatenate([uniques,np.array([item])])
    return uniques


def _chk_asarray(a, axis):
    if axis is None:
        a = np.ravel(a)
        outaxis = 0
    else:
        a = np.asarray(a)
        outaxis = axis
    return a, outaxis


def _chk2_asarray(a, b, axis):
    if axis is None:
        a = np.ravel(a)
        b = np.ravel(b)
        outaxis = 0
    else:
        a = np.asarray(a)
        b = np.asarray(b)
        outaxis = axis
    return a, b, outaxis
