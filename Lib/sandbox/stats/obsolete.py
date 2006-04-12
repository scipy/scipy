def shellsort(a):
    # fixme: remove; obsolete
    """
Shellsort algorithm.  Sorts a 1D-array.

Returns: sorted-a, sorting-index-vector (for original array)

Use fastsort for speed.
"""
    a = asarray(a)
    n = len(a)
    svec = a*1.0
    ivec = range(n)
    gap = n/2   # integer division needed
    while gap > 0:
        for i in range(gap,n):
            for j in range(i-gap,-1,-gap):
                while j>=0 and svec[j]>svec[j+gap]:
                    temp        = svec[j]
                    svec[j]     = svec[j+gap]
                    svec[j+gap] = temp
                    itemp       = ivec[j]
                    ivec[j]     = ivec[j+gap]
                    ivec[j+gap] = itemp
        gap = gap / 2  # integer division needed
#    svec is now sorted input vector, ivec has the order svec[i] = vec[ivec[i]]
    return array(svec), array(ivec)


def summult(array1, array2, axis=0):
    """
Multiplies elements in array1 and array2, element by element, and
returns the sum (along 'axis') of all resulting multiplications.
Axis can equal None (ravel array first), or an integer (the
axis over which to operate),
"""
    array1, array2 = map(asarray, (array1, array2))
    if axis is None:
        array1 = ravel(array1)
        array2 = ravel(array2)
        axis = 0
    return sum(array1*array2,axis)

def sumdiffsquared(a, b, axis=0):
    """
Takes pairwise differences of the values in arrays a and b, squares
these differences, and returns the sum of these squares.  Axis
can equal None (ravel array first), an integer (the axis over
which to operate).

Returns: sum[(a-b)**2]
"""

    a, b = _chk2_asarray(a, b, axis)
    return sum((a-b)**2,axis)


