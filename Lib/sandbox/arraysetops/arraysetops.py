"""
Set operations for 1D numeric arrays based on sort() function.

Contains:
  ediff1d,
  unique1d,
  intersect1d,
  intersect1d_nu,
  setxor1d,
  setmember1d,
  union1d,
  setdiff1d

Concerning the speed, test_unique1d_speed() reveals that up to 10000000
elements unique1d() is about 10 times faster than the standard dictionary-based
numpy.unique().

Limitations: Except unique1d, union1d and intersect1d_nu, all functions expect
inputs with unique elements. Speed could be gained in some operations by an
implementaion of sort(), that can provide directly the permutation vectors,
avoiding thus calls to argsort().

To do: Optionally return indices analogously to unique1d for all functions.

Author: Robert Cimrman
"""
# 02.11.2005, c
import time
import scipy

##
# 03.11.2005, c
def ediff1d( ar1, toEnd = None, toBegin = None ):
    """Array difference with prefixed and/or appended value."""
    dar1 = ar1[1:] - ar1[:-1]
    if toEnd and toBegin:
        shape = (ar1.shape[0] + 1,) + ar1.shape[1:]
        ed = numpy.empty( shape, dtype = ar1.dtype )
        ed[0], ed[-1] = toBegin, toEnd
        ed[1:-1] = dar1
    elif toEnd:
        ed = numpy.empty( ar1.shape, dtype = ar1.dtype )
        ed[-1] = toEnd
        ed[:-1] = dar1
    elif toBegin:
        ed = numpy.empty( ar1.shape, dtype = ar1.dtype )
        ed[0] = toBegin
        ed[1:] = dar1
    else:
        ed = dar1

    return ed


##
# 01.11.2005, c
# 02.11.2005
def unique1d( ar1, retIndx = False ):
    """Unique elements of 1D array. When retIndx is True, return also the
    indices indx such that ar1[indx] is the resulting array of unique
    elements."""
    ar = numpy.array( ar1 ).ravel()
    if retIndx:
        perm = numpy.argsort( ar )
        aux = numpy.take( ar, perm )
        flag = ediff1d( aux, 1 ) != 0
        return numpy.compress( flag, perm ), numpy.compress( flag, aux )
    else:
        aux = numpy.sort( ar )
        return numpy.compress( ediff1d( aux, 1 ) != 0, aux ) 

##
# 01.11.2005, c
def intersect1d( ar1, ar2 ):
    """Intersection of 1D arrays with unique elements."""
    aux = numpy.sort( numpy.concatenate( (ar1, ar2 ) ) )
    return numpy.compress( (aux[1:] - aux[:-1]) == 0, aux )

##
# 01.11.2005, c
def intersect1d_nu( ar1, ar2 ):
    """Intersection of 1D arrays with any elements."""
    # Might be faster then unique1d( intersect1d( ar1, ar2 ) )?
    aux = numpy.sort( numpy.concatenate( (unique1d( ar1 ),
                                          unique1d( ar2  )) ) )
    return numpy.compress( (aux[1:] - aux[:-1]) == 0, aux )

##
# 01.11.2005, c
def setxor1d( ar1, ar2 ):
    """Set exclusive-or of 1D arrays with unique elements."""
    aux = numpy.sort( numpy.concatenate( (ar1, ar2 ) ) )
    flag = ediff1d( aux, toEnd = 1, toBegin = 1 ) == 0
    flag2 = ediff1d( flag, 0 ) == 0
    return numpy.compress( flag2, aux )

##
# 03.11.2005, c
# 05.01.2006
def setmember1d( ar1, ar2 ):
    """Return an array of shape of ar1 containing 1 where the elements of
    ar1 are in ar2 and 0 otherwise."""
    ar = numpy.concatenate( (ar1, ar2 ) )
    tt = numpy.concatenate( (numpy.zeros_like( ar1 ),
                             numpy.zeros_like( ar2 ) + 1) )
    perm = numpy.argsort( ar )
    aux = numpy.take( ar, perm )
    aux2 = numpy.take( tt, perm )
    flag = ediff1d( aux, 1 ) == 0

    ii = numpy.where( flag * aux2 )
    aux = perm[ii+1]
    perm[ii+1] = perm[ii]
    perm[ii] = aux

    indx = numpy.argsort( perm )[:len( ar1 )]

    return numpy.take( flag, indx )

##
# 03.11.2005, c
def union1d( ar1, ar2 ):
    """Union of 1D arrays with unique elements."""
    return unique1d( numpy.concatenate( (ar1, ar2) ) )

##
# 03.11.2005, c
def setdiff1d( ar1, ar2 ):
    """Set difference of 1D arrays with unique elements."""
    aux = setmember1d( ar1, ar2 )
    return numpy.compress( aux == 0, ar1 )

##
# 03.11.2005, c
def test_unique1d():

    a = numpy.array( [5, 7, 1, 2, 1, 5, 7] )

    ec = numpy.array( [1, 2, 5, 7] )
    c = unique1d( a )
    assert numpy.alltrue( c == ec )

##
# 03.11.2005, c
def test_intersect1d():

    a = numpy.array( [5, 7, 1, 2] )
    b = numpy.array( [2, 4, 3, 1, 5] )

    ec = numpy.array( [1, 2, 5] )
    c = intersect1d( a, b )
    assert numpy.alltrue( c == ec )

##
# 03.11.2005, c
def test_intersect1d_nu():

    a = numpy.array( [5, 5, 7, 1, 2] )
    b = numpy.array( [2, 1, 4, 3, 3, 1, 5] )

    ec = numpy.array( [1, 2, 5] )
    c = intersect1d_nu( a, b )
    assert numpy.alltrue( c == ec )

##
# 03.11.2005, c
def test_setxor1d():

    a = numpy.array( [5, 7, 1, 2] )
    b = numpy.array( [2, 4, 3, 1, 5] )

    ec = numpy.array( [3, 4, 7] )
    c = setxor1d( a, b )
    assert numpy.alltrue( c == ec )

    a = numpy.array( [1, 2, 3] )
    b = numpy.array( [6, 5, 4] )

    ec = numpy.array( [1, 2, 3, 4, 5, 6] )
    c = setxor1d( a, b )
    assert numpy.alltrue( c == ec )

    a = numpy.array( [1, 8, 2, 3] )
    b = numpy.array( [6, 5, 4, 8] )

    ec = numpy.array( [1, 2, 3, 4, 5, 6] )
    c = setxor1d( a, b )
    assert numpy.alltrue( c == ec )


##
# 03.11.2005, c
def test_setmember1d():

    a = numpy.array( [5, 7, 1, 2] )
    b = numpy.array( [2, 4, 3, 1, 5] )

    ec = numpy.array( [True, False, True, True] )
    c = setmember1d( a, b )
    assert numpy.alltrue( c == ec )

    a[0] = 8
    ec = numpy.array( [False, False, True, True] )
    c = setmember1d( a, b )
    assert numpy.alltrue( c == ec )

    a[0], a[3] = 4, 8
    ec = numpy.array( [True, False, True, False] )
    c = setmember1d( a, b )
    assert numpy.alltrue( c == ec )

##
# 03.11.2005, c
def test_union1d():

    a = numpy.array( [5, 4, 7, 1, 2] )
    b = numpy.array( [2, 4, 3, 3, 2, 1, 5] )

    ec = numpy.array( [1, 2, 3, 4, 5, 7] )
    c = union1d( a, b )
    assert numpy.alltrue( c == ec )

##
# 03.11.2005, c
def test_setdiff1d():

    a = numpy.array( [6, 5, 4, 7, 1, 2] )
    b = numpy.array( [2, 4, 3, 3, 2, 1, 5] )

    ec = numpy.array( [6, 7] )
    c = setdiff1d( a, b )
    assert numpy.alltrue( c == ec )

##
# 03.11.2005, c
def test_manyways():

    nItem = 100
    a = numpy.fix( nItem / 10 * numpy.random.random( nItem ) )
    b = numpy.fix( nItem / 10 * numpy.random.random( nItem ) )

    c1 = intersect1d_nu( a, b )
    c2 = unique1d( intersect1d( a, b ) )    
    assert numpy.alltrue( c1 == c2 )

    a = numpy.array( [5, 7, 1, 2, 8] )
    b = numpy.array( [9, 8, 2, 4, 3, 1, 5] )

    c1 = setxor1d( a, b )
    aux1 = intersect1d( a, b )
    aux2 = union1d( a, b )
    c2 = setdiff1d( aux2, aux1 )
    assert numpy.alltrue( c1 == c2 )

##
# 02.11.2005, c
def test_unique1d_speed( plotResults = False ):
#    exponents = numpy.linspace( 2, 7, 9 )
    exponents = numpy.linspace( 2, 6, 9 )
    ratios = []
    nItems = []
    dt1s = []
    dt2s = []
    for ii in exponents:

        nItem = 10 ** ii
        print 'using %d items:' % nItem
        a = numpy.fix( nItem / 10 * numpy.random.random( nItem ) )

        print 'dictionary:'
        tt = time.clock() 
        b = numpy.unique( a )
        dt1 = time.clock() - tt
        print dt1

        print 'array:'
        tt = time.clock() 
        c = unique1d( a )
        dt2 = time.clock() - tt
        print dt2


        if dt1 < 1e-8:
            ratio = 'ND'
        else:
            ratio = dt2 / dt1
        print 'ratio:', ratio
        print 'nUnique: %d == %d\n' % (len( b ), len( c ))

        nItems.append( nItem )
        ratios.append( ratio )
        dt1s.append( dt1 )
        dt2s.append( dt2 )

        assert numpy.alltrue( b == c )


    print nItems
    print dt1s
    print dt2s
    print ratios

    if plotResults:
        import pylab

        def plotMe( fig, fun, nItems, dt1s, dt2s ):
            pylab.figure( fig )
            fun( nItems, dt1s, 'g-o', linewidth = 2, markersize = 8 )
            fun( nItems, dt2s, 'b-x', linewidth = 2, markersize = 8 )
            pylab.legend( ('dictionary', 'array' ) )
            pylab.xlabel( 'nItem' )
            pylab.ylabel( 'time [s]' )

        plotMe( 1, pylab.loglog, nItems, dt1s, dt2s )
        plotMe( 2, pylab.plot, nItems, dt1s, dt2s )
        pylab.show()

if (__name__ == '__main__'):

    test_unique1d()
    test_intersect1d()
    test_intersect1d_nu()
    test_setxor1d()
    test_setmember1d()
    test_union1d()
    test_setdiff1d()
    test_manyways()
    test_unique1d_speed( plotResults = True )
