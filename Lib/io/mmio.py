## Automatically adapted for scipy Oct 19, 2005 by convertcode.py

"""
  Matrix Market I/O in Python.
"""
#
# Author: Pearu Peterson <pearu@cens.ioc.ee>
# Created: October, 2004
#
# References:
#  http://math.nist.gov/MatrixMarket/
#
# TODO: support for sparse matrices, need spmatrix.tocoo().

import os
from numpy import asarray, real, imag, conj, zeros, ndarray

__all__ = ['mminfo','mmread','mmwrite']

def mminfo(source):
    """ Queries the contents of the Matrix Market file 'filename' to
    extract size and storage information.

    Inputs:

      source     - Matrix Market filename (extension .mtx) or open file object

    Outputs:

      rows,cols  - number of matrix rows and columns
      entries    - number of non-zero entries of a sparse matrix
                   or rows*cols for a dense matrix
      rep        - 'coordinate' | 'array'
      field      - 'real' | 'complex' | 'pattern' | 'integer'
      symm       - 'general' | 'symmetric' | 'skew-symmetric' | 'hermitian'
    """
    close_it = 0
    if type(source) is type(''):
        if not os.path.isfile(source):
            if source[-4:] != '.mtx':
                source = source + '.mtx'
        source = open(source,'r')
        close_it = 1
    line = source.readline().split()
    if not line[0].startswith('%%MatrixMarket'):
        raise ValueError,'source is not in Matrix Market format'

    assert len(line)==5,`line`

    assert line[1].strip().lower()=='matrix',`line`

    rep = line[2].strip().lower()
    if rep=='dense': rep='array'
    elif rep=='sparse': rep='coordinate'

    field = line[3].strip().lower()

    symm = line[4].strip().lower()

    while line:
        line = source.readline()
        if line.startswith('%'):
            continue
        line = line.split()
        if rep=='array':
            assert len(line)==2,`line`
            rows,cols = map(eval,line)
            entries = rows*cols
        else:
            assert len(line)==3,`line`
            rows,cols,entries = map(eval,line)
        break

    if close_it:
        source.close()
    return (rows,cols,entries,rep,field,symm)

def mmread(source):
    """ Reads the contents of a Matrix Market file 'filename' into a matrix.

    Inputs:

      source    - Matrix Market filename (extensions .mtx, .mtz.gz)
                  or open file object.

    Outputs:

      a         - sparse or full matrix
    """
    close_it = 0
    if type(source) is type(''):
        if not os.path.isfile(source):
            if os.path.isfile(source+'.mtx'):
                source = source + '.mtx'
            elif os.path.isfile(source+'.mtx.gz'):
                source = source + '.mtx.gz'
        if source[-3:] == '.gz':
            import gzip
            source = gzip.open(source)
        else:
            source = open(source,'r')
        close_it = 1

    rows,cols,entries,rep,field,symm = mminfo(source)

    try:
        from scipy.sparse import coo_matrix
    except ImportError:
        coo_matrix = None

    if field=='integer':
        dtype='i'
    elif field=='real':
        dtype='d'
    elif field=='complex':
        dtype='D'
    elif field=='pattern':
        raise NotImplementedError,`field`
    else:
        raise ValueError,`field`

    has_symmetry = symm in ['symmetric','skew-symmetric','hermitian']
    is_complex = field=='complex'
    is_skew = symm=='skew-symmetric'
    is_herm = symm=='hermitian'

    if rep == 'array':
        a = zeros((rows,cols),dtype=dtype)
        line = 1
        i,j = 0,0
        while line:
            line = source.readline()
            if not line or line.startswith('%'):
                continue
            if is_complex:
                aij = complex(*map(float,line.split()))
            else:
                aij = float(line)
            a[i,j] = aij
            if has_symmetry and i!=j:
                if is_skew:
                    a[j,i] = -aij
                elif is_herm:
                    a[j,i] = conj(aij)
                else:
                    a[j,i] = aij
            if i<rows-1:
                i = i + 1
            else:
                j = j + 1
                if not has_symmetry:
                    i = 0
                else:
                    i = j
        assert i in [0,j] and j==cols,`i,j,rows,cols`

    elif rep=='coordinate' and coo_matrix is None:
        # Read sparse matrix to dense when coo_matrix is not available.
        a = zeros((rows,cols), dtype=dtype)
        line = 1
        k = 0
        while line:
            line = source.readline()
            if not line or line.startswith('%'):
                continue
            l = line.split()
            i,j = map(int,l[:2])
            i,j = i-1,j-1
            if is_complex:
                aij = complex(*map(float,l[2:]))
            else:
                aij = float(l[2])
            a[i,j] = aij
            if has_symmetry and i!=j:
                if is_skew:
                    a[j,i] = -aij
                elif is_herm:
                    a[j,i] = conj(aij)
                else:
                    a[j,i] = aij
            k = k + 1
        assert k==entries,`k,entries`

    elif rep=='coordinate':
        k = 0
        data,row,col = [],[],[]
        row_append = row.append
        col_append = col.append
        data_append = data.append
        line = '%'
        while line:
            if not line.startswith('%'):
                l = line.split()
                i = int(l[0])-1
                j = int(l[1])-1
                if is_complex:
                    aij = complex(*map(float,l[2:]))
                else:
                    aij = float(l[2])
                row_append(i)
                col_append(j)
                data_append(aij)
                if has_symmetry and i!=j:
                    if is_skew:
                        aij = -aij
                    elif is_herm:
                        aij = conj(aij)
                    row_append(j)
                    col_append(i)
                    data_append(aij)
                k += 1
            line = source.readline()
        assert k==entries,`k,entries`
        a = coo_matrix((data, (row, col)), dims=(rows, cols), dtype=dtype)
    else:
        raise NotImplementedError,`rep`

    if close_it:
        source.close()
    return a

def mmwrite(target,a,comment='',field=None,precision=None):
    """ Writes the sparse or dense matrix A to a Matrix Market formatted file.

    Inputs:

      target    - Matrix Market filename (extension .mtx) or open file object
      a         - sparse or full matrix
      comment   - comments to be prepended to the Matrix Market file
      field     - 'real' | 'complex' | 'pattern' | 'integer'
      precision - Number of digits to display for real or complex values.
    """
    close_it = 0
    if type(target) is type(''):
        if target[-4:] != '.mtx':
            target = target + '.mtx'
        target = open(target,'w')
        close_it = 1

    if isinstance(a, list) or isinstance(a, ndarray) or isinstance(a, tuple) or hasattr(a,'__array__'):
        rep = 'array'
        a = asarray(a)
        if len(a.shape) != 2:
            raise ValueError, 'expected matrix'
        rows,cols = a.shape
        entries = rows*cols
        typecode = a.dtype.char
        if field is not None:
            if field=='integer':
                a = a.astype('i')
            elif field=='real':
                if typecode not in 'fd':
                    a = a.astype('d')
            elif field=='complex':
                if typecode not in 'FD':
                    a = a.astype('D')
            elif field=='pattern':
                pass
            else:
                raise ValueError,'unknown field '+field
        typecode = a.dtype.char
    else:
        rep = 'coordinate'
        from scipy.sparse import spmatrix
        if not isinstance(a,spmatrix):
            raise ValueError,'unknown matrix type ' + `type(a)`
        rows,cols = a.shape
        entries = a.getnnz()
        typecode = a.gettypecode()

    if precision is None:
        if typecode in 'fF':
            precision = 8
        else:
            precision = 16
    if field is None:
        if typecode in 'li':
            field = 'integer'
        elif typecode in 'df':
            field = 'real'
        elif typecode in 'DF':
            field = 'complex'
        else:
            raise TypeError,'unexpected typecode '+typecode

    if rep == 'array':
        symm = _get_symmetry(a)
    else:
        symm = 'general'

    target.write('%%%%MatrixMarket matrix %s %s %s\n' % (rep,field,symm))

    for line in comment.split('\n'):
        target.write('%%%s\n' % (line))

    if field in ['real','integer']:
        if field=='real':
            format = '%%.%ie\n' % precision
        else:
            format = '%i\n'
    elif field=='complex':
        format = '%%.%ie %%.%ie\n' % (precision,precision)

    if rep == 'array':
        target.write('%i %i\n' % (rows,cols))
        if field in ['real','integer']:
            if symm=='general':
                for j in range(cols):
                    for i in range(rows):
                        target.write(format % a[i,j])
            else:
                for j in range(cols):
                    for i in range(j,rows):
                        target.write(format % a[i,j])
        elif field=='complex':
            if symm=='general':
                for j in range(cols):
                    for i in range(rows):
                        aij = a[i,j]
                        target.write(format % (real(aij),imag(aij)))
            else:
                for j in range(cols):
                    for i in range(j,rows):
                        aij = a[i,j]
                        target.write(format % (real(aij),imag(aij)))
        elif field=='pattern':
            raise ValueError,'Pattern type inconsisted with dense matrix'
        else:
            raise TypeError,'Unknown matrix type '+`field`
    else:
        format = '%i %i ' + format
        target.write('%i %i %i\n' % (rows,cols,entries))
        assert symm=='general',`symm`
        if field in ['real','integer']:
            for i in range(entries):
                target.write(format % (a.rowcol(i)+(a.getdata(i),)))
        elif field=='complex':
            for i in range(entries):
                value = a.getdata(i)
                target.write(format % ((a.rowcol(i))+(real(value),imag(value))))
        elif field=='pattern':
            raise NotImplementedError,`field`
        else:
            raise TypeError,'Unknown matrix type '+`field`

    if close_it:
        target.close()
    else:
        target.flush()
    return

def _get_symmetry(a):
    m,n = a.shape
    if m!=n:
        return 'general'
    issymm = 1
    isskew = 1
    isherm = a.dtype.char in 'FD'
    for j in range(n):
        for i in range(j+1,n):
            aij,aji = a[i][j],a[j][i]
            if issymm and aij != aji:
                issymm = 0
            if isskew and aij != -aji:
                isskew = 0
            if isherm and aij != conj(aji):
                isherm = 0
            if not (issymm or isskew or isherm):
                break
    if issymm: return 'symmetric'
    if isskew: return 'skew-symmetric'
    if isherm: return 'hermitian'
    return 'general'

if __name__ == '__main__':
    import sys
    import time
    for filename in sys.argv[1:]:
        print 'Reading',filename,'...',
        sys.stdout.flush()
        t = time.time()
        mmread(filename)
        print 'took %s seconds' % (time.time() - t)
