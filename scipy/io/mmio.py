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

import os
from numpy import asarray, real, imag, conj, zeros, ndarray, concatenate, \
                  ones, ascontiguousarray, vstack, savetxt, fromfile, fromstring
from numpy.compat import asbytes, asstr

__all__ = ['mminfo','mmread','mmwrite', 'MMFile']


#-------------------------------------------------------------------------------
def mminfo(source):
    """
    Queries the contents of the Matrix Market file 'filename' to
    extract size and storage information.

    Parameters
    ----------

    source : file
        Matrix Market filename (extension .mtx) or open file object

    Returns
    -------

    rows,cols : int
       Number of matrix rows and columns
    entries : int
        Number of non-zero entries of a sparse matrix
        or rows*cols for a dense matrix

    format : {'coordinate', 'array'}

    field : {'real', 'complex', 'pattern', 'integer'}

    symm : {'general', 'symmetric', 'skew-symmetric', 'hermitian'}

    """
    return MMFile.info(source)

#-------------------------------------------------------------------------------
def mmread(source):
    """
    Reads the contents of a Matrix Market file 'filename' into a matrix.

    Parameters
    ----------

    source : file
        Matrix Market filename (extensions .mtx, .mtz.gz)
        or open file object.

    Returns
    -------
    a:
        Sparse or full matrix

    """
    return MMFile().read(source)

#-------------------------------------------------------------------------------
def mmwrite(target, a, comment='', field=None, precision=None):
    """
    Writes the sparse or dense matrix A to a Matrix Market formatted file.

    Parameters
    ----------

    target : file
        Matrix Market filename (extension .mtx) or open file object
    a : array like
        Sparse or full matrix
    comment : str
        comments to be prepended to the Matrix Market file

    field : {'real', 'complex', 'pattern', 'integer'}, optional

    precision :
        Number of digits to display for real or complex values.

    """
    MMFile().write(target, a, comment, field, precision)


################################################################################
class MMFile (object):
    __slots__ = (
      '_rows',
      '_cols',
      '_entries',
      '_format',
      '_field',
      '_symmetry')

    @property
    def rows(self): return self._rows
    @property
    def cols(self): return self._cols
    @property
    def entries(self): return self._entries
    @property
    def format(self): return self._format
    @property
    def field(self): return self._field
    @property
    def symmetry(self): return self._symmetry

    @property
    def has_symmetry(self):
        return self._symmetry in (self.SYMMETRY_SYMMETRIC,
          self.SYMMETRY_SKEW_SYMMETRIC, self.SYMMETRY_HERMITIAN)

    # format values
    FORMAT_COORDINATE = 'coordinate'
    FORMAT_ARRAY = 'array'
    FORMAT_VALUES = (FORMAT_COORDINATE, FORMAT_ARRAY)

    @classmethod
    def _validate_format(self, format):
        if format not in self.FORMAT_VALUES:
            raise ValueError('unknown format type %s, must be one of %s' %
                                (format, self.FORMAT_VALUES))

    # field values
    FIELD_INTEGER = 'integer'
    FIELD_REAL    = 'real'
    FIELD_COMPLEX = 'complex'
    FIELD_PATTERN = 'pattern'
    FIELD_VALUES = (FIELD_INTEGER, FIELD_REAL, FIELD_COMPLEX, FIELD_PATTERN)

    @classmethod
    def _validate_field(self, field):
        if field not in self.FIELD_VALUES:
            raise ValueError('unknown field type %s, must be one of %s' %
                                (field, self.FIELD_VALUES))

    # symmetry values
    SYMMETRY_GENERAL        = 'general'
    SYMMETRY_SYMMETRIC      = 'symmetric'
    SYMMETRY_SKEW_SYMMETRIC = 'skew-symmetric'
    SYMMETRY_HERMITIAN      = 'hermitian'
    SYMMETRY_VALUES = ( SYMMETRY_GENERAL,        SYMMETRY_SYMMETRIC,
                        SYMMETRY_SKEW_SYMMETRIC, SYMMETRY_HERMITIAN)

    @classmethod
    def _validate_symmetry(self, symmetry):
        if symmetry not in self.SYMMETRY_VALUES:
            raise ValueError('unknown symmetry type %s, must be one of %s' %
                                (symmetry, self.SYMMETRY_VALUES))

    DTYPES_BY_FIELD = {
      FIELD_INTEGER: 'i',
      FIELD_REAL:    'd',
      FIELD_COMPLEX: 'D',
      FIELD_PATTERN: 'd'}

    #---------------------------------------------------------------------------
    @staticmethod
    def reader(): pass

    #---------------------------------------------------------------------------
    @staticmethod
    def writer(): pass

    #---------------------------------------------------------------------------
    @classmethod
    def info(self, source):
        source, close_it = self._open(source)

        try:

            # read and validate header line
            line = source.readline()
            mmid, matrix, format, field, symmetry  = \
              [asstr(part.strip().lower()) for part in line.split()]
            if not mmid.startswith('%%matrixmarket'):
                raise ValueError('source is not in Matrix Market format')
            if not matrix == 'matrix':
                raise ValueError("Problem reading file header: " + line)

            # ??? Is this necessary?  I don't see 'dense' or 'sparse' in the spec
            # http://math.nist.gov/MatrixMarket/formats.html
            if format == 'dense':
                format = self.FORMAT_ARRAY
            elif format == 'sparse':
                format = self.FORMAT_COORDINATE

            # skip comments
            while line.startswith(asbytes('%')): line = source.readline()

            line = line.split()
            if format == self.FORMAT_ARRAY:
                if not len(line) == 2:
                    raise ValueError("Header line not of length 2: " + line)
                rows,cols = map(float, line)
                entries = rows*cols
            else:
                if not len(line) == 3:
                    raise ValueError("Header line not of length 3: " + line)
                rows, cols, entries = map(float, line)

            return (rows, cols, entries, format, field, symmetry)

        finally:
            if close_it: source.close()

    #---------------------------------------------------------------------------
    @staticmethod
    def _open(filespec, mode='rb'):
        """
        Return an open file stream for reading based on source.  If source is
        a file name, open it (after trying to find it with mtx and gzipped mtx
        extensions).  Otherwise, just return source.
        """
        close_it = False
        if type(filespec) is type(''):
            close_it = True

            # open for reading
            if mode[0] == 'r':

                # determine filename plus extension
                if not os.path.isfile(filespec):
                    if os.path.isfile(filespec+'.mtx'):
                        filespec = filespec + '.mtx'
                    elif os.path.isfile(filespec+'.mtx.gz'):
                        filespec = filespec + '.mtx.gz'
                    elif os.path.isfile(filespec+'.mtx.bz2'):
                        filespec = filespec + '.mtx.bz2'
                # open filename
                if filespec.endswith('.gz'):
                    import gzip
                    stream = gzip.open(filespec, mode)
                elif filespec.endswith('.bz2'):
                    import bz2
                    stream = bz2.BZ2File(filespec, 'rb')
                else:
                    stream = open(filespec, mode)

            # open for writing
            else:
                if filespec[-4:] != '.mtx':
                    filespec = filespec + '.mtx'
                stream = open(filespec, mode)
        else:
            stream = filespec

        return stream, close_it

    #---------------------------------------------------------------------------
    @staticmethod
    def _get_symmetry(a):
        m,n = a.shape
        if m!=n:
            return MMFile.SYMMETRY_GENERAL
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
        if issymm: return MMFile.SYMMETRY_SYMMETRIC
        if isskew: return MMFile.SYMMETRY_SKEW_SYMMETRIC
        if isherm: return MMFile.SYMMETRY_HERMITIAN
        return MMFile.SYMMETRY_GENERAL

    #---------------------------------------------------------------------------
    @staticmethod
    def _field_template(field, precision):
        return {
          MMFile.FIELD_REAL: '%%.%ie\n' % precision,
          MMFile.FIELD_INTEGER: '%i\n',
          MMFile.FIELD_COMPLEX: '%%.%ie %%.%ie\n' % (precision,precision)
        }.get(field, None)

    #---------------------------------------------------------------------------
    def __init__(self, **kwargs): self._init_attrs(**kwargs)

    #---------------------------------------------------------------------------
    def read(self, source):
        stream, close_it = self._open(source)

        try:
            self._parse_header(stream)
            return self._parse_body(stream)

        finally:
            if close_it: stream.close()

    #---------------------------------------------------------------------------
    def write(self, target, a, comment='', field=None, precision=None):
        stream, close_it = self._open(target, 'wb')

        try:
            self._write(stream, a, comment, field, precision)

        finally:
            if close_it: stream.close()
            else: stream.flush()

    #---------------------------------------------------------------------------
    def _init_attrs(self, **kwargs):
        """
        Initialize each attributes with the corresponding keyword arg value
        or a default of None
        """
        attrs = self.__class__.__slots__
        public_attrs = [attr[1:] for attr in attrs]
        invalid_keys = set(kwargs.keys()) - set(public_attrs)

        if invalid_keys:
            raise ValueError('found %s invalid keyword arguments, please only use %s' %
                                (tuple(invalid_keys), public_attrs))

        for attr in attrs: setattr(self, attr, kwargs.get(attr[1:], None))

    #---------------------------------------------------------------------------
    def _parse_header(self, stream):
        rows, cols, entries, format, field, symmetry = \
          self.__class__.info(stream)
        self._init_attrs(rows=rows, cols=cols, entries=entries, format=format,
          field=field, symmetry=symmetry)

    #---------------------------------------------------------------------------
    def _parse_body(self, stream):
        rows, cols, entries, format, field, symm = \
          (self.rows, self.cols, self.entries, self.format, self.field, self.symmetry)

        try:
            from scipy.sparse import coo_matrix
        except ImportError:
            coo_matrix = None

        dtype = self.DTYPES_BY_FIELD.get(field, None)

        has_symmetry = self.has_symmetry
        is_complex = field == self.FIELD_COMPLEX
        is_skew = symm == self.SYMMETRY_SKEW_SYMMETRIC
        is_herm = symm == self.SYMMETRY_HERMITIAN
        is_pattern = field == self.FIELD_PATTERN

        if format == self.FORMAT_ARRAY:
            a = zeros((rows,cols),dtype=dtype)
            line = 1
            i,j = 0,0
            while line:
                line = stream.readline()
                if not line or line.startswith(asbytes('%')):
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
            if not (i in [0,j] and j == cols):
                raise ValueError("Parse error, did not read all lines.")

        elif format == self.FORMAT_COORDINATE and coo_matrix is None:
            # Read sparse matrix to dense when coo_matrix is not available.
            a = zeros((rows,cols), dtype=dtype)
            line = 1
            k = 0
            while line:
                line = stream.readline()
                if not line or line.startswith(asbytes('%')):
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
            if not k == entries:
                ValueError("Did not read all entries")

        elif format == self.FORMAT_COORDINATE:
            # Read sparse COOrdinate format

            if entries == 0:
                # empty matrix
                return coo_matrix((rows, cols), dtype=dtype)

            try:
                # fromfile works for normal files
                flat_data = fromfile(stream, sep=' ')
            except:
                # fallback - fromfile fails for some file-like objects
                flat_data = fromstring(stream.read(), sep=' ')

                # TODO use iterator (e.g. xreadlines) to avoid reading
                # the whole file into memory

            if is_pattern:
                flat_data = flat_data.reshape(-1,2)
                I = ascontiguousarray(flat_data[:,0], dtype='intc')
                J = ascontiguousarray(flat_data[:,1], dtype='intc')
                V = ones(len(I), dtype='int8')  # filler
            elif is_complex:
                flat_data = flat_data.reshape(-1,4)
                I = ascontiguousarray(flat_data[:,0], dtype='intc')
                J = ascontiguousarray(flat_data[:,1], dtype='intc')
                V = ascontiguousarray(flat_data[:,2], dtype='complex')
                V.imag = flat_data[:,3]
            else:
                flat_data = flat_data.reshape(-1,3)
                I = ascontiguousarray(flat_data[:,0], dtype='intc')
                J = ascontiguousarray(flat_data[:,1], dtype='intc')
                V = ascontiguousarray(flat_data[:,2], dtype='float')

            I -= 1 #adjust indices (base 1 -> base 0)
            J -= 1

            if has_symmetry:
                mask = (I != J)       #off diagonal mask
                od_I = I[mask]
                od_J = J[mask]
                od_V = V[mask]

                I = concatenate((I,od_J))
                J = concatenate((J,od_I))

                if is_skew:
                    od_V *= -1
                elif is_herm:
                    od_V = od_V.conjugate()

                V = concatenate((V,od_V))

            a = coo_matrix((V, (I, J)), shape=(rows, cols), dtype=dtype)
        else:
            raise NotImplementedError(format)

        return a

    #---------------------------------------------------------------------------
    def _write(self, stream, a, comment='', field=None, precision=None):

        if isinstance(a, list) or isinstance(a, ndarray) or isinstance(a, tuple) or hasattr(a,'__array__'):
            rep = self.FORMAT_ARRAY
            a = asarray(a)
            if len(a.shape) != 2:
                raise ValueError('expected matrix')
            rows,cols = a.shape
            entries = rows*cols

            if field is not None:

                if field == self.FIELD_INTEGER:
                    a = a.astype('i')
                elif field == self.FIELD_REAL:
                    if a.dtype.char not in 'fd':
                        a = a.astype('d')
                elif field == self.FIELD_COMPLEX:
                    if a.dtype.char not in 'FD':
                        a = a.astype('D')

        else:
            from scipy.sparse import spmatrix
            if not isinstance(a,spmatrix):
                raise ValueError('unknown matrix type: %s' % type(a))
            rep = 'coordinate'
            rows, cols = a.shape
            entries = a.getnnz()

        typecode = a.dtype.char

        if precision is None:
            if typecode in 'fF':
                precision = 8
            else:
                precision = 16

        if field is None:
            kind = a.dtype.kind
            if kind == 'i':
                field = 'integer'
            elif kind == 'f':
                field = 'real'
            elif kind == 'c':
                field = 'complex'
            else:
                raise TypeError('unexpected dtype kind ' + kind)

        if rep == self.FORMAT_ARRAY:
            symm = self._get_symmetry(a)
        else:
            symm = self.SYMMETRY_GENERAL

        # validate rep, field, and symmetry
        self.__class__._validate_format(rep)
        self.__class__._validate_field(field)
        self.__class__._validate_symmetry(symm)

        # write initial header line
        stream.write(asbytes('%%%%MatrixMarket matrix %s %s %s\n' % (rep,field,symm)))

        # write comments
        for line in comment.split('\n'):
            stream.write(asbytes('%%%s\n' % (line)))


        template = self._field_template(field, precision)

        # write dense format
        if rep == self.FORMAT_ARRAY:

            # write shape spec
            stream.write(asbytes('%i %i\n' % (rows,cols)))

            if field in (self.FIELD_INTEGER, self.FIELD_REAL):

                if symm == self.SYMMETRY_GENERAL:
                    for j in range(cols):
                        for i in range(rows):
                            stream.write(asbytes(template % a[i,j]))
                else:
                    for j in range(cols):
                        for i in range(j,rows):
                            stream.write(asbytes(template % a[i,j]))

            elif field == self.FIELD_COMPLEX:

                if symm == self.SYMMETRY_GENERAL:
                    for j in range(cols):
                        for i in range(rows):
                            aij = a[i,j]
                            stream.write(asbytes(template % (real(aij),imag(aij))))
                else:
                    for j in range(cols):
                        for i in range(j,rows):
                            aij = a[i,j]
                            stream.write(asbytes(template % (real(aij),imag(aij))))

            elif field == self.FIELD_PATTERN:
                raise ValueError('pattern type inconsisted with dense format')

            else:
                raise TypeError('Unknown field type %s' % field)

        # write sparse format
        else:

            if symm != self.SYMMETRY_GENERAL:
                raise NotImplementedError('symmetric matrices not yet supported')

            coo = a.tocoo() # convert to COOrdinate format

            # write shape spec
            stream.write(asbytes('%i %i %i\n' % (rows, cols, coo.nnz)))

            fmt = '%%.%dg' % precision

            if field == self.FIELD_PATTERN:
                IJV = vstack((coo.row, coo.col)).T
            elif field in [ self.FIELD_INTEGER, self.FIELD_REAL ]:
                IJV = vstack((coo.row, coo.col, coo.data)).T
            elif field == self.FIELD_COMPLEX:
                IJV = vstack((coo.row, coo.col, coo.data.real, coo.data.imag)).T
            else:
                raise TypeError('Unknown field type %s' % field)

            IJV[:,:2] += 1 # change base 0 -> base 1

            savetxt(stream, IJV, fmt=fmt)

#-------------------------------------------------------------------------------
if __name__ == '__main__':
    import sys
    import time
    for filename in sys.argv[1:]:
        print 'Reading',filename,'...',
        sys.stdout.flush()
        t = time.time()
        mmread(filename)
        print 'took %s seconds' % (time.time() - t)
